"""
AIMNet2 plugin for Q|R.
Based on the AIMNet2 ASE calculator by the authors of the AIMNet2 model.
"""

from __future__ import print_function
from __future__ import absolute_import
import numpy as np
import ase.units as ase_units

import torch
from torch import nn, Tensor
from typing import Union, Dict, Any
import os
import requests
from ase.calculators.calculator import Calculator, all_changes
import warnings

import numba
import numpy as np
try:
    import numba.cuda
    assert numba.cuda.is_available()
    _numba_cuda = True
except:
    _numba_cuda = False


# model registry aliases
model_registry_aliases = {}
model_registry_aliases['aimnet2-qr'] = 'aimnet2-qr/aimnet2-qr_b97md4_qzvp_2'


def get_model_path(s: str):
    # direct file path
    if os.path.isfile(s):
        print('Found model file:', s)
        return s
    # check aliases
    if s in model_registry_aliases:
        s = model_registry_aliases[s]
    # add jpt extension
    if not s.endswith('.jpt'):
        s = s + '.jpt'
    sdir = os.path.dirname(s)
    os.makedirs(os.path.join(os.path.dirname(__file__), 'assets', sdir), exist_ok=True)
    s_local = os.path.join(os.path.dirname(__file__), 'assets', s)
    if os.path.isfile(s_local):
        print('Found model file:', s_local)
    else:
        url = f'https://github.com/zubatyuk/aimnet-model-zoo/raw/main/{s}'
        print('Downloading model file from', url)
        r = requests.get(url)
        r.raise_for_status()
        with open(s_local, 'wb') as f:
            f.write(r.content)
        print('Saved to ', s_local)
    return s_local


def calc_nbmat(coord, cutoff: float, max_nb: int = 128):
    """ Calculate neighbor list matrix for a given set of coordinates.
    Args:
    coord (torch.Tensor): Tensor of shape (N, 3) with atomic coordinates.
    cutoff (float): Cutoff radius for the neighbor list.
    max_nb (int): Maximum number of neighbors to consider.
    Returns:
    torch.Tensor: Neighbor list matrix of shape (N, max_nb) with neighbor indices.
    """
    device = coord.device
    coord = coord.detach()
    if coord.is_cuda and not _numba_cuda:
        warnings.warn('Numba CUDA is not available, falling back to CPU.')
        coord = coord.cpu()
    if coord.is_cuda:
        N = coord.shape[0]
        nbmat = torch.full((N+1, max_nb), N, dtype=torch.int32, device=device)
        nnb = torch.zeros(N, dtype=torch.int32, device='cuda')
        n_coord = numba.cuda.as_cuda_array(coord)
        n_nbmat = numba.cuda.as_cuda_array(nbmat)
        n_nnb = numba.cuda.as_cuda_array(nnb)
        threadsperblock = 32
        blockspergrid = (N + (threadsperblock - 1)) // threadsperblock
        _nbmat_kernel_cuda[blockspergrid, threadsperblock](n_coord, n_nbmat, n_nnb, cutoff**2, max_nb)
        return nbmat.to(torch.long), nnb.to(torch.long)
    else:
        coord = coord.cpu().numpy()
        nbmat, nnb = _nbmat_kernel_cpu(coord, cutoff, max_nb)
        return torch.as_tensor(nbmat, dtype=torch.long, device=device), torch.as_tensor(nnb, dtype=torch.long, device=device)


if _numba_cuda:
    @numba.cuda.jit
    def _nbmat_kernel_cuda(coord, nbmat, nnb, cufoff_squared, maxnb):
        N = coord.shape[0]
        i = numba.cuda.grid(1)

        c0 = coord[i, 0]
        c1 = coord[i, 1]
        c2 = coord[i, 2]

        for j in range(i+1, N):
            d0 = c0 - coord[j, 0]
            d1 = c1 - coord[j, 1]
            d2 = c2 - coord[j, 2]
            dist_squared = d0 * d0 + d1 * d1 + d2 * d2
            if dist_squared > cufoff_squared:
                continue

            pos = numba.cuda.atomic.add(nnb, i, 1)
            if pos < maxnb:
                nbmat[i, pos] = j
            pos = numba.cuda.atomic.add(nnb, j, 1)
            if pos < maxnb:
                nbmat[j, pos] = i


@numba.njit(parallel=True)
def _nbmat_kernel_cpu(coord, cutoff, maxnb):
    N = coord.shape[0]
    nbmat = np.full((N+1, maxnb), N, dtype=np.int32)
    nnb = np.zeros(N, dtype=np.int32)
    cutoff2 = cutoff * cutoff
    for i in numba.prange(N):
        c_i = coord[i]
        pos = 0
        for j in range(i+1, N):
            c_j = coord[j]
            diff = c_i - c_j
            dist2 = (diff * diff).sum(-1)
            if dist2 < cutoff2:
                nbmat[i, pos] = j
                pos += 1
                if pos < maxnb:
                    nnb[i] += 1
    nnb_half = nnb.copy()
    for i in range(N):
        for m in range(nnb_half[i]):
            j = nbmat[i, m]
            pos = nnb[j]
            nbmat[j, pos] = i
            nnb[j] += 1

    return nbmat, nnb


class BaseAIMNet2Calculator:
    """ Base AIMNet2 calculator 
    A helper class to load AIMNet2 models and perform inference.
    """

    keys_in = {
        'coord': torch.float,
        'numbers': torch.int,
        'charge': torch.float
        }
    keys_out = ['energy', 'charges', 'forces']
    atom_feature_keys = ['coord', 'numbers', 'charges', 'forces']
    
    def __init__(self, model: Union[str, torch.nn.Module] = 'aimnet2'):
        if torch.cuda.is_available():
            self.device = torch.device('cuda')
        elif torch.backends.mps.is_available():
            self.device = torch.device('mps')
        else:
            self.device = torch.device('cpu')

        _x = torch.zeros(1, device=self.device)
        print('Running AIMNet2/PyTorch on device:', str(_x.device).upper())
        _numba_device = 'CUDA' if _x.is_cuda else 'CPU'
        print('Running Numba on device:', _numba_device)

        if isinstance(model, str):
            p = get_model_path(model)
            self.model = torch.jit.load(p, map_location=self.device)
        elif isinstance(model, nn.Module):
            self.model = model.to(self.device)
        else:
            raise AttributeError('Invalid model type/name.')
        # indicator if input was flattened
        self._batch = None
        # placeholder for tensors that require grad
        self._saved_for_grad = None        
        self.cutoff = self.model.cutoff

    def __call__(self, *args, **kwargs):
        return self.eval(*args, **kwargs)

    def eval(self, data: Dict[str, Any], forces=False) -> Dict[str, Tensor]:
        data = self.prepare_input(data)
        data = self.set_grad_tensors(data, forces=forces)
        with torch.jit.optimized_execution(False):
            data = self.model(data)
        data = self.get_derivatives(data, forces=forces)
        data = self.process_output(data)
        return data
        
    def prepare_input(self, data: Dict[str, Any]) -> Dict[str, Tensor]:
        data = self.to_input_tensors(data)
        data = self.mol_flatten(data)
        data = self.make_nbmat(data)
        data = self.pad_input(data)
        return data
    
    def process_output(self, data: Dict[str, Tensor]) -> Dict[str, Tensor]:
        data = self.unpad_output(data)
        data = self.mol_unflatten(data)
        data = self.keep_only(data)
        return data

    def to_input_tensors(self, data: Dict[str, Any]) -> Dict[str, Tensor]:
        ret = dict()
        for k in self.keys_in:
            assert k in data, f'Missing key {k} in the input data'
            # always detach !!
            ret[k] = torch.as_tensor(data[k], device=self.device, dtype=self.keys_in[k]).detach()
        # convert any scalar tensors to shape (1,) tensors
        for k, v in ret.items():
            if v.ndim == 0:
                ret[k] = v.unsqueeze(0)
        return ret

    def mol_flatten(self, data: Dict[str, Tensor]) -> Dict[str, Tensor]:
        assert data['coord'].ndim in {2, 3}, 'Expected 2D or 3D tensor for coord'
        if data['coord'].ndim == 3:
            B, N = data['coord'].shape[:2]
            self._batch = B
            data['mol_idx'] = torch.repeat_interleave(torch.arange(0, B, device=self.device), torch.full((B,), N, device=self.device))
            for k, v in data.items():
                if k in self.atom_feature_keys:
                    assert v.ndim >= 2, f'Expected at least 2D tensor for {k}, got {v.ndim}D'
                    data[k] = v.flatten(0, 1)
        else:
            self._batch = None
            if 'mol_idx' not in data:
                data['mol_idx'] = torch.zeros(data['coord'].shape[0], device=self.device)
        return data
    
    def mol_unflatten(self, data: Dict[str, Tensor], batch=None) -> Dict[str, Tensor]:
        batch = batch or self._batch
        if batch is not None:
            for k, v in data.items():
                if k in self.atom_feature_keys:
                    data[k] = v.view(self._batch, -1, *v.shape[1:])
        return data
    
    def make_nbmat(self, data: Dict[str, Tensor]) -> Dict[str, Tensor]:        
        if 'nbmat' not in data:
            maxnb = 128
            while True:
                nbmat, nnb = calc_nbmat(data['coord'], self.cutoff, maxnb)
                if (nnb > maxnb-2).any():
                    maxnb *= 2
                else:
                    break
            data['nbmat'] = nbmat
        return data
    
    def pad_input(self, data: Dict[str, Tensor]) -> Dict[str, Tensor]:
        N = data['nbmat'].shape[0]
        data['coord'] = maybe_pad_dim0(data['coord'], N)
        data['numbers'] = maybe_pad_dim0(data['numbers'], N)
        data['mol_idx'] = maybe_pad_dim0(data['mol_idx'], N, value=data['mol_idx'][-1])
        return data

    def unpad_output(self, data: Dict[str, Tensor]) -> Dict[str, Tensor]:
        N = data['nbmat'].shape[0] - 1
        for k, v in data.items():
            if k in self.atom_feature_keys:
                data[k] = maybe_unpad_dim0(v, N)
        return data
    
    def set_grad_tensors(self, data: Dict[str, Tensor], forces=False) -> Dict[str, Tensor]:
        self._saved_for_grad = dict() 
        if forces:
            data['coord'].requires_grad_(True)
            self._saved_for_grad['coord'] = data['coord']
        return data
    
    def keep_only(self, data: Dict[str, Tensor]) -> Dict[str, Tensor]:
        ret = dict()
        for k, v in data.items():
            if k in self.keys_out or (k.endswith('_std') and k[:-4] in self.keys_out):
                ret[k] = v
        return ret
    
    def get_derivatives(self, data: Dict[str, Tensor], forces=False) -> Dict[str, Tensor]:
        if forces:
            tot_energy = data['energy'].sum()
            deriv = torch.autograd.grad(tot_energy, self._saved_for_grad['coord'])
            data['forces'] = - deriv[0]
        return data

    
def maybe_pad_dim0(a: Tensor, N: int, value=0.0) -> Tensor:
    _shape_diff = N - a.shape[0]
    assert _shape_diff == 0 or _shape_diff == 1, 'Invalid shape'
    if _shape_diff == 1:
        a = pad_dim0(a, value=value)
    return a


def pad_dim0(a: Tensor, value=0.0) -> Tensor:
    shapes = [0] * ((a.ndim - 1)*2) + [0, 1]
    a = torch.nn.functional.pad(a, shapes, mode='constant', value=value)
    return a


def maybe_unpad_dim0(a: Tensor, N: int) -> Tensor:
    _shape_diff = a.shape[0] - N
    assert _shape_diff == 0 or _shape_diff == 1, 'Invalid shape'
    if _shape_diff == 1:
        a = a[:-1]
    return a


class AIMNet2ASE(Calculator):
    implemented_properties = ['energy', 'forces', 'free_energy', 'charges']
    def __init__(self, base_calc: Union[BaseAIMNet2Calculator, str] = 'aimnet2', charge=0):
        super().__init__()
        if isinstance(base_calc, str):
            base_calc = BaseAIMNet2Calculator(base_calc)
        self.base_calc = base_calc
        self.charge = charge
        self.do_reset()
        # list of implemented species
        if hasattr(base_calc, 'implemented_species'):
            self.implemented_species = base_calc.implemented_species.cpu().numpy()
        else:
            self.implemented_species = None
        

    def do_reset(self):
        self._t_numbers = None
        self._t_charge = None
        self._t_mol_idx = None
        self.charge = 0.0

    def set_atoms(self, atoms):
        if self.implemented_species is not None and not np.in1d(atoms.numbers, self.implemented_species).all():
            raise ValueError('Some species are not implemented in the AIMNet2Calculator')
        self.atoms = atoms
        self.do_reset()

    def set_charge(self, charge):
        self.charge = charge

    def uptade_tensors(self):
        if self._t_numbers is None:
            self._t_numbers = torch.tensor(self.atoms.numbers, dtype=torch.int64, device=self.base_calc.device)
        if self._t_charge is None:
            self._t_charge = torch.tensor(self.charge, dtype=torch.float32, device=self.base_calc.device)
        if self._t_mol_idx is None:
            self.mol_idx = torch.zeros(len(self.atoms), dtype=torch.int64, device=self.base_calc.device)            

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.uptade_tensors()

        results = self.base_calc({
            'coord': torch.tensor(self.atoms.positions, dtype=torch.float32, device=self.base_calc.device),
            'numbers': self._t_numbers,
            'charge': self._t_charge,
        }, forces='forces' in properties)
        for k, v in results.items():
            results[k] = v.detach().cpu().numpy()

        self.results['energy'] = results['energy']
        self.results['charges'] = results['charges']
        if 'forces' in properties:
            self.results['forces'] = results['forces']


class AIMNet2Calculator(AIMNet2ASE):
    """ Modification of the AIMNet2ASE class to work with Q|R.
    """
    def run_qr(self, atoms, coordinates, charge, pointcharges, define_str=None):
        """
        This method is called every time an energy and forces are needed.
        The Q|R code calls this method at each step of LBFGS.
        Args:
        atoms (ase.atoms.Atoms) an updated set of atoms.
        TODO: coordinates (numpy.ndarray) an updated set of coordinates. are these even being used?
        charge (int) the charge on the molecular system.
        TODO: pointcharges (?) are these even being used?
        command (str) not used in this calculator
        define_str() not used in this calculator, legacy from Turbomole?
        """
        atoms.calc = self
        self.atoms = atoms
        self.set_charge(charge)
        unit_convert = ase_units.kcal / ase_units.mol
        self.calculate(atoms, properties=['energy', 'forces'])
        self.energy_free = self.results['energy'][0] * unit_convert
        self.forces = self.results['forces'].astype(np.float64) * unit_convert

        if 0: # we need debugging flag here to switch on and off.
            # print(("AIMNet2: ",self.method))
            print(('Energy:',self.energy_free))
            print(('Force:', self.forces))
