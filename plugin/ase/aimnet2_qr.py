"""
The AIMNet2 calculator can compute an energy and force for a molecular system using a Artificial Neural Network.
"""

from __future__ import print_function
from __future__ import absolute_import
import os
import numpy as np
import ase.units as ase_units
from ase.calculators.calculator import Calculator, all_changes
import torch
import torch.nn.functional as F
import numba
from numba import cuda
import libtbx.load_env
from libtbx import easy_run
from scitbx.array_family import flex
from libtbx import easy_mp
import traceback
from libtbx.utils import Sorry

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_aimnet_models = os.path.join(qrefine, "plugin","ase","ani")
MODEL_FILE=os.path.join(qr_aimnet_models, "aimnet2_qr_pbeh_cpcm_230419.jpt")
#MODEL_FILE=os.path.join(qr_aimnet_models, "aimnet2_qr_b97m_cpcm_230419.jpt")
DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
MODEL = torch.jit.load(MODEL_FILE).to(DEVICE)


### dense neighbor matrix kernels
@numba.njit(cache=True, parallel=True)
def _cpu_dense_nb_mat_sft(conn_matrix):
    N, S = conn_matrix.shape[:2]
    # figure out max number of neighbors
    _s_flat_conn_matrix = conn_matrix.reshape(N, -1)
    maxnb = np.max(np.sum(_s_flat_conn_matrix, axis=-1))
    M = maxnb
    # atom idx matrix
    mat_idxj = np.full((N + 1, M), -1, dtype=np.int_)
    # padding matrix
    mat_pad = np.ones((N + 1, M), dtype=np.bool_)
    # shitfs matrix
    mat_S_idx = np.zeros((N + 1, M), dtype=np.int_)
    for _n in numba.prange(N):
        _i = 0
        for _s in range(S):
            for _m in range(N):
                if conn_matrix[_n, _s, _m] == True:
                    mat_idxj[_n, _i] = _m
                    mat_pad[_n, _i] = False
                    mat_S_idx[_n, _i] = _s
                    _i += 1
    return mat_idxj, mat_pad, mat_S_idx


@numba.njit(cache=True, parallel=True)
def _cpu_dense_nb_mat(conn_matrix):
    N = conn_matrix.shape[0]
    # figure out max number of neighbors
    _s_flat_conn_matrix = conn_matrix.reshape(N, -1)
    maxnb = np.max(np.sum(_s_flat_conn_matrix, axis=-1))
    M = maxnb
    # atom idx matrix
    mat_idxj = np.full((N + 1, M), -1, dtype=np.int_)
    # padding matrix
    mat_pad = np.ones((N + 1, M), dtype=np.bool_)
    for _n in numba.prange(N):
        _i = 0
        for _m in range(N):
            if conn_matrix[_n, _m] == True:
                mat_idxj[_n, _i] = _m
                mat_pad[_n, _i] = False
                _i += 1
    return mat_idxj, mat_pad


@cuda.jit(cache=True)
def _cuda_dense_nb_mat_sft(conn_matrix, mat_idxj, mat_pad, mat_S_idx):
    i = cuda.grid(1)
    if i < conn_matrix.shape[0]:
        k = 0
        for s in range(conn_matrix.shape[1]):
            for j in range(conn_matrix.shape[2]):
                if conn_matrix[i, s, j] > 0:
                    mat_idxj[i, k] = j
                    mat_pad[i, k] = 0
                    mat_S_idx[i, k] = s
                    k += 1


@cuda.jit(cache=True)
def _cuda_dense_nb_mat(conn_matrix, mat_idxj, mat_pad):
    i = cuda.grid(1)
    if i < conn_matrix.shape[0]:
        k = 0
        for j in range(conn_matrix.shape[1]):
            if conn_matrix[i, j] > 0:
                mat_idxj[i, k] = j
                mat_pad[i, k] = 0
                k += 1


class AIMNet2Calculator(Calculator):
    """ ASE calculator for AIMNet2 model
    Arguments:
        model (:class:`torch.nn.Module`): AIMNet2 model
        cutoff (float): Short-range cutoff. Usually about 5 Angs.
            Should not be smaller then specified within the model. Default: guess from the model from `aev.rc_s` parameter.
        use_coulomb (bool): provide long-range neighbour lists to the model. Default: False
        coulomb_cutoff (float): Long-range cutoff. Default: 15 Angs
    """

    implemented_properties = ['energy', 'forces', 'free_energy', 'charges']

    def __init__(self, cutoff=None, use_coulomb=False, use_pbc=False, coulomb_cutoff=15.0):
        super().__init__()
        self.device = DEVICE  
        self.model = MODEL
        if cutoff is None:
            cutoff = max(v.item() for k, v in self.model.state_dict().items() if k.endswith('aev.rc_s'))
        self.cutoff = float(cutoff)
        self.use_coulomb = use_coulomb
        self.coulomb_cutoff = float(coulomb_cutoff)
        self.species = np.array([1, 6, 7, 8, 16, 34])
        self.method = 'qr-ef'
        self.use_pbc = use_pbc 
        self.do_reset()

    def do_reset(self):
        self._nblist = dict()
        self._t_numbers = None
        self._t_charge = None
        if self.use_coulomb:
            self._nblist['coul_cutoff'] = torch.tensor(self.coulomb_cutoff, dtype=torch.float, device=self.device)
        self.charge = 0.0

    def set_charge(self, charge):
        self.charge = float(charge)

    def set_species(self, species):
        self.species = np.array([int(x) for x in species])

    def _make_input(self):
        coord = torch.as_tensor(self.atoms.positions).to(torch.float).to(self.device)
        if self._t_numbers is None:
            self._t_numbers = torch.as_tensor(self.atoms.numbers).to(torch.long).to(self.device)
            self._t_numbers = F.pad(self._t_numbers, (0, 1))
            self._t_charge = torch.tensor([self.charge], dtype=torch.float, device=self.device)
        d = dict(coord=coord, numbers=self._t_numbers, charge=self._t_charge)
        if self.use_pbc and self.atoms.pbc.any():
            cell = torch.as_tensor(self.atoms.cell.array, dtype=torch.float32, device=self.device)
            d['cell'] = cell
            dn, d['coord'] = self.calc_nblist_pbc(coord, cell)
            d.update(dn)
        else:
            d.update(self.calc_nblist(coord))
        d['coord'] = F.pad(d['coord'], (0, 0, 0, 1))
        return d

    def calc_nblist(self, coord):
        coord = coord.detach()
        dmat = torch.cdist(coord, coord)
        conn_mat = dmat < self.cutoff
        conn_mat.fill_diagonal_(False)
        if self.device.type == 'cuda':
            self._nblist['idx_j'], self._nblist['nb_pad_mask'] = self._nblist_cuda(conn_mat)
        else:
            self._nblist['idx_j'], self._nblist['nb_pad_mask'] = self._nblist_cpu(conn_mat)
        if self.use_coulomb:
            conn_mat = dmat < self.coulomb_cutoff
            conn_mat.fill_diagonal_(False)
            if self.device.type == 'cuda':
                self._nblist['idx_j_coul'], self._nblist['nb_pad_mask_coul'] = self._nblist_cuda(conn_mat)
            else:
                self._nblist['idx_j_coul'], self._nblist['nb_pad_mask_coul'] = self._nblist_cpu(conn_mat)
        return self._nblist

    def _nblist_cpu(self, conn_mat):
        conn_mat = conn_mat.cpu().numpy()
        idx_j, mat_pad = _cpu_dense_nb_mat(conn_mat)
        idx_j = torch.as_tensor(idx_j).to(self.device)
        mat_pad = torch.as_tensor(mat_pad).to(self.device)
        return idx_j, mat_pad

    def _nblist_cuda(self, conn_mat):
        N = conn_mat.shape[0]
        M = conn_mat.sum(-1).max()
        threadsperblock = 8
        blockspergrid = (N + (threadsperblock - 1)) // threadsperblock
        idx_j = torch.full((N + 1, M), -1, dtype=torch.int64, device=self.device)
        mat_pad = torch.ones((N + 1, M), dtype=torch.int8, device=self.device)
        conn_mat = conn_mat.to(torch.int8)
        _conn_mat = cuda.as_cuda_array(conn_mat)
        _idx_j = cuda.as_cuda_array(idx_j)
        _mat_pad = cuda.as_cuda_array(mat_pad)
        _cuda_dense_nb_mat[blockspergrid, threadsperblock](_conn_mat, _idx_j, _mat_pad)
        mat_pad = mat_pad.to(torch.bool)
        return idx_j, mat_pad

    def calc_nblist_pbc(self, coord, cell):
        coord, cell = coord.detach(), cell.detach()
        reciprocal_cell = cell.inverse().t()
        coord = (coord @ cell.inverse()).fmod(1) @ cell
        inv_distances = reciprocal_cell.norm(2, -1)
        shifts = self._calc_shifts(inv_distances, self.cutoff)
        d = torch.cdist(coord.unsqueeze(0), coord.unsqueeze(0) + (shifts @ cell).unsqueeze(1))
        conn_mat = ((d < self.cutoff) & (d > 0.1)).transpose(0, 1).contiguous()
        if self.device.type == 'cuda':
            mat_idxj, mat_pad, mat_S = self._nblist_pbc_cuda(conn_mat, shifts)
        else:
            mat_idxj, mat_pad, mat_S = self._nblist_pbc_cpu(conn_mat, shifts)
        self._nblist['idx_j'], self._nblist['nb_pad_mask'], self._nblist['shifts'] = mat_idxj, mat_pad, mat_S
        if self.use_coulomb:
            shifts = self._calc_shifts(inv_distances, self.coulomb_cutoff)
            d = torch.cdist(coord.unsqueeze(0), coord.unsqueeze(0) + (shifts @ cell).unsqueeze(1))
            conn_mat = ((d < self.cutoff) & (d > 0.1)).transpose(0, 1).contiguous()
            if self.device.type == 'cuda':
                mat_idxj, mat_pad, mat_S = self._nblist_pbc_cuda(conn_mat, shifts)
            else:
                mat_idxj, mat_pad, mat_S = self._nblist_pbc_cpu(conn_mat, shifts)
            self._nblist['idx_j_coul'], self._nblist['nb_pad_mask_coul'], self._nblist['shifts_coul'] = mat_idxj, mat_pad, mat_S
        return self._nblist, coord

    def _calc_shifts(self, inv_distances, cutoff):
        num_repeats = torch.ceil(cutoff * inv_distances).to(torch.long)
        dc = [torch.arange(-num_repeats[i], num_repeats[i] + 1, device=self.device) for i in range(len(num_repeats))]
        shifts = torch.cartesian_prod(*dc).to(torch.float)
        return shifts

    def _nblist_pbc_cpu(self, conn_mat, shifts):
        conn_mat = conn_mat.cpu().numpy()
        mat_idxj, mat_pad, mat_S_idx = _cpu_dense_nb_mat_sft(conn_mat)
        mat_idxj = torch.from_numpy(mat_idxj).to(self.device)
        mat_pad = torch.from_numpy(mat_pad).to(self.device)
        mat_S_idx = torch.from_numpy(mat_S_idx).to(self.device)
        mat_S = shifts[mat_S_idx]
        return mat_idxj, mat_pad, mat_S

    def _nblist_pbc_cuda(self, conn_mat, shifts):
        N = conn_mat.shape[0]
        M = conn_mat.view(N, -1).sum(-1).max()
        threadsperblock = 32
        blockspergrid = (N + (threadsperblock - 1)) // threadsperblock
        idx_j = torch.full((N + 1, M), -1, dtype=torch.int64, device=self.device)
        mat_pad = torch.ones((N + 1, M), dtype=torch.int8, device=self.device)
        S_idx = torch.zeros((N + 1, M), dtype=torch.int64, device=self.device)
        conn_mat = conn_mat.to(torch.int8)
        _conn_mat = cuda.as_cuda_array(conn_mat)
        _idx_j = cuda.as_cuda_array(idx_j)
        _mat_pad = cuda.as_cuda_array(mat_pad)
        _S_idx = cuda.as_cuda_array(S_idx)
        _cuda_dense_nb_mat_sft[blockspergrid, threadsperblock](_conn_mat, _idx_j, _mat_pad, _S_idx)
        mat_pad = mat_pad.to(torch.bool)
        return idx_j, mat_pad, shifts[S_idx]

    def _eval_model(self, d, forces=True):
        prev = torch.is_grad_enabled()
        torch._C._set_grad_enabled(forces)
        if forces:
            d['coord'].requires_grad_(True)
        _out = self.model(d)
        ret = dict(energy=_out['energy'].item(), charges=_out['charges'].detach()[:-1].cpu().numpy())
        if forces:
            if 'forces' in _out:
                f = _out['forces'][:-1]
            else:
                f = - torch.autograd.grad(_out['energy'], d['coord'])[0][:-1]
            ret['forces'] = f.detach().cpu().numpy()
        torch._C._set_grad_enabled(prev)
        return ret

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        _in = self._make_input()
        do_forces = 'forces' in properties
        _out =  self._eval_model(_in, do_forces)

        self.results['energy'] = _out['energy']
        self.results['charges'] = _out['charges']
        if do_forces:
            self.results['forces'] = _out['forces']

    # QR API functions
    def check_trained_atoms(self):
         if not np.in1d(self.atoms.numbers, self.species).all():
            raise NotImplementedError("Unfortunately, we do not have a trained model for all elements in your system.")

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
        self.atoms = atoms
        self.check_trained_atoms()
        atoms.calc = self
        self.do_reset()
        self.set_charge(charge)
        unit_convert = ase_units.kcal / ase_units.mol
        self.calculate(atoms, properties=['energy', 'forces'])
        self.energy_free = self.results['energy'] * unit_convert
        self.forces = self.results['forces'].astype(np.float64) * unit_convert

        if 0: # we need debugging flag here to switch on and off.
            print(("AIMNet2: ",self.method))
            print(('Energy:',self.energy_free))
            print(('Force:', self.forces))
            
