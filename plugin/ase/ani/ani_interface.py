import torch
import numpy as np
from ase.calculators import calculator
from ase import units

# from torchani
def map2central(cell, coordinates, pbc):
    """Map atoms outside the unit cell into the cell using PBC.
    Arguments:
        cell (:class:`torch.Tensor`): tensor of shape (3, 3) of the three
            vectors defining unit cell:
            .. code-block:: python
                tensor([[x1, y1, z1],
                        [x2, y2, z2],
                        [x3, y3, z3]])
        coordinates (:class:`torch.Tensor`): Tensor of shape
            ``(molecules, atoms, 3)``.
        pbc (:class:`torch.Tensor`): boolean vector of size 3 storing
            if pbc is enabled for that direction.
    Returns:
        :class:`torch.Tensor`: coordinates of atoms mapped back to unit cell.
    """
    # Step 1: convert coordinates from standard cartesian coordinate to unit
    # cell coordinates
    inv_cell = torch.inverse(cell)
    coordinates_cell = torch.matmul(coordinates, inv_cell)
    # Step 2: wrap cell coordinates into [0, 1)
    coordinates_cell -= coordinates_cell.floor() * pbc
    # Step 3: convert from cell coordinates back to standard cartesian
    # coordinate
    return torch.matmul(coordinates_cell, cell)

class ANI_interface(object):
    def __init__(self, module):
        self.module = module
        self.device = next(module.parameters()).device


    def run(self, flag, R, Z, pbc, cell):
        R, Z, pbc, cell = self.prepare_in(R, Z, pbc, cell)

        arrays_out = self.calculate(flag, R, Z, pbc, cell)
        return arrays_out
        
    def prepare_in(self, R, Z, pbc, cell):
    
        natom = Z.shape[-1]
        R = R.reshape(-1, natom, 3)
        Z = Z.reshape(-1, natom)
        nbatch = R.shape[0]
        assert Z.shape[0] == nbatch
        pbc = pbc.reshape(3)
        cell = cell.reshape(3, 3)

        return R, Z, pbc, cell

    def calculate(self, flag, R, Z, pbc, cell):
        pbc_enabled = pbc.any()

        tensors_out = list()

        R = torch.as_tensor(R, dtype=torch.float32, device=self.device)
        Z = torch.as_tensor(Z, dtype=torch.int64, device=self.device)
        
        if flag > 0:
            R.requires_grad_(True)

        if pbc_enabled:
            pbc = torch.as_tensor(pbc, dtype=torch.bool, device=self.device)
            cell = torch.as_tensor(
                cell, dtype=torch.float32, device=self.device)
            R = map2central(cell, R, pbc)
            if flag == 2:
                assert R.shape[0] == 1, 'Batch stress is not implemented.'
                scaling = torch.eye(3, requires_grad=True,
                                    dtype=torch.float32, device=self.device)
                cell = cell @ scaling
                R = R @ scaling
            E = self.module((Z, R), cell=cell, pbc=pbc).energies
        else:
            E = self.module((Z, R)).energies
        tensors_out.append(E)

        if flag > 0:
            F = - torch.autograd.grad(E, R, retain_graph=flag == 2)[0]
            tensors_out.append(F)

        if flag == 2:
            V = cell.det().abs()
            stress = torch.autograd.grad(E, scaling)[0] / V
            stress = stress.squeeze(0)
            tensors_out.append(stress)

        arrays = [t.detach().cpu().numpy() for t in tensors_out]
        return arrays
        

class ANIRPCCalculator(calculator.Calculator):
    implemented_properties = ['energy', 'forces', 'stress', 'free_energy']

    def __init__(self,model):
        super(ANIRPCCalculator, self).__init__()
        self.interface = ANI_interface(model)

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=calculator.all_changes):
        super(ANIRPCCalculator, self).calculate(
            atoms, properties, system_changes)

        calc_type = 0
        if 'forces' in properties:
            calc_type = 1
        if 'stress' in properties:
            calc_type = 2
        
    
        R = self.atoms.get_positions()[None, ...].astype(np.float32)
        Z = self.atoms.get_atomic_numbers()[None, ...].astype(np.uint8)
        pbc = self.atoms.get_pbc().astype(np.bool_)
        cell = self.atoms.get_cell(complete=True).astype(np.float32)
        

        results = self.interface.run(calc_type, R, Z, pbc, cell)
        self.results['energy'] = self.results['free_energy'] = results[0][0] * units.Hartree
        if calc_type == 1:
            self.results['forces'] = results[1][0] * units.Hartree
        if calc_type == 2:
            self.results['stress'] = results[2] * units.Hartree

