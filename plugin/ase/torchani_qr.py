import os
import numpy as np
import torch
import torchani
from ase.calculators.general import Calculator

device = torch.device('cpu')

anipath = os.path.dirname(__file__)
const_file = anipath + '/ani/ani-1x_dft_x8ens/rHCNO-5.2R_16-3.5A_a4-8.params'
sae_file = anipath + '/ani/ani-1x_dft_x8ens/sae_linfit.dat'
network_dir = anipath + '/ani/ani-1x_dft_x8ens/train'


class TorchAni(Calculator):

  def __init__(self,label="ase",atoms=None,coordinates='tmp_ase.pdb',**kwargs):
    self.label=label
    coordinates=os.path.dirname(label)+"/"+ coordinates
    self.coordinates=coordinates
    self.atoms = atoms
    self.energy_free = None
    self.forces = []
    self.aev_computer = torchani.SortedAEV(const_file=const_file, device=device)
    self.nn = torchani.ModelOnAEV(self.aev_computer, derivative=True,
                         from_nc=network_dir, ensemble=8)
    self.shift_energy = torchani.EnergyShifter(sae_file)

  def run_qr(self, atoms, coordinates, charge, pointcharges, command=None, define_str=None):
    self.atoms = atoms
    self.coordinates = coordinates
    self.charge = charge
    self.pointcharges = pointcharges
    atoms_symbols = self.atoms.get_chemical_symbols()
    xyz = self.atoms.get_positions()
    coords = torch.tensor([xyz], dtype=self.aev_computer.dtype, device=self.aev_computer.device)
    energy, derivative = self.nn(coords,atoms_symbols)
    energy = self.shift_energy.add_sae(energy, atoms_symbols)
    force = -derivative
    self.forces = force.squeeze().numpy().astype(np.float64)
    self.energy_free = energy.item()
    if 0:
        print(" RUNNING Torch ANI")
        print('Energy:',self.energy_free)
        print('Force:', self.forces)

  def get_command(self):
        return "TorchANI"

  def set_label(self, label):
      self.label = label

