"""
The TorchANI calculator can compute an energy and force for a molecular system using a Artificial Neural Network .

TorchANI is developed by Xiang Gao (@zasdfgbnm) and can be found at:
https://github.com/aiqm/torchani

TorchANI is based on PyTorch which is a open source library for deep neural networks.
https://pytorch.org/

We are using the ASE calculator interface:
https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html

"""
import os
import numpy as np
import torch
import torchani
from ase.calculators.general import Calculator

device_str = 'cuda' if torch.cuda.is_available() else 'cpu'
device = torch.device(device_str )

anipath = os.path.dirname(__file__)
const_file  = anipath + '/ani/ani-1x_dft_x8ens/rHCNO-5.2R_16-3.5A_a4-8.params'
sae_file    = anipath + '/ani/ani-1x_dft_x8ens/sae_linfit.dat'
network_dir = anipath + '/ani/ani-1x_dft_x8ens/train'

class TorchAni(Calculator):
  """
  This class is one strategy to compute the energy and gradients using the ANI neural network.

  If you have a GPU available, then PyTorch will use it, otherwise the CPU is used.


  Attributes:
      label (str) a useful name for the ASE calculator.
      atoms (ase.atoms.Atoms) of atom
      coordinates(numpy.ndarray) get the positions as x,y and z coordinates in Angstroem.
      energy_free (float) gets the energy as a sum of atomic contributions.
      forces (numpy.ndarray) stores the derivative of the energy with respect to xyz coordinates.
  """
  def __init__(self,label="ase",atoms=None,coordinates='tmp_ase.pdb',**kwargs):
    self.label  = label
    self.atoms  = atoms
    coordinates = os.path.dirname(label)+"/"+ coordinates
    self.coordinates = coordinates


    self.aev_computer = torchani.SortedAEV(const_file=const_file, device=device)
    self.nn = torchani.ModelOnAEV(self.aev_computer, derivative=True,
                         from_nc=network_dir, ensemble=1)
    self.shift_energy = torchani.EnergyShifter(sae_file)

    self.energy_free = None
    self.forces = []

  def run_qr(self, atoms, coordinates, charge, pointcharges, command=None, define_str=None):
    """
    This method is called every time an energy and forces are needed.
    The Q|R code calls this method at each step of LBFGS.

    Args:
      atoms (ase.atoms.Atoms) an updated set of atoms.
      TODO: coordinates (numpy.ndarray) an updated set of coordinates. are these even being used?
      charge (int) the charge on the molecular system.
      TODO: pointcharges (?) are these even being used?
      command (str) not used in this calcualtor
      define_str() not used in this calculator, legacy from Turbomole?

    """
    self.atoms = atoms
    self.coordinates = coordinates
    self.charge = charge
    self.pointcharges = pointcharges
    atoms_symbols = self.atoms.get_chemical_symbols()
    xyz = self.atoms.get_positions()
    coords = torch.tensor([xyz], dtype=self.aev_computer.dtype, device=self.aev_computer.device)
    energy, derivative = self.nn(coords,atoms_symbols)
    energy = self.shift_energy.add_sae(energy, atoms_symbols)
    self.energy_free = energy.item()
    force = -derivative
    self.forces = force.squeeze().numpy().astype(np.float64)

    if 0: # we need debugging flag here to switch on and off.
        print(" RUNNING Torch ANI")
        print('Energy:',self.energy_free)
        print('Force:', self.forces)

  def get_command(self):
    """
    This command is not used for running the ASE calcualtor.
    It is simply used here to satisfy the base class.
    :return:
        command (str) just the name of the calculator, not the command to actually run as seen in other calculators.
    """
    return "TorchANI"

  def set_label(self, label):

    self.label = label

