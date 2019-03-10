"""
The TorchANI calculator can compute an energy and force for a molecular system using a Artificial Neural Network .

TorchANI is developed by Xiang Gao (@zasdfgbnm) and can be found at:
https://github.com/aiqm/torchani

TorchANI is based on PyTorch which is a open source library for deep neural networks.
https://pytorch.org/

TorchANI is a pytorch implementation of the ANI models:
https://github.com/isayev/ASE_ANI

We are using the ASE calculator interface:
https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html

"""
import os
import numpy as np
import torch
import torchani
from sets import Set
from ase.calculators.general import Calculator


# If you have a GPU available, then PyTorch will use it, otherwise the CPU is used.
device_str = 'cuda' if torch.cuda.is_available() else 'cpu'
device = torch.device(device_str )


class TorchAni(Calculator):
  """
  This class interfaces to the TorchANI package to compute the energy and gradients using the ANI neural network.

  Models available:
  - ani-1ccx_8x    https://chemrxiv.org/articles/Outsmarting_Quantum_Chemistry_Through_Transfer_Learning/6744440/1
  - ani-1x_8x      https://aip.scitation.org/doi/10.1063/1.5023802
  - ani-1          https://pubs.rsc.org/en/content/articlelanding/2017/sc/c6sc05720a#!divAbstract

  Attributes:
      label (str) a useful name for the ASE calculator.
      atoms (ase.atoms.Atoms) set of atoms in the molecular system.
      coordinates(numpy.ndarray) get the positions as x,y and z coordinates in Angstrom.
      energy_free (float) gets the energy as a sum of atomic contributions.
      forces (numpy.ndarray) stores the derivative of the energy with respect to the x,y,z coordinates.
  """
  def __init__(self,method='ani-1x_8x',label="ase",atoms=None,coordinates='tmp_ase.pdb',**kwargs):

    self.label  = label
    self.atoms  = atoms

    coordinates = os.path.dirname(label)+"/"+ coordinates
    self.coordinates = coordinates

    self.method = method

    if self.method == 'ani-1ccx':
        self.model = torchani.models.ANI1ccx()
        # TODO: There is a new api for this in master branch of torchani.
        # calculator = torchani.models.ANI1ccx().ase()
        # self.atoms.set_calculator(calculator)
    else: #'ani-1x_8x'
        self.model = torchani.models.ANI1x()

    self.energy_free = None
    self.forces = []


  def check_trained_atoms(self):
      """
      The ANN models are only trained on H,C,N,O
      """
      trained_atoms = Set(['C', 'H', 'N', 'O'])
      if not Set(self.atoms.get_chemical_symbols()).issubset(trained_atoms):
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
    self.coordinates = coordinates
    self.charge = charge
    self.pointcharges = pointcharges
    self.check_trained_atoms()

    atoms_symbols = self.atoms.get_chemical_symbols()
    atoms_symbols = self.model.species_to_tensor(atoms_symbols).to(device).unsqueeze(0)

    xyz = self.atoms.get_positions().tolist()
    coords = torch.tensor([xyz],requires_grad = True, device = device)

    _, energy = self.model((atoms_symbols, coords))
    derivative = torch.autograd.grad(energy.sum(), coords)[0]
    force = -derivative

    self.energy_free = energy.item()
    self.forces = force.squeeze().numpy().astype(np.float64)

    if 0: # we need debugging flag here to switch on and off.
        print("Torch ANI: ",self.method)
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

  def set_method(self, method):
    self.method = method

