from __future__ import absolute_import
from qrefine.tests.unit import run_tests
from ase import Atoms
import warnings
import unittest
import os

class atoms_present_test(unittest.TestCase):

    def setUp(self):
        self.atoms = Atoms('CHHHH', [[0.03192167, 0.00638559, 0.01301679],
                            [-0.83140486, 0.39370209, -0.26395324],
                            [-0.66518241, -0.84461308, 0.20759389],
                            [0.45554739, 0.54289633, 0.81170881],
                            [0.66091919, -0.16799635, -0.91037834]])

        # This has a S atom in the system.
        # This should cause the calculator to raise an error.
        self.new_atoms = Atoms('CHHHS', [[0.03192167, 0.00638559, 0.01301679],
                                     [-0.83140486, 0.39370209, -0.26395324],
                                     [-0.66518241, -0.84461308, 0.20759389],
                                     [0.45554739, 0.54289633, 0.81170881],
                                     [0.66091919, -0.16799635, -0.91037834]])

        self.calculator = TorchAni(method='ani-1x_8x', label="ase", atoms=self.atoms)


    def test_good(self,):
        self.calculator.run_qr(atoms=self.atoms, coordinates=self.atoms.get_positions(), charge=0, pointcharges=None)

    def test_bad(self,):
       with self.assertRaises(NotImplementedError):
           self.calculator.run_qr(atoms=self.new_atoms, coordinates=self.atoms.get_positions(), charge=0, pointcharges=None)

def run(prefix):
    """
    Assert that TorchANI calculator raises an error when model can't work on atom types in molecular system.
    """
    unittest.main()


if(__name__ == "__main__"):
  torchani_installed = False
  try:
    import torch
    device = torch.device('cpu')
    import torchnani
    from qrefine.plugin.ase.torchani_qr import TorchAni
    torchani_installed = True
  except ImportError:
    pass
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=not torchani_installed)
