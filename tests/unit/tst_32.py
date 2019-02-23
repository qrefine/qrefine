import run_tests
from ase import Atoms
from qrefine.plugin.ase.torchani_qr import TorchAni
import warnings
import unittest

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
  except ImportError:
      with warnings.catch_warnings():
          warnings.simplefilter("ignore")
          import torchani
          torchani_installed = True

  prefix="tst_32"
  if torchani_installed:
      rc = run_tests.runner(function=run, prefix=prefix, disable=False)
  else:
      rc = run_tests.runner(function=run, prefix=prefix, disable=True)
  assert not rc, '%s rc: %s' % (prefix, rc)