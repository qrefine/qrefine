import run_tests
import os
from ase import Atoms
from qrefine.plugin.ase.pyscf_qr import Pyscf
import warnings
from libtbx.test_utils import approx_equal

def run(prefix):
    """
    Assert that PySCF can compute an Energy and Gradient for a given molecular system.
    """


    atoms = Atoms('CHHHH', [[0.03192167, 0.00638559, 0.01301679],
                        [-0.83140486, 0.39370209, -0.26395324],
                        [-0.66518241, -0.84461308, 0.20759389],
                        [0.45554739, 0.54289633, 0.81170881],
                        [0.66091919, -0.16799635, -0.91037834]])

    calculator = Pyscf(method='hf', label="ase", atoms=atoms)

    calculator.run_qr(atoms, coordinates=None, charge=None, pointcharges=None, command=None, define_str=None)

    assert approx_equal(calculator.energy_free , -24894.129755)

    assert  approx_equal(calculator.forces.tolist(),[[  10.3539855,  -109.02533792,  -52.22553528],
 [ -89.1104867 ,  135.541275 ,   -63.22351153],
 [  81.34525894,  -32.53859855,   36.0316453 ],
 [  17.81046893 ,   0.39223847,   28.47319708],
 [ -20.39922667 ,   5.630423 ,    50.94420443]])



if(__name__ == "__main__"):
  pyscf_installed = False
  try:
    import pyscf
    from pyscf import gto, scf, grad, dft
    pyscf_installed = True
  except ImportError:
    print "Pyscf not installed. Please run: phenix.python pip -m  install pyscf"
  prefix = os.path.basename(__file__).replace(".py","")
  if pyscf_installed:
    run_tests.runner(function=run, prefix=prefix, disable=False)
  else:
    run_tests.runner(function=run, prefix=prefix, disable=True)
