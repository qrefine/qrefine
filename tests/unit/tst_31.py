import run_tests
from ase import Atoms
import warnings
from libtbx.test_utils import approx_equal
import os

def run(prefix):
    """
    Exercise TorchANI calculator works as expected for both ANI1x and ANI1ccx
    """

    atoms = Atoms('CHHHH', [[0.03192167, 0.00638559, 0.01301679],
                            [-0.83140486, 0.39370209, -0.26395324],
                            [-0.66518241, -0.84461308, 0.20759389],
                            [0.45554739, 0.54289633, 0.81170881],
                            [0.66091919, -0.16799635, -0.91037834]])

    #Test the ANIx model
    calculatorANIx = TorchAni(method='ani-1x',label="ase",atoms=atoms )

    calculatorANIx.run_qr(atoms=atoms,coordinates=atoms.get_positions(),charge=0,pointcharges=None )

    assert calculatorANIx.method == 'ani-1x'
    assert approx_equal(calculatorANIx.energy_free,  -40.459022522)
    assert approx_equal(calculatorANIx.forces , [[ 0.0306, -0.1316, -0.0527],
        [-0.1293,  0.1639, -0.0774],
        [ 0.0856, -0.0429,  0.0408],
        [ 0.0268,  0.0060,  0.0381],
        [-0.0138,  0.0046,  0.0511]],0.001)

    #Test the ANIccx model

    calculatorANIccx = TorchAni(method='ani-1ccx', label="ase", atoms=atoms)

    calculatorANIccx.run_qr(atoms=atoms, coordinates=atoms.get_positions(), charge=0, pointcharges=None)

    assert calculatorANIccx.method == 'ani-1ccx'
    assert approx_equal(calculatorANIccx.energy_free, -40.4256210327)
    assert approx_equal(calculatorANIccx.forces, [[0.0312, -0.1272, -0.0511],
                                                [-0.1200, 0.1628, -0.0761],
                                                [0.0856, -0.0448, 0.0407],
                                                [0.0219, 0.0044, 0.0343],
                                                [-0.0187, 0.0049, 0.0521]], 0.001)

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
