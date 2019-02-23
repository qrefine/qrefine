import run_tests
from ase import Atoms
from qrefine.plugin.ase.torchani_qr import TorchAni

try:
    import torch
    import torchnani
    torchani_installed = True
except ImportError:
    torchani_installed = False

device = torch.device('cpu')

def run(prefix):
    """
    Exercise TorchANI calculator works as expected.
    """

    atoms = Atoms('CHHHH', [[0.03192167, 0.00638559, 0.01301679],
                            [-0.83140486, 0.39370209, -0.26395324],
                            [-0.66518241, -0.84461308, 0.20759389],
                            [0.45554739, 0.54289633, 0.81170881],
                            [0.66091919, -0.16799635, -0.91037834]])

    calculator = TorchAni(label="ase",atoms=atoms )



    print "atoms", atoms
    print "coordinate", atoms.get_positions()

    calculator.run_qr(atoms=atoms,coordinates=atoms.get_positions(),charge=0,pointcharges=None )

    print "energy is :", calculator.energy_free
    print "forces are :", calculator.forces

    assert calculator.energy_free == -40.425621032714844
    assert calculator.forces == [[-0.00555373, -0.00059946,  0.00076646],
                                 [ 0.01750624, -0.01134875,  0.00614653],
                                 [ 0.00449754,  0.01362024, -0.00396855],
                                 [-0.00801547, -0.00424718, -0.00878889],
                                 [-0.00843458,  0.00257515,  0.00584446]]

if(__name__ == "__main__"):
  prefix="tst_31"
  if torchani_installed:
      rc = run_tests.runner(function=run, prefix=prefix, disable=True)
  else:
      rc = run_tests.runner(function=run, prefix=prefix, disable=False)
  assert not rc, '%s rc: %s' % (prefix, rc)