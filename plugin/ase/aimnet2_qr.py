"""
AIMNet2 plugin for Q|R.
Based on the AIMNet2 ASE calculator by the authors of the AIMNet2 model.

Two AIMNet2 python packages exist:

  * ``aimnet``      -- isayevlab/aimnetcentral, the current package (preferred)
  * ``aimnet2calc`` -- zubatyuk/aimnet2calc, the original package (legacy)

Both are supported.  The AQuaRef model ``aimnet2-qr`` is a legacy TorchScript
(.jpt) model that is not in the ``aimnet`` model registry, so this module
resolves it itself: it looks for the weights in the ``aimnet`` cache directory
(``$AIMNET_CACHE_DIR`` or ``~/.cache/aimnet``) and downloads them from the
original model zoo if they are missing, exactly as ``aimnet2calc`` used to.
"""
from __future__ import absolute_import, division, print_function
import os
import sys
import warnings
import numpy as np
import ase.units as ase_units
from libtbx.utils import Sorry

warnings.filterwarnings('ignore', module='aimnet')
warnings.filterwarnings('ignore', module='aimnet2calc')

try:
  from aimnet.calculators import AIMNet2ASE
  from aimnet.calculators import AIMNet2Calculator as _AIMNet2Base
  from aimnet.calculators.model_registry import get_cache_dir
  AIMNET_BACKEND = 'aimnet'
except ModuleNotFoundError:
  from aimnet2calc import AIMNet2ASE  # legacy package
  AIMNET_BACKEND = 'aimnet2calc'

# Legacy (.jpt) models that the ``aimnet`` package does not know about:
#   name -> (file name in the aimnet cache directory, download URL)
LEGACY_MODELS = {
  'aimnet2-qr': (
    'aimnet2-qr_b97md4_qzvp_2.jpt',
    'https://github.com/zubatyuk/aimnet-model-zoo/raw/main/aimnet2-qr/aimnet2-qr_b97md4_qzvp_2.jpt'),
}

def _ensure_ptxas():
  """
  The ``aimnet`` neighbour-list kernels are JIT-compiled with Triton, which
  needs the CUDA ``ptxas`` binary.  Phenix ships it in its conda environment
  but that ``bin`` directory is not on PATH, so point Triton at it.
  """
  if os.environ.get('TRITON_PTXAS_PATH'): return
  ptxas = os.path.join(sys.prefix, 'bin', 'ptxas')
  if os.path.isfile(ptxas):
    os.environ['TRITON_PTXAS_PATH'] = ptxas

def resolve_legacy_model(model):
  """
  Map a legacy model name (e.g. ``aimnet2-qr``) to a local .jpt file,
  downloading it into the ``aimnet`` cache directory on first use.
  Anything else (registry names, Hugging Face ids, explicit paths) is passed
  through to the ``aimnet`` resolver unchanged.
  """
  if model not in LEGACY_MODELS: return model
  file_name, url = LEGACY_MODELS[model]
  path = os.path.join(get_cache_dir(), file_name)
  if os.path.isfile(path): return path
  print("AIMNet2 model '%s' not found in %s, downloading from %s" % (
    model, get_cache_dir(), url))
  try:
    from urllib.request import urlretrieve
    tmp = path + '.part'
    urlretrieve(url, tmp)
    os.replace(tmp, path)
  except Exception as e:
    raise Sorry(
      "Could not download the AIMNet2 model '%s' (%s).\n"
      "Place the file %s in %s (or set AIMNET_CACHE_DIR) and try again." % (
        model, str(e), file_name, get_cache_dir()))
  return path


class AIMNet2Calculator(AIMNet2ASE):
    """ Modification of the AIMNet2ASE class to work with Q|R.
    """
    def __init__(self, model='aimnet2-qr', charge=0, mult=1):
        if AIMNET_BACKEND == 'aimnet':
          _ensure_ptxas()
          # nb_threshold=0 forces the explicit neighbour-list code path for
          # systems of any size.  Legacy .jpt models require an explicit
          # neighbour matrix, and with the default threshold the ``aimnet``
          # calculator does not build one for systems of <= 120 atoms (which
          # Q|R clusters routinely are).
          base = _AIMNet2Base(resolve_legacy_model(model), nb_threshold=0)
          super(AIMNet2Calculator, self).__init__(base, charge=charge, mult=mult)
        else:
          super(AIMNet2Calculator, self).__init__(model, charge=charge, mult=mult)

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
        atoms.calc = self
        self.atoms = atoms
        self.set_charge(charge)
        unit_convert = ase_units.kcal / ase_units.mol
        self.calculate(atoms, properties=['energy', 'forces'])
        # ``aimnet`` stores the energy as a scalar, ``aimnet2calc`` as a
        # 1-element array; .item() handles both.
        self.energy_free = np.asarray(self.results['energy']).item() * unit_convert
        self.forces = self.results['forces'].astype(np.float64) * unit_convert

        if 0: # we need debugging flag here to switch on and off.
            # print(("AIMNet2: ",self.method))
            print(('Energy:',self.energy_free))
            print(('Force:', self.forces))
