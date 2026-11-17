"""
AIMNet2 plugin for Q|R.
Based on the AIMNet2 ASE calculator by the authors of the AIMNet2 model.

Two AIMNet2 python packages exist:

  * ``aimnet``      -- isayevlab/aimnetcentral, the current package (preferred)
  * ``aimnet2calc`` -- zubatyuk/aimnet2calc, the original package (legacy)

Both are supported.  The AQuaRef model ``aimnet2-qr`` is a legacy TorchScript
(.jpt) model that is not in the ``aimnet`` model registry, so this module
resolves it itself: it looks for the weights in the ``aimnet`` cache directory
(``$AIMNET_CACHE_DIR`` or ``~/.cache/aimnet``), then in
``${PREFIX}/share/aimnet`` of the conda environment (where the Phenix package
installs them on Linux and macOS).  If neither has them, with ``aimnet`` this
module downloads them from the original model zoo into the cache directory;
with ``aimnet2calc``, aimnet2calc looks in its own assets folder and downloads
them there.
"""
from __future__ import absolute_import, division, print_function
import os
import sys
import warnings
import numpy as np
import ase.units as ase_units
from libtbx.env_config import get_conda_prefix
from libtbx.utils import Sorry

warnings.filterwarnings('ignore', module='aimnet')
warnings.filterwarnings('ignore', module='aimnet2calc')

try:
  from aimnet.calculators import AIMNet2ASE
  from aimnet.calculators import AIMNet2Calculator as _AIMNet2Base
  AIMNET_BACKEND = 'aimnet'
except ModuleNotFoundError:
  from aimnet2calc import AIMNet2ASE  # legacy package
  AIMNET_BACKEND = 'aimnet2calc'

# Legacy (.jpt) models that the ``aimnet`` package does not know about:
#   name -> (file name in cache_dir() or share_dir(), download URL)
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

def cache_dir():
  """
  The ``aimnet`` model cache: ``$AIMNET_CACHE_DIR``, else ``~/.cache/aimnet``
  (the rule of aimnet's get_cache_dir, except that an empty value means the
  default; aimnet fails on it).  Not created here, so that a read-only or
  missing home directory does not stop the lookup in share_dir().
  """
  return (os.environ.get('AIMNET_CACHE_DIR') or
          os.path.join(os.path.expanduser('~'), '.cache', 'aimnet'))

def share_dir():
  """
  ``${PREFIX}/share/aimnet`` of the conda environment, where the Phenix
  package installs the legacy model files on Linux and macOS;
  ``Library\\share\\aimnet`` on Windows (no files installed there yet), as for
  ``share/cctbx`` in libtbx.env_config.get_installed_path.
  get_conda_prefix() also handles macOS python.app (the Phenix GUI).
  """
  if sys.platform == 'win32':
    return os.path.join(sys.prefix, 'Library', 'share', 'aimnet')
  return os.path.join(get_conda_prefix(), 'share', 'aimnet')

def resolve_legacy_model(model):
  """
  Map a legacy model name (e.g. ``aimnet2-qr``) to a local .jpt file:
  cache_dir(), then share_dir(); if neither has it, the ``aimnet`` backend
  downloads it into cache_dir() and the ``aimnet2calc`` backend gets the name
  back (aimnet2calc then uses its own assets folder or downloads there).
  Anything else (registry names, Hugging Face ids, explicit paths) is passed
  through unchanged.
  """
  if model not in LEGACY_MODELS: return model
  file_name, url = LEGACY_MODELS[model]
  folders = [cache_dir(), share_dir()]
  for folder in folders:
    path = os.path.join(folder, file_name)
    if os.path.isfile(path): return path
  if AIMNET_BACKEND == 'aimnet2calc': return model
  path = os.path.join(folders[0], file_name)
  print("AIMNet2 model '%s' not found in %s, downloading from %s" % (
    model, " or ".join(folders), url))
  try:
    from urllib.request import urlretrieve
    os.makedirs(folders[0], exist_ok=True)
    tmp = path + '.part'
    urlretrieve(url, tmp)
    os.replace(tmp, path)
  except Exception as e:
    raise Sorry(
      "Could not download the AIMNet2 model '%s' (%s).\n"
      "Place the file %s in %s (or set AIMNET_CACHE_DIR) or in %s and try "
      "again." % (model, str(e), file_name, folders[0], folders[1]))
  return path


class AIMNet2Calculator(AIMNet2ASE):
    """ Modification of the AIMNet2ASE class to work with Q|R.
    """
    def __init__(self, model='aimnet2-qr', charge=0, mult=1):
        model = resolve_legacy_model(model)
        if AIMNET_BACKEND == 'aimnet':
          _ensure_ptxas()
          # nb_threshold=0 forces the explicit neighbour-list code path for
          # systems of any size.  Legacy .jpt models require an explicit
          # neighbour matrix, and with the default threshold the ``aimnet``
          # calculator does not build one for systems of <= 120 atoms (which
          # Q|R clusters routinely are).
          base = _AIMNet2Base(model, nb_threshold=0)
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
