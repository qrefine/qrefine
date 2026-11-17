from __future__ import absolute_import, division, print_function
import importlib, os, shutil, sys, tempfile, types
import urllib.request
from libtbx.utils import Sorry

# Where plugin/ase/aimnet2_qr.py looks for the legacy AIMNet2 weights: the
# aimnet cache ($AIMNET_CACHE_DIR or ~/.cache/aimnet), then
# ${PREFIX}/share/aimnet (installed by the Phenix package); if neither has
# them, Q|R downloads them into the cache on the aimnet backend and leaves
# the lookup to aimnet2calc on the aimnet2calc backend. Both packages are
# replaced by stubs, so this runs without torch and without network access.

PLUGIN = "qrefine.plugin.ase.aimnet2_qr"
FILE = "aimnet2-qr_b97md4_qzvp_2.jpt"

class Recorder(object):
  """Stands in for the AIMNet2 calculator classes: records its arguments."""
  def __init__(self, *args, **kwargs):
    self.args = args
    self.kwargs = kwargs

def aimnet_get_cache_dir():
  """aimnet.calculators.model_registry.get_cache_dir (aimnet 0.2.x)."""
  cache_dir = os.environ.get("AIMNET_CACHE_DIR")
  if cache_dir is None:
    cache_dir = os.path.join(os.path.expanduser("~"), ".cache", "aimnet")
  os.makedirs(cache_dir, exist_ok=True)
  return cache_dir

def stub(name, **attrs):
  m = types.ModuleType(name)
  m.__dict__.update(attrs)
  return m

def load_plugin(backend, share=None):
  """Import the plugin afresh on a stub backend; share_dir() returns share."""
  for k in list(sys.modules):
    if k == PLUGIN or k.split(".")[0] in ("aimnet", "aimnet2calc"):
      del sys.modules[k]
  if backend == "aimnet":
    sys.modules["aimnet"] = stub("aimnet")
    sys.modules["aimnet.calculators"] = stub("aimnet.calculators",
      AIMNet2ASE=Recorder, AIMNet2Calculator=Recorder)
    sys.modules["aimnet.calculators.model_registry"] = stub(
      "aimnet.calculators.model_registry", get_cache_dir=aimnet_get_cache_dir)
  else:
    sys.modules["aimnet"] = None # import aimnet raises ModuleNotFoundError
    sys.modules["aimnet2calc"] = stub("aimnet2calc", AIMNet2ASE=Recorder)
  p = importlib.import_module(PLUGIN)
  assert p.AIMNET_BACKEND == backend, p.AIMNET_BACKEND
  if share is not None: p.share_dir = lambda: share
  return p

def put(folder):
  if not os.path.isdir(folder): os.makedirs(folder)
  path = os.path.join(folder, FILE)
  with open(path, "w") as f: f.write(folder)
  return path

class fake_urlretrieve(object):
  """Replaces urllib.request.urlretrieve; records the URLs it is given."""
  def __init__(self, error=None):
    self.error = error
    self.urls = []
  def __call__(self, url, file_name):
    self.urls.append(url)
    if self.error is not None: raise self.error
    with open(file_name, "w") as f: f.write("weights")
  def __enter__(self):
    self.real = urllib.request.urlretrieve
    urllib.request.urlretrieve = self
    return self
  def __exit__(self, *exc):
    urllib.request.urlretrieve = self.real

def expect_sorry(f, *has):
  try: f()
  except Sorry as e:
    for s in has: assert s in str(e), (s, str(e))
  else: raise AssertionError("expected Sorry")

def check_search_order(backend, root):
  cache = os.path.join(root, backend, "cache") # not created by the lookup
  share = os.path.join(root, backend, "share")
  os.environ["AIMNET_CACHE_DIR"] = cache
  p = load_plugin(backend, share)
  # registry names, Hugging Face ids and paths pass through
  for name in ("aimnet2", "isayevlab/aimnet2-wb97m-d3", "/x/model.jpt"):
    assert p.resolve_legacy_model(name) == name
  # share only; the lookup does not create the cache folder
  in_share = put(share)
  assert p.resolve_legacy_model("aimnet2-qr") == in_share
  assert not os.path.exists(cache)
  # cache and share: the cache wins
  in_cache = put(cache)
  assert p.resolve_legacy_model("aimnet2-qr") == in_cache
  # the calculator is built on the resolved file
  c = p.AIMNet2Calculator("aimnet2-qr")
  if backend == "aimnet":
    assert c.args[0].args == (in_cache,)
    assert c.args[0].kwargs == {"nb_threshold": 0}
  else:
    assert c.args == (in_cache,)
  os.remove(in_cache); os.remove(in_share)
  return p, cache, share

def check_cache_not_creatable(backend, root):
  """$HOME read-only or missing: the share copy is still found."""
  blocker = os.path.join(root, backend + "_blocker")
  with open(blocker, "w") as f: f.write("a file, so no folder below it")
  cache = os.path.join(blocker, "cache")
  share = os.path.join(root, backend, "share_only")
  os.environ["AIMNET_CACHE_DIR"] = cache
  p = load_plugin(backend, share)
  in_share = put(share)
  assert p.resolve_legacy_model("aimnet2-qr") == in_share
  os.remove(in_share)
  if backend == "aimnet":
    with fake_urlretrieve() as fake:
      expect_sorry(lambda: p.resolve_legacy_model("aimnet2-qr"),
        FILE, cache, share)
    assert fake.urls == []
  else:
    assert p.resolve_legacy_model("aimnet2-qr") == "aimnet2-qr"

def test_aimnet(root):
  p, cache, share = check_search_order("aimnet", root)
  # in neither folder: downloaded into the cache (created now), not share
  with fake_urlretrieve() as fake:
    assert p.resolve_legacy_model("aimnet2-qr") == os.path.join(cache, FILE)
  assert fake.urls == [p.LEGACY_MODELS["aimnet2-qr"][1]]
  assert sorted(os.listdir(cache)) == [FILE]
  assert os.listdir(share) == []
  # a failed download names the file and both folders
  os.remove(os.path.join(cache, FILE))
  with fake_urlretrieve(error=IOError("offline")):
    expect_sorry(lambda: p.resolve_legacy_model("aimnet2-qr"),
      "offline", FILE, cache, share)
  check_cache_not_creatable("aimnet", root)

def test_aimnet2calc(root):
  p, cache, share = check_search_order("aimnet2calc", root)
  # in neither folder: aimnet2calc resolves the name itself (its assets
  # folder, else its own download); nothing is written here
  with fake_urlretrieve() as fake:
    assert p.resolve_legacy_model("aimnet2-qr") == "aimnet2-qr"
  assert fake.urls == []
  assert os.listdir(cache) == [] and os.listdir(share) == []
  check_cache_not_creatable("aimnet2calc", root)

def test_cache_dir():
  """aimnet's rule, except that an empty value means the default."""
  default = os.path.join(os.path.expanduser("~"), ".cache", "aimnet")
  for backend in ("aimnet", "aimnet2calc"):
    p = load_plugin(backend)
    os.environ.pop("AIMNET_CACHE_DIR", None)
    assert p.cache_dir() == default
    os.environ["AIMNET_CACHE_DIR"] = ""
    assert p.cache_dir() == default
    os.environ["AIMNET_CACHE_DIR"] = "/x/y"
    assert p.cache_dir() == "/x/y"

def test_share_dir(root):
  """${PREFIX}/share/aimnet for fake conda prefixes on each platform."""
  p = load_plugin("aimnet2calc")
  prefix = os.path.join(root, "conda_prefix")
  os.makedirs(os.path.join(prefix, "conda-meta"))
  share = os.path.join(root, "conda_prefix", "share", "aimnet")
  saved = sys.prefix, sys.platform
  try:
    sys.platform = "linux"; sys.prefix = prefix
    assert p.share_dir() == share, p.share_dir()
    # the Phenix GUI on macOS runs ${PREFIX}/python.app/Contents/.../python
    sys.platform = "darwin"
    sys.prefix = os.path.join(prefix, "python.app", "Contents")
    assert os.path.normpath(p.share_dir()) == share, p.share_dir()
    sys.platform = "win32"; sys.prefix = prefix
    assert p.share_dir() == os.path.join(
      root, "conda_prefix", "Library", "share", "aimnet"), p.share_dir()
  finally:
    sys.prefix, sys.platform = saved

def run():
  saved_env = dict(os.environ)
  saved_modules = dict(sys.modules)
  root = tempfile.mkdtemp()
  try:
    test_aimnet(root)
    test_aimnet2calc(root)
    test_cache_dir()
    test_share_dir(root)
  finally:
    shutil.rmtree(root)
    os.environ.clear(); os.environ.update(saved_env)
    for k in list(sys.modules):
      if k not in saved_modules: del sys.modules[k]
    sys.modules.update(saved_modules)

if(__name__ == "__main__"):
  run()
  print("OK")
