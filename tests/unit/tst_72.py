from __future__ import absolute_import, division, print_function
import os
from qrefine.tests.unit import tst_56 as t

# Split from tst_56 (api.refine): where it writes (working folder only,
# intermediate folder, called twice, error copy of a failed cluster run).
# Helpers and fixture are in tst_56.

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  t.setup_fixture()
  t.test_no_writes_outside()
  t.test_intermediate_folder()
  t.test_called_twice()
  t.test_cluster_error_copy()

if(__name__ == "__main__"):
  run()
  print("OK")
