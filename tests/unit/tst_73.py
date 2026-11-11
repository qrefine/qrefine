from __future__ import absolute_import, division, print_function
import os
from qrefine.tests.unit import tst_56 as t

# Split from tst_56 (api.refine): no files with working_folder=None, a removed
# temporary folder for MOPAC, opt and real-space modes, raw params, mlhl.
# Helpers and fixture are in tst_56.

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  t.setup_fixture()
  t.test_no_files_at_all()
  t.test_no_files_at_all_clustering()
  t.test_mopac_tempdir_removed()
  t.test_opt_mode()
  t.test_real_space_mode()
  t.test_raw_params_refused()
  t.test_auto_target_mlhl()

if(__name__ == "__main__"):
  run()
  print("OK")
