from __future__ import division
import os
import time
import shutil
import iotbx.pdb
import libtbx.load_env
from libtbx import easy_run
from scitbx.array_family import flex

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests/unit")

def setup_helix_example(pdb_name = "data_files/m00_good.pdb",
                        mtz_name = "data_files/m00_good.mtz"):
  # Read good model and compute data from it.
  xrs_good = iotbx.pdb.input(
    file_name = os.path.join(qr_unit_tests,
    pdb_name )).xray_structure_simple()
  f_obs = abs(xrs_good.structure_factors(d_min=1.5).f_calc())
  r_free_flags_data = flex.bool()
  for i in xrange(f_obs.data().size()):
    if(i%10 == 0): r_free_flags_data.append(True)
    else:          r_free_flags_data.append(False)
  r_free_flags = f_obs.array(data = r_free_flags_data)
  mtz_dataset = f_obs.as_mtz_dataset(column_root_label="F-obs")
  mtz_dataset.add_miller_array(
    miller_array      = r_free_flags,
    column_root_label = "R-free-flags")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = os.path.join(qr_unit_tests,
    mtz_name))
  # Read poor model
  xrs_poor = iotbx.pdb.input(
    file_name = os.path.join(qr_unit_tests,
      "data_files/m00_poor.pdb")).xray_structure_simple()
  return xrs_good,xrs_poor,f_obs,r_free_flags

def run_cmd(prefix,args,pdb_name = "data_files/m00_poor.pdb",
            mtz_name = "data_files/m00_good.mtz"):
  test_folder_name = "./%s/"%prefix
  cmd = ["cctbx.python",
        os.path.join(qrefine,"qr.py"),
        os.path.join(qr_unit_tests,mtz_name),
        os.path.join(qr_unit_tests,pdb_name)]
  for arg in args:
      cmd.append(arg)
  cmd.append("output_folder_name=%s"%test_folder_name)
  cmd.append("> %s.log"%prefix)
  #print " ".join(cmd)
  return easy_run.go(" ".join(cmd))

def clean_up(prefix,mtz_name = "data_files/m00_good.mtz"):
  test_folder_name = "./%s/"%prefix
  print "Removing:", test_folder_name
  if(os.path.exists(test_folder_name)):
    shutil.rmtree(test_folder_name)
  os.remove(os.path.join(qr_unit_tests,mtz_name))
  os.remove("%s.log"%prefix)

def run():
  tests = [
    "tst_00.py",
    "tst_01.py",
    "tst_02.py",
    "tst_03.py",
    "tst_04.py",
    "tst_05.py",
    "tst_06.py",
    "tst_07.py",
    "tst_08.py",
    "tst_09.py",
    "tst_10.py",
    "tst_11.py",
    "tst_12.py",
    "tst_13.py"
  ]
  for file_name in tests:
    print "Running test:", file_name
    easy_run.call("cctbx.python %s"%(
      os.path.join(qr_unit_tests,file_name)))

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print "Total time (all tests): %6.2f"%(time.time()-t0)
  print "OK"
