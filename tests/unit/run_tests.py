from __future__ import division
from __future__ import print_function
import os
import time
import shutil
import iotbx.pdb
import libtbx.load_env
from libtbx import easy_run
from scitbx.array_family import flex
from libtbx import easy_mp
import traceback
from libtbx.utils import Sorry

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")
qr_unit_tests_data = os.path.join(qr_unit_tests,"data_files")

def setup_helix_example(pdb_name = "m00_good.pdb",
                        mtz_name = "m00_good.mtz"):
  # Read good model and compute data from it.
  xrs_good = iotbx.pdb.input(
    file_name = os.path.join(qr_unit_tests_data,
    pdb_name )).xray_structure_simple()
  f_obs = abs(xrs_good.structure_factors(d_min=1.5).f_calc())
  r_free_flags_data = flex.bool()
  for i in range(f_obs.data().size()):
    if(i%10 == 0): r_free_flags_data.append(True)
    else:          r_free_flags_data.append(False)
  r_free_flags = f_obs.array(data = r_free_flags_data)
  mtz_dataset = f_obs.as_mtz_dataset(column_root_label="F-obs")
  mtz_dataset.add_miller_array(
    miller_array      = r_free_flags,
    column_root_label = "R-free-flags")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = mtz_name)
  xrs_poor = iotbx.pdb.input(
    file_name = os.path.join(qr_unit_tests_data,
      "m00_poor.pdb")).xray_structure_simple()
  return xrs_good,xrs_poor,f_obs,r_free_flags

def run_cmd(prefix,args,pdb_name = "m00_poor.pdb",
            mtz_name = "m00_good.mtz"):
  test_folder_name = prefix
  cmd = [
    "qr.refine",
    mtz_name,
    os.path.join(qr_unit_tests_data,pdb_name)]
  for arg in args:
    cmd.append(arg)
  cmd.append("output_folder_name=%s"%test_folder_name)
  cmd.append("> %s.log"%prefix)
  if(1): print(" ".join(cmd))
  rc = easy_run.go(" ".join(cmd))
  if rc.return_code != 0:
    rc.show_stderr()
    rc.show_stdout()
    raise SystemExit(f"A command within a test failed!")
  return rc.return_code

def runner(function, prefix, disable=False):
  prefix = "qrefine_"+prefix
  import sys
  rc = 0
  try:
    if(disable):
      print(prefix + ": Skipped (not recommended, do something about it!)")
    else:
      t0 = time.time()
      os.mkdir(prefix)
      os.chdir(prefix)
      function(prefix = prefix)
      print(prefix + ":  OK  " + "Time: %6.2f (s)" % (time.time() - t0))
  except Exception as e:
      print(prefix, "FAILED", str(e))
      print(traceback.format_exc())
      rc=1
  os.chdir('..')
  assert not rc, "%s rc: %s" % (prefix, rc)
  return rc

def run(nproc=1,
        only_i=None,
        non_mopac_only=False):
  cwd = os.getcwd()
  assert cwd.find(' ')==-1, 'test do not work in directory with a space " "'
  t0=time.time()
  print('Running tests on %d processors' % nproc)
  # Individual test runner
  def _run_test(file_name, in_separate_directory=True):
    t00=time.time()
    if in_separate_directory:
      fn = file_name.split('.')[0]
      if not os.path.exists(fn):
        os.mkdir(fn)
      os.chdir(fn)
    full_test_file_name = os.path.join(qr_unit_tests,file_name)
    print("Running test: %s in folder: %s"%(full_test_file_name,fn))
    rc = easy_run.go("qrefine.python %s"%(full_test_file_name))
    if rc.return_code != 0:
      rc.show_stderr()
      rc.show_stdout()
    if in_separate_directory:
      os.chdir('..')
    print('%sTime (this test): %6.2f (s)' % (' '*7, time.time()-t00))
    return rc.return_code
  # Collect test files
  tests = []
  for fn in os.listdir(qr_unit_tests):
    if(fn.startswith("tst_") and fn.endswith(".py")):
      i_test = fn[:].replace("tst_","").replace(".py","")
      i_test = i_test[0].replace('0','')+i_test[1]
      if(only_i is not None):
        if(only_i == i_test):
          tests.append(fn)
      else:
        tests.append(fn)
  tests.sort()
  if non_mopac_only:
    remove = []
    for i, file_name in enumerate(tests):
      f=open(os.path.join(qr_unit_tests, file_name), 'r')
      lines=f.read()
      del f
      if lines.lower().find('mopac')>-1:
        remove.append(i)
    if remove:
      remove.reverse()
      for r in remove:
        print('Removing test %s from list' % tests[r])
        del tests[r]
  print("Following tests will be executed:")
  print(" ".join(tests))
  #
  failed = 0
  in_separate_directory=True # not(nproc==1)
  for i, file_name in enumerate(tests):
    tests[i]=tuple([file_name, in_separate_directory])
  #
  for args, res, err_str in easy_mp.multi_core_run( _run_test,
                                                    tests,
                                                    nproc,
                                                    ):
    print('%sTime (total)    : %6.2f (s)' % (' '*7, time.time()-t0))
    if err_str:
      print('Error output from %s' % args)
      print(err_str)
      print('_'*80)
    if res:
      print('\n\t %s %s\n' % (args,res))
      failed += 1
  if failed:
    print('Failed tests : %d' % failed)
    return 1
  return 0

if(__name__ == "__main__"):
  t0 = time.time()
  args = sys.argv[1:]
  nproc=1
  if args: nproc=int(args[0])
  rc = run(nproc=nproc)
  print("Total time (all tests): %6.2f"%(time.time()-t0))
  if rc:
    assert not rc
  else:
    print("OK")
