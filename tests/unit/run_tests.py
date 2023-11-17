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

def assert_folder_is_empty(prefix):
  if(len(os.listdir("."))>0):
    print('-'*80)
    print("Folder is not empty: test prefix:", prefix)
    print("Remove before proceeding:", os.listdir("."))
    print('-'*80)
    raise Sorry("Folder is not empty: test prefix:", prefix)

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
  return easy_run.go(" ".join(cmd))

def clean_up(prefix, mtz_name = None):
  test_folder_name = prefix
  if(os.path.exists(test_folder_name)):
    shutil.rmtree(test_folder_name)
  if(mtz_name is not None): os.remove(mtz_name)
  try: os.remove("%s.log"%prefix)
  except: pass
  try: os.remove("m00_good.mtz")
  except: pass
  files_to_remove1 = [
    'c_terminal_capping.pdb', 'c_terminal_capping_capping.pdb', 'cluster.xml',
    'helix.pdb', 'helix_complete.pdb', 'helix_readyset_input.pdb', 'qmmm.xml',
    'test_10_capping.pdb', 'test_cys_hg_capping.pdb',
    'test_point_charges.pdb', 'q.xyz', "expansion.pdb",
    'test_cys_hg_capping_capping.pdb', 'test_original_pdb.pdb',
    'test_original_pdb_capping.pdb', 'test_short_gap.pdb',
    'test_short_gap_capping.pdb', 'entire_qm.pdb', 'cluster_qm.pdb',
    '%s.pdb'%prefix, '%s.log'%prefix, '%s_complete.pdb'%prefix,
    '%s_readyset_input.pdb'%prefix]
  files_to_remove2 = [
    'tst_14.pdb', 'tst_14_p1.pdb', 'tst_14_super_cell.pdb', 'super_cell.pdb',
    'tst_14_super_sphere.pdb','test_zn_his_charge.pdb']
  files_to_remove3 = ['cluster_false.pkl', 'cluster_true.pkl','1-20.npy']
  files_to_remove4 = ['A.pdb', 'B.pdb', 'W.pdb']
  for f in files_to_remove1+files_to_remove2+files_to_remove3+files_to_remove4:
    try: os.remove(f)
    except: pass
  for folder in ['ase', 'ase_error', 'pdb', '1-20','frag_pdbs']:
    try:
      shutil.rmtree(folder)
    except: pass

def runner(function, prefix, disable=False):
  import sys
  assert_folder_is_empty(prefix=prefix)
  rc = 0
  try:
    if(disable):
      print(prefix + ": Skipped (not recommended, do something about it!)")
    else:
      t0 = time.time()
      function(prefix = prefix)
      print(prefix + ":  OK  " + "Time: %6.2f (s)" % (time.time() - t0))
      # Turning this off to get some output from CI pipeline
      if os.environ['CLEANUP_TESTS'] == 'FALSE':
          print("not cleaning up tests")
      else:  
        clean_up(prefix)
  except Exception as e:
      print(prefix, "FAILED", str(e))
      exc_type, exc_value, exc_traceback = sys.exc_info() 
      traceback_template = ''' ** qrefine exception handler: **
      %(type)s => File "%(filename)s" \n line %(lineno)s, in %(name)s: \n %(message)s
        \n'''
      traceback_details = {
                        'filename': exc_traceback.tb_frame.f_code.co_filename,
                        'lineno'  : exc_traceback.tb_lineno,
                        'name'    : exc_traceback.tb_frame.f_code.co_name,
                        'type'    : exc_type.__name__,
                        'message' : exc_value.args[0]
                      }
      del(exc_type, exc_value, exc_traceback)
      print(traceback.format_exc())
      print(traceback_template % traceback_details)
    # traceback.print_exc()
      rc=1
  # clean_up(prefix)
  assert not rc, "%s rc: %s" % (prefix, rc)
  return rc

def run(nproc=6,
        only_i=None,
        non_mopac_only=False):
  cwd = os.getcwd()
  assert cwd.find(' ')==-1, 'test do not work in directory with a space " "'
  t0=time.time()
  print('Running tests on %d processors' % nproc)
  def _run_test(file_name, in_separate_directory=True):
    if in_separate_directory:
      fn = file_name.split('.')[0]
      if not os.path.exists(fn):
        os.mkdir(fn)
      os.chdir(fn)
    rc = easy_run.call("qrefine.python %s"%(
      os.path.join(qr_unit_tests,file_name)))
    if in_separate_directory:
      os.chdir('..')
    return rc
  # Collect test files
  tests = []
  for fn in os.listdir(qr_unit_tests):
    if(fn.startswith("tst_") and fn.endswith(".py")):
      # i_test = int(fn[:].replace("tst_","").replace(".py",""))
      i_test = fn[:].replace("tst_","").replace(".py","")
      i_test = i_test[0].replace('0','')+i_test[1]
      # print(i_test,only_i)
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
  for args, res, err_str in easy_mp.multi_core_run( _run_test,
                                                    tests,
                                                    nproc,
                                                    ):
    print('%sTotal time: %6.2f (s)' % (' '*7, time.time()-t0))
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
