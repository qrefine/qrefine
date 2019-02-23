from __future__ import division
import os
import time
import shutil
import iotbx.pdb
import libtbx.load_env
from libtbx import easy_run
from scitbx.array_family import flex
from libtbx import easy_mp
import traceback

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")
qr_unit_tests_data = os.path.join(qr_unit_tests,"data_files")

def assert_folder_is_empty(prefix):
  if(len(os.listdir("."))>0):
    print "Folder is not empty: test prefix:", prefix
    print "Remove before proceeding:", os.listdir(".")
    STOP()

def setup_helix_example(pdb_name = "m00_good.pdb",
                        mtz_name = "m00_good.mtz"):
  # Read good model and compute data from it.
  xrs_good = iotbx.pdb.input(
    file_name = os.path.join(qr_unit_tests_data,
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
  if(1): print " ".join(cmd)
  return easy_run.go(" ".join(cmd))

def clean_up(prefix,mtz_name = None):
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
    'test_short_gap_capping.pdb', 'entire_qm.pdb', 'cluster_qm.pdb']
  files_to_remove2 = [
    'tst_14.pdb', 'tst_14_p1.pdb', 'tst_14_super_cell.pdb', 'super_cell.pdb',
    'tst_14_super_sphere.pdb','test_zn_his_charge.pdb']
  files_to_remove3 = ['cluster_false.pkl', 'cluster_true.pkl']
  files_to_remove4 = ['A.pdb', 'B.pdb', 'W.pdb']
  for f in files_to_remove1+files_to_remove2+files_to_remove3+files_to_remove4:
    try: os.remove(f)
    except: pass
  try:
    shutil.rmtree('ase')
    shutil.rmtree('ase_error')
  except: pass

def runner(function, prefix, disable=False):
  assert_folder_is_empty(prefix=prefix)
  rc = 0
  try:
    if(disable):
      print prefix + ": Skipped (not recommended, enable ASAP or remove!)"
    else:
      t0 = time.time()
      function(prefix = prefix)
      print prefix + ":  OK  " + "Time: %6.2f (s)" % (time.time() - t0)
  except Exception as e:
    print prefix, "FAILED", str(e)
    traceback.print_exc()
    rc=1
  clean_up(prefix)
  return rc

def run(nproc=6, only_i=None):
  try:
    only_i=int(only_i)
    nproc=1
  except: only_i=None
  t0=time.time()
  print 'Running tests on %d processors' % nproc
  def _run_test(file_name, in_separate_directory=True):
    if in_separate_directory:
      fn = file_name.split('.')[0]
      if not os.path.exists(fn):
        os.mkdir(fn)
      os.chdir(fn)
    rc = easy_run.call("cctbx.python %s"%(
      os.path.join(qr_unit_tests,file_name)))
    if in_separate_directory:
      os.chdir('..')
    return rc
  tests = [
    "tst_00.py",
    "tst_01.py",
    "tst_03.py",
    "tst_05.py",
    "tst_06.py",
    "tst_07.py",
    "tst_08.py",
    "tst_09.py",
    "tst_10.py",
    #"tst_11.py",
    "tst_12.py",
    "tst_13.py",
    "tst_14.py",
    "tst_15.py",
    "tst_16.py",
    "tst_17.py",
    "tst_18.py",
    "tst_19.py",
    "tst_20.py",
    "tst_21.py",
    "tst_22.py",
    "tst_23.py",
    "tst_24.py",
    'tst_25.py',
    'tst_26.py',
    'tst_27.py',
    'tst_28.py',
    'tst_29.py',
    'tst_30.py',
    'tst_31.py',
    'tst_32.py',
    'tst_33.py',
  ]
  failed = 0
  in_separate_directory = not(nproc==1)
  remove=[]
  for i, file_name in enumerate(tests):
    if only_i is not None:
      j = file_name[4:6]
      if int(j)!=only_i:
        remove.append(i)
        continue
      print 'adding',file_name
    tests[i]=tuple([file_name, in_separate_directory])
  if remove:
    remove.reverse()
    for r in remove: del tests[r]
    print tests
  for args, res, err_str in easy_mp.multi_core_run( _run_test,
                                                    tests,
                                                    nproc,
                                                    ):
    print '%sTotal time: %6.2f (s)' % (' '*7, time.time()-t0)
    if err_str:
      print 'Error output from %s' % args
      print err_str
      print '_'*80
    if res: failed += 1
  if failed:
    print 'Failed tests : %d' % failed
    return 1
  return 0

if(__name__ == "__main__"):
  t0 = time.time()
  args = sys.argv[1:]
  nproc=1
  if args: nproc=int(args[0])
  rc = run(nproc=nproc)
  print "Total time (all tests): %6.2f"%(time.time()-t0)
  if rc:
    assert not rc
  else:
    print "OK"
