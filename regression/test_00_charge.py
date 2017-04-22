import os, sys
import time
from StringIO import StringIO
import iotbx
import test_data
from iotbx import pdb
import libtbx.load_env
from libtbx import easy_run
from libtbx import easy_mp
import qrefine.core.finalise
import qrefine.core.charges as charges
import qrefine.core.completion
from qrefine.core.utils import hierarchy_utils

qr_repo_parent = libtbx.env.find_in_repositories("qrefine")

#is this a unit test?

def test_qxyz_non_zero():
  def _check_non_zero_charge(filename):
    for line in open(filename, 'rb').readlines():
      tmp = line.split()
      if len(tmp)<2: continue
      assert float(tmp[0])!=0, 'no partial charge %s' % line
  for residue in ['pro',
                  'gly',
                  ]:
    tf='%s.pdb' % residue
    f=file(tf, "wb")
    f.write(test_data.pdbs[residue])
    f.close()
    pdb_inp = pdb.input(tf)
    hierarchy = pdb_inp.construct_hierarchy()
    charges.write_pdb_hierarchy_qxyz_file(hierarchy,
                                               'test_%s.dat' % residue,
                                             )
    _check_non_zero_charge('test_%s.dat' % residue)

def test_qxyz_xyzq():
  tf='water.pdb'
  f=file(tf, "wb")
  f.write(test_data.pdbs["water"])
  f.close()
  pdb_inp = pdb.input(tf)
  hierarchy = pdb_inp.construct_hierarchy()
  if  os.path.exists('test_water.dat'): os.remove('test_water.dat')
  charges.write_pdb_hierarchy_qxyz_file(hierarchy,
                                             'test_water.dat',
                                             )
  assert not os.path.exists('test_water.dat')
  charges.write_pdb_hierarchy_qxyz_file(hierarchy,
                                             'test_water.dat',
                                             exclude_water=False,
                                             )
  tst_str = '''\
3  
  
-0.408  -0.21  0.0  -0.296    
0.204  0.733  0.0  -0.296    
0.204  -0.524  0.0  0.593'''
  lines = open('test_water.dat', 'rb').read()
  assert lines.strip()==tst_str, '%s %s' % (tst_str, lines)
  os.remove('test_water.dat')
  charges.write_pdb_hierarchy_xyzq_file(hierarchy,
                                             'test_water.dat',
                                             exclude_water=False,
                                             )
  tst_str = '''\
-0.21  0.0  -0.296  -0.408    
0.733  0.0  -0.296  0.204    
-0.524  0.0  0.593  0.204'''
  lines = open('test_water.dat', 'rb').read()
  assert lines.strip()==tst_str, '%s'% (lines)
  os.remove('test_water.dat')

def test_1yjp_charge():
  tf='%s/qr-core/regression/datasets/1yjp.pdb' % qr_repo_parent
  try: os.symlink(tf, os.path.basename(tf))
  except: pass
  tf = os.path.basename(tf)
  pdb_inp = pdb.input(tf)
  hierarchy = pdb_inp.construct_hierarchy()
  try:
    charge = charges.calculate_pdb_hierarchy_charge(hierarchy)
    assert 0
  except Exception, e:
    assert e.message.find('no hydrogens')>-1
  cmd = 'iotbx.python %s/qr-core/finalise.py %s' % (qr_repo_parent, tf)
  easy_run.call(cmd)
  pdb_inp = pdb.input(tf.replace('.pdb', '_complete.pdb'))
  hierarchy = pdb_inp.construct_hierarchy()
  charge = charges.calculate_pdb_hierarchy_charge(hierarchy, verbose=1)
  assert charge==0, 'charge of 1yjp should be zero not %s' % charge

def test_terminal_charge(residue, charge=0):
  # must run after other PRO
  tf = '%s_terminal_complete.pdb' % residue
  ppf = hierarchy_utils .get_processed_pdb(pdb_filename=tf)
  inter_residue_bonds = charges.get_inter_residue_bonds(ppf)
  # should the hierarchy come from ppf???
  pdb_inp = iotbx.pdb.input(tf)
  hierarchy = pdb_inp.construct_hierarchy()
  hetero_charges = charges.get_hetero_charges(pdb_inp)
  if not hetero_charges:
    # some defaults
    hetero_charges = charges.default_ion_charges
  total_charge = charges.calculate_pdb_hierarchy_charge(
    hierarchy,
    hetero_charges=hetero_charges,
    inter_residue_bonds=inter_residue_bonds,
    verbose=True,
  )
  assert total_charge==charge, "total_charge: %d, charge:%d"%(total_charge,charge)

def test_PRO_terminal_charge():
  test_terminal_charge('PRO')

def test_GLY_terminal_charge():
  test_terminal_charge('GLY')

def test_helix():
  tf = 'helix.pdb'
  f=file(tf, "wb")
  f.write(test_data.pdbs["helix"])
  f.close()
  pdb_inp=pdb.input(tf)
  hierarchy = pdb_inp.construct_hierarchy()
  charge = charges.calculate_pdb_hierarchy_charge(hierarchy, verbose=1)
  assert charge==0, 'charge of helix should be zero not %s' % charge
  cmd = 'iotbx.python %s/qr-core/finalise.py %s' % (qr_repo_parent, tf)
  easy_run.call(cmd)
  pdb_inp = pdb.input(tf.replace('.pdb', '_complete.pdb'))
  hierarchy = pdb_inp.construct_hierarchy()
  must_find = ['H1', 'H2', 'H3', 'OXT']
  for atom in hierarchy.atoms():
    if atom.name.strip() in must_find:
      must_find.remove(atom.name.strip())
  assert not must_find
  pdb_inp=pdb.input(tf.replace('.pdb', '_complete.pdb'))
  hierarchy = pdb_inp.construct_hierarchy()
  charge = charges.calculate_pdb_hierarchy_charge(hierarchy, verbose=1)
  assert charge==1, 'charge of helix should be one not %s' % charge

def test_charge_for_charmm_pdbs(only_i=None):
  from qrefine.core.charges import calculate_pdb_hierarchy_charge
  charge_dict = {'3kyi': -12,
                 '2oy0': 16,
                 '1y1l': -8,
                 '3dtj': 0,
                 '3tz9': -10,
                 '4rnf': -13,
                 '2jee': -32,
                 '4k2r': -1,
                 '5d12': -11,
                 '1il5': 0,
                 '1va7': -8,
                 '2oeq': -5,
                 '4drw': -16,
                 '4xa1': -45,
                 '2ghj': 46,
                 '1ok9': -14, 
		 '3nak': 0, 
		 '2x10': -22, 
		 '3oe9': 27, 
		 '4fsx': -42, 
		 '4p7h': -19,
		 '5diz': -25, 
		 '3uds': -8, 
		 '3uj4': 0, 
		 '4ctd': -44,
                 '1byz': -4, 
                 '1lzt': 8, 
		 '1vfy': -1, 
		 '1m24': -2, 
		 '1i07': 4, 
		 '1opd': -6, 
  		 '1vbw': 8, 
		 '1rfs': -4, 
		 '1ly2': 7, 
		 '1a7y': 0, 
		 '1v7s': 8,
	         '2wpz': 0, 
		 '2akf': -3,
		 '4lzt': 8,
		 '3ovj': 4,
		 '4lzl': -5,
		 '4w71': 0,
		 '4uiv': 4,
		 '4rp6': -1,
		 '3u29': 0,
		 '2omq': -4,
		 '2ol9': 0,
		 '4xfo': 0,
		 '2xmu': -10,
		 '2xmt': -10,
		 '4uit': 4,
		 '4uiu': 4,
		 '4onk': 0,
		 '2f2n': 8,
		 '4uiw': 4,
		 '5cgo': 6,
		 '2y3j': 0,
		 '2i1u': -1,
		 '4w67': 0, 
		 '3osm': 9,
		 '4wxt': -4,
		 '3o2h': -1,
		 '2f30': 8,
		 '4w5y': 0,
		 '5e5v': 0,
		 '2ona': 0, 
		 '4itk': -9,
		 '5e61': 0,
		 '4kdw': -17,
		 '4z0w': -2,
		 '4uby': 0,
		 '2w9r': -1,
		 '3t4f': 0,
		 '3pzz': 0,
		 '2y3k': 0}
  pdb_dir = '%s/qr-core/regression/charmm' % qr_repo_parent
  pdb_files = os.listdir(pdb_dir)
  for i, pdb_file in enumerate(pdb_files):
    if only_i is not None and i!=only_i: continue
    if pdb_file[:-4] not in charge_dict: continue
    if pdb_file.endswith(".pdb"):
      if pdb_file in ['2i1u.pdb',
                      '4wxt.pdb',
      ]:
        # disuphide bridge = in 2i1u 2.94 Phenix says yes, Charmm says no
        continue
      pdb_file_path = os.path.join(pdb_dir, pdb_file)
      charge = charges.get_total_charge_from_pdb(pdb_file_path)
      assert charge==charge_dict[pdb_file[:-4]], \
        '%s charge is %d, charmm charge is %d,  no matchy matchy' % (
          pdb_file[:-4],
          charge,
          charge_dict[pdb_file[:-4]],
        )

def test_charge_of_neutral_terminal():
  charge = 0
  tf = "%s/qr-core/regression/babel/clusters/neutral_nterminal.pdb" % qr_repo_parent
  charge_neutral_nterminal = charges.get_total_charge_from_pdb(tf)
  assert charge == charge_neutral_nterminal, 'no match %s %s' % (
    charge,
    charge_neutral_nterminal,
    )
  tf = "%s/qr-core/regression/babel/clusters/neutral_cterminal.pdb" % qr_repo_parent
  charge_neutral_cterminal = charges.get_total_charge_from_pdb(tf)
  assert charge == charge_neutral_cterminal, 'no match %s %s' % (
    charge,
    charge_neutral_cterminal,
    )


def run(prefix = "tst_reg_01_charge"):
  """
  Exercise structure preparation including charge, capping, completion
  """
  tests = [
    [test_charge_of_neutral_terminal, 1],
    [test_qxyz_non_zero, 1],
    [test_helix, 1],
    [test_qxyz_xyzq, 1],
    [test_1yjp_charge, 1],
    [test_charge_for_charmm_pdbs, 80],
    # need files from alt loc test
    [test_GLY_terminal_charge, 1],
    [test_PRO_terminal_charge, 1],
    ]
  def get_test(i, j):
    func = tests[i][0]
    if j is not None: return func(j)
    else: return func()
  ####
  argss = []
  for i, (func, j) in enumerate(tests):
    if j==1: argss.append((i,None))
    else:
      for p in range(j):
        argss.append((i,p))
  # for testing the regression
  if 0:
    argss = [
      [9, 60],
      ]
  #
  passed=0
  failed=0
  print dir(easy_mp)
  for args, res, errstr in easy_mp.multi_core_run(get_test, argss, 6):
    if errstr:
      print '-'*80
      print args
      print 'RESULT - ERROR   : %s %s' % (tests[args[0]][0].func_name, args)
      print errstr
      print '-'*80
      failed+=1
    else:
      print 'RESULT - SUCCESS : %s %s' % (tests[args[0]][0].func_name, args)
      passed+=1
  print '\n\tpassed : %d\n\tfailed : %s' % (passed, failed)

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print "Time: %6.2f"%(time.time()-t0)
  print "OK"
