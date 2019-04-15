from __future__ import division

import time, os
import iotbx.pdb
from libtbx.test_utils import approx_equal
import iotbx.pdb
from scitbx.array_family import flex
import string
from cctbx import uctbx
from cctbx import crystal
from libtbx.utils import null_out
import mmtbx
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.restraints
from mmtbx import monomer_library
import libtbx.load_env
import run_tests

mon_lib_srv = mmtbx.monomer_library.server.server()
ener_lib    = mmtbx.monomer_library.server.ener_lib()

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine,"tests","unit","data_files")

def get_grm(ph, cs):
  # XXX Nth copy-paste
  params = mmtbx.monomer_library.pdb_interpretation.master_params.extract()
  params.use_neutron_distances = True
  params.restraints_library.cdl = False
  params.sort_atoms = False
  processed_pdb_file = mmtbx.monomer_library.pdb_interpretation.process(
    mon_lib_srv              = mon_lib_srv,
    ener_lib                 = ener_lib,
    params                   = params,
    crystal_symmetry         = cs,
    pdb_inp                  = ph.as_pdb_input(),
    strict_conflict_handling = False,
    force_symmetry           = True,
    log                      = null_out())
  xrs = processed_pdb_file.xray_structure()
  sctr_keys = xrs.scattering_type_registry().type_count_dict().keys()
  has_hd = "H" in sctr_keys or "D" in sctr_keys
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies                = False,
    assume_hydrogens_all_missing = not has_hd,
    plain_pairs_radius           = 5.0)
  return mmtbx.restraints.manager(
     geometry = geometry, normalization = False)

def get_grads(sel_f_str, sel_buffer_str, file_name):
  pdb_inp = iotbx.pdb.input(file_name=file_name)
  ph = pdb_inp.construct_hierarchy()
  grm = get_grm(ph=ph, cs=pdb_inp.crystal_symmetry())
  es = grm.energies_sites(sites_cart=ph.atoms().extract_xyz(),
    compute_gradients=True)
  asc = ph.atom_selection_cache()
  sel_f = asc.selection(sel_f_str)
  sel_buffer = asc.selection(sel_buffer_str)
  return es.gradients, sel_f, sel_buffer

def run3(prefix):
  file_name=os.path.join(qr_unit_tests_data,"h_altconf_2.pdb")

  s_b_W1_A_str = "altloc A or resseq 95:97"
  s_f_W1_A_str = "resseq 87:93 or altloc A or resseq 95:97"

  s_b_W1_B_str = "altloc B or resseq 95:97"
  s_f_W1_B_str = "resseq 87:93 or altloc B or resseq 95:97"

  s_b_W1_AB_str = "resseq 95:97"
  s_f_W1_AB_str = "resseq 87:93 or resseq 95:97"

  s_b_W2_str = "altloc B or resseq 91:93"
  s_f_W2_str = "resseq 95:99 or altloc B or resseq 91:93"

  s_b_A_str = "resseq 90:93 or resseq 95:96"
  s_f_A_str = "altloc A or resseq 90:93 or resseq 95:96"

  s_b_B_str = "resseq 90:93 or resseq 95:96"
  s_f_B_str = "altloc B or resseq 90:93 or resseq 95:96"

  g    , junk1, junk2           = get_grads(file_name=file_name, sel_f_str = "all",     sel_buffer_str = "not all")
  gA_f , s_A_f,s_A_b            = get_grads(file_name=file_name, sel_f_str =s_f_A_str,  sel_buffer_str = s_b_A_str)
  gB_f , s_B_f,s_B_b            = get_grads(file_name=file_name, sel_f_str =s_f_B_str,  sel_buffer_str = s_b_B_str)
  gW1_A_f, s_W1_A_f,s_W1_A_b    = get_grads(file_name=file_name, sel_f_str =s_f_W1_A_str, sel_buffer_str = s_b_W1_A_str)
  gW1_B_f, s_W1_B_f,s_W1_B_b    = get_grads(file_name=file_name, sel_f_str =s_f_W1_B_str, sel_buffer_str = s_b_W1_B_str)
  gW1_AB_f, s_W1_AB_f,s_W1_AB_b = get_grads(file_name=file_name, sel_f_str =s_f_W1_AB_str, sel_buffer_str = s_b_W1_AB_str)
  gW2_f, s_W2_f,s_W2_b          = get_grads(file_name=file_name, sel_f_str =s_f_W2_str, sel_buffer_str = s_b_W2_str)

  g = flex.vec3_double(g.size())
  g_W1_A = g.set_selected(s_W1_A_f, gW1_A_f)
  g_W1_A = g_W1_A.set_selected(s_W1_A_b, [0,0,0])

  g = flex.vec3_double(g.size())
  g_W1_B = g.set_selected(s_W1_B_f, gW1_B_f)
  g_W1_B = g_W1_B.set_selected(s_W1_B_b, [0,0,0])

  g = flex.vec3_double(g.size())
  g_W1_AB = g.set_selected(s_W1_AB_f, gW1_AB_f)
  g_W1_AB = g_W1_AB.set_selected(s_W1_AB_b, [0,0,0])

  g = flex.vec3_double(g.size())
  g_W2 = g.set_selected(s_W2_f, gW2_f)
  g_W2 = g_W2.set_selected(s_W2_b,  [0,0,0])

  g = flex.vec3_double(g.size())
  g_A = g.set_selected(s_A_f, gA_f)
  g_A = g_A.set_selected(s_A_b,[0,0,0])

  g = flex.vec3_double(g.size())
  g_B = g.set_selected(s_B_f, gB_f)
  g_B = g_B.set_selected(s_B_b,[0,0,0])

  g1, junk1,junk2 = get_grads(file_name=file_name, sel_f_str = "all", sel_buffer_str = "not all")
  g2 = g_A + g_B + g_W1_A + g_W1_B + g_W2  -g_W1_AB

  diff = g1-g2
  assert approx_equal(diff.max(), [0,0,0])
  if(0):
    print diff.max()
    for d in diff:
      print d

def run2(prefix):
  file_name=os.path.join(qr_unit_tests_data,"h_altconf.pdb")

  s_b_W1_str = "altloc A or resseq 95:97"
  s_f_W1_str = "resseq 87:93 or altloc A or resseq 95:97"

  s_b_W2_str = "altloc B or resseq 91:93"
  s_f_W2_str = "resseq 95:99 or altloc B or resseq 91:93"

  s_b_A_str = "resseq 90:93 or resseq 95:96"
  s_f_A_str = "altloc A or resseq 90:93 or resseq 95:96"

  s_b_B_str = "resseq 90:93 or resseq 95:96"
  s_f_B_str = "altloc B or resseq 90:93 or resseq 95:96"

  g    , junk1, junk2  = get_grads(file_name=file_name, sel_f_str = "all",     sel_buffer_str = "not all")
  gA_f , s_A_f, s_A_b  = get_grads(file_name=file_name, sel_f_str =s_f_A_str,  sel_buffer_str = s_b_A_str)
  gB_f , s_B_f, s_B_b  = get_grads(file_name=file_name, sel_f_str =s_f_B_str,  sel_buffer_str = s_b_B_str)
  gW1_f, s_W1_f,s_W1_b = get_grads(file_name=file_name, sel_f_str =s_f_W1_str, sel_buffer_str = s_b_W1_str)
  gW2_f, s_W2_f,s_W2_b = get_grads(file_name=file_name, sel_f_str =s_f_W2_str, sel_buffer_str = s_b_W2_str)

  g = flex.vec3_double(g.size())
  g_W1 = g.set_selected(s_W1_f, gW1_f)
  g_W1 = g_W1.set_selected(s_W1_b, [0,0,0])

  g = flex.vec3_double(g.size())
  g_W2 = g.set_selected(s_W2_f, gW2_f)
  g_W2 = g_W2.set_selected(s_W2_b,  [0,0,0])

  g = flex.vec3_double(g.size())
  g_A = g.set_selected(s_A_f, gA_f)
  g_A = g_A.set_selected(s_A_b,[0,0,0])

  g = flex.vec3_double(g.size())
  g_B = g.set_selected(s_B_f, gB_f)
  g_B = g_B.set_selected(s_B_b,[0,0,0])

  g1, junk1,junk2 = get_grads(file_name=file_name,sel_f_str = "all", sel_buffer_str = "not all")
  g2 = g_A + g_B + g_W1 + g_W2

  diff = g1-g2
  assert approx_equal(diff.max(), [0,0,0])
  if(0):
    print diff.max()
    for d in diff:
      print d

def run1(prefix):
  pdb_inp = iotbx.pdb.input(file_name=os.path.join(qr_unit_tests_data,"m.pdb"))
  ph = pdb_inp.construct_hierarchy()
  grm = get_grm(ph=ph, cs=pdb_inp.crystal_symmetry())
  es = grm.energies_sites(sites_cart=ph.atoms().extract_xyz(),
    compute_gradients=True)
  #
  asc = ph.atom_selection_cache()
  sA = asc.selection("altloc A or altloc ' '")
  if 0: print sA.count(True)
  phA = ph.select(sA)
  phA.write_pdb_file("A.pdb")

  sB = asc.selection("altloc B or altloc ' '")
  if 0: print sB.count(True)
  phB = ph.select(sB)
  phB.write_pdb_file("B.pdb")

  sW  = asc.selection("altloc ' '")
  if 0: print sW.count(True)
  phW = ph.select(sW)
  phW.write_pdb_file("W.pdb")
  #
  grmA = grm.select(sA)
  esA = grmA.energies_sites(sites_cart=phA.atoms().extract_xyz(),
    compute_gradients=True)

  grmB = grm.select(sB)
  esB = grmB.energies_sites(sites_cart=phB.atoms().extract_xyz(),
    compute_gradients=True)

  grmW = grm.select(sW)
  esW = grmW.energies_sites(sites_cart=phW.atoms().extract_xyz(),
    compute_gradients=True)
  #
  if 0:
    print list(phW.atoms())[0].name
    print esW.gradients[0]
    print

    print list(ph.atoms())[0].name
    print es.gradients[0]
    print esA.gradients[0]
    print esB.gradients[0]
    print
    print list(phA.atoms())[5].name, list(ph.atoms())[5].name
    print es.gradients[5]
    print esA.gradients[5]
    print
    print list(phB.atoms())[5].name, list(ph.atoms())[12].name
    print es.gradients[12]
    print esB.gradients[5]
    print
    print "-------"
  result = list(flex.double(esA.gradients[0])+
                flex.double(esB.gradients[0])-
                flex.double(esW.gradients[0]))
  assert approx_equal(result, es.gradients[0])

if(__name__ == '__main__'):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run1, prefix=prefix, disable=False)
  run_tests.runner(function=run2, prefix=prefix, disable=False)
  run_tests.runner(function=run3, prefix=prefix, disable=False)
