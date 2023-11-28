from __future__ import division
from __future__ import absolute_import

import os
from qrefine.tests.unit import run_tests
import libtbx.load_env
import iotbx.pdb
import mmtbx.restraints
from libtbx.test_utils import approx_equal
from qrefine import restraints
from scitbx.array_family import flex
import mmtbx.model
from qrefine.command_line import granalyse

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")

def run(prefix):
  """
  CCTBX gradients in presence of altlocs. method = "average"
  
  XXX No assertions yet.
  """
  method = "average"
  for file_name in ["gly2_1.pdb", "gly2_2.pdb", "h_altconf_complete.pdb",
                    "h_altconf_2_complete.pdb", "altlocs2.pdb", "altlocs.pdb"]:


    file_name = os.path.join(qr_unit_tests,"data_files",file_name)
    print(file_name)
    pi = iotbx.pdb.input(file_name = file_name)
    
    model = restraints.model_from_hierarchy(
      pdb_hierarchy    = pi.construct_hierarchy(), 
      crystal_symmetry = pi.crystal_symmetry())
      
    rm = restraints.from_cctbx(
      restraints_manager = model.get_restraints_manager())
    _, g0 = rm.target_and_gradients(sites_cart = model.get_sites_cart())

    # This uses primitive ad-hoc code to get gradients
    ph = model.get_hierarchy()
    cs = model.crystal_symmetry()
    # Gradients from CCTBX (directly)
    g1 = restraints.get_cctbx_gradients(ph = ph, cs = cs).gradients
    # Gradients from CCTBX (via decomposing hierarchy into concormers A, B, blanc)
    g2 = restraints.from_cctbx_altlocs(ph = ph, cs = cs, option=1, method=method)
    g3 = restraints.from_cctbx_altlocs(ph = ph, cs = cs, option=2, method=method)

    # This does the same but using the standard code (expansion, as qr.refine).
    from qrefine import qr
    params = qr.get_default_params()
    params.restraints="cctbx"
    params.cluster.clustering=False

    _, gX = restraints.from_altlocs2(model = model, method=method,
      params = params).target_and_gradients(sites_cart = ph.atoms().extract_xyz())
     
    #
    #assert approx_equal(g1, gX)
    #assert approx_equal(g0, g1)
    #assert approx_equal(g1, g2)
    #assert approx_equal(g2, g3)
    
    print( "(g1, gX):", max(granalyse.get_grad_wdelta(ref=g1, g=gX)) )
    print( "(g0, g1):", max(granalyse.get_grad_wdelta(ref=g0, g=g1)) )
    print( "(g1, g2):", max(granalyse.get_grad_wdelta(ref=g1, g=g2)) )
    print( "(g2, g3):", max(granalyse.get_grad_wdelta(ref=g2, g=g3)) )
    

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
