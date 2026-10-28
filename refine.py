"""This is the entry point into Q|R. It takes user inputs,
   and then constructs all of the objects
   needed to carry out the quantum refinement.

   - fmodel (crystallographic information)
  - calculator (composite object)
  - restraints_manager (computes energy and gradients using either qm codes or
    cctbx (standard))
  - geometry_restraints_manager (analyses geometry e.g. bond RMSDs)
  - weights (scale factors needed to scale up or down data versus restraints
    contributions)
  Then we process them by the refinement/optimization engine, driver.py:

   driver.refine(params   = params,
               fmodel     = fmodel,
               calculator = calculator_manager,
               results    = results)
   results_manager (store all reportable information, and write it out as a log,
   and also write our final pdb structure.)
   """
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import time
import pickle

from libtbx.utils import Sorry, null_out
import iotbx.pdb.hierarchy
import mmtbx.model
from cctbx.crystal.super_cell import get_siiu
from libtbx import easy_pickle
from qrefine.fragment import fragments

from qrefine import calculator
from qrefine import driver
from qrefine import restraints
from qrefine import cluster_restraints
from qrefine import conformers
from qrefine import results

from qrefine.super_cell import expand
from libtbx import Auto
from cctbx import maptbx
from scitbx.array_family import flex

def show_cc(map_data, xray_structure, log=None):
  from mmtbx.maps.correlation import five_cc
  xrs = xray_structure
  xrs.scattering_type_registry(table="electron")
  d99 = maptbx.d99(
    map = map_data, crystal_symmetry = xrs.crystal_symmetry()).result.d99
  if(log is not None): print("Resolution of map is: %6.4f" %d99, file=log)
  result = five_cc(
    map               = map_data,
    xray_structure    = xrs,
    d_min             = d99,
    compute_cc_box    = False,
    compute_cc_mask   = True,
    compute_cc_peaks  = False,
    compute_cc_volume = False).result
  if(log is not None):
    print("Map-model correlation coefficient (CC)", file=log)
    print("  CC_mask  : %6.4f" %result.cc_mask, file=log)
  return result.cc_mask

def create_fragment_manager(
      model,
      params,
      working_folder = "ase"):
  if(not params.cluster.clustering): return None
  return fragments(
    cif_objects                = model._restraint_objects,
    working_folder             = working_folder + "/",
    clustering_method          = params.cluster.clustering_method,
    maxnum_residues_in_cluster = params.cluster.maxnum_residues_in_cluster,
    charge_embedding           = params.cluster.charge_embedding,
    two_buffers                = params.cluster.two_buffers,
    distance_buffer            = params.cluster.distance_buffer,
    distance_buffer_radius     = params.cluster.distance_buffer_radius,
    pdb_hierarchy              = model.get_hierarchy(),
    qm_engine_name             = params.quantum.engine_name,
    crystal_symmetry           = model.crystal_symmetry(),
    debug                      = params.debug,
    charge_cutoff              = params.cluster.charge_cutoff,
    save_clusters              = params.cluster.save_clusters,
    select_within_radius       = params.cluster.select_within_radius,
    cap_fragments              = params.restraints == "qm")

def exchangeable_hd_pairs(pdb_hierarchy, tolerance=0.01):
  """
  Exchangeable H/D sites as [[i_seq_H_or_D], [i_seq_D_or_H]]: the pairs of
  pdb_hierarchy.exchangeable_hd_selections() (H and D in atom groups with
  different altlocs, at least one of them H/D-only), plus H and D atoms with
  equal name.strip()[1:] at the same position (within tolerance) in two atom
  groups with different altlocs that both have heavy atoms (neutron models
  where the H/D labels are those of a disordered side chain).
  """
  result = list(pdb_hierarchy.exchangeable_hd_selections())
  if(len(result) == 0):
    elements = set(e.strip() for e in pdb_hierarchy.atoms().extract_element())
    if(not ("H" in elements and "D" in elements)): return result
  known = set([(p[0][0], p[1][0]) for p in result])
  for rg in pdb_hierarchy.residue_groups():
    ags = [ag for ag in rg.atom_groups() if ag.altloc.strip() != ""]
    heavy = [ag for ag in ags if any(
      [a.element.strip() not in ("H", "D") for a in ag.atoms()])]
    for k, ag1 in enumerate(heavy):
      for ag2 in heavy[k+1:]:
        if(ag1.altloc == ag2.altloc): continue
        for a1 in ag1.atoms():
          e1 = a1.element.strip()
          if(e1 not in ("H", "D")): continue
          for a2 in ag2.atoms():
            e2 = a2.element.strip()
            if(e2 not in ("H", "D") or e1 == e2): continue
            if(a1.name.strip()[1:] != a2.name.strip()[1:]): continue
            if(a1.distance(a2) > tolerance): continue
            key = (a1.i_seq, a2.i_seq)
            if(key in known): continue
            known.add(key)
            result.append([[a1.i_seq], [a2.i_seq]])
  return result

class hd_mapper(object):
  """
  Utility function that allows to use a model with exhangeable H/D sites as
  single protonated model by re-mapping 'short' gradients into the full-size
  gradients.

  The single model drops the D of every exchangeable pair and moves its H into
  the blank atom group of its residue group (the parent atom of an exchangeable
  site is blank, so the H belongs to every conformer), then de-deuterates the
  rest. Atom order changes: self.order[k] is the full-model index of atom k of
  the single model.
  """

  def __init__(self, model):
    self.model = model
    ph = model.get_hierarchy().deep_copy()
    ph.atoms().reset_i_seq()
    self.ehd = exchangeable_hd_pairs(ph)
    self.elements = [e.strip().upper() for e in ph.atoms().extract_element()]
    self.size = model.size()
    keep = flex.bool(self.size, True)
    self.pairs = [] # (kept H, dropped D)
    for pair in self.ehd:
      i,j = pair[0][0], pair[1][0]
      if(self.elements[i] == "D"): i,j = j,i
      assert self.elements[i] == "H" and self.elements[j] == "D", [i, j]
      keep[j] = False
      self.pairs.append((i, j))
    kept_h = set([i for i, j in self.pairs])
    for a in ph.atoms():
      a.tmp = a.i_seq
    # Without copy_atoms the selection keeps the tmp of the atoms of ph (copies
    # do not); deep_copy below detaches it from ph.
    h = ph.select(keep)
    # Move the kept H of every pair into the blank atom group.
    for rg in h.residue_groups():
      moves = []
      for ag in rg.atom_groups():
        if(ag.altloc.strip() == ""): continue
        for a in ag.atoms():
          if(a.tmp in kept_h): moves.append((ag, a))
      if(len(moves) == 0): continue
      blank = [ag for ag in rg.atom_groups() if ag.altloc.strip() == ""]
      if(len(blank) > 0):
        blank = blank[0]
      else:
        blank = iotbx.pdb.hierarchy.atom_group(altloc="",
          resname=moves[0][0].resname)
        rg.insert_atom_group(0, blank)
      for ag, a in moves:
        a_new = a.detached_copy()
        a_new.tmp = a.tmp
        a_new.occ = 1.0
        ag.remove_atom(a)
        blank.append_atom(a_new)
      for ag in rg.atom_groups():
        if(ag.atoms_size() == 0): rg.remove_atom_group(ag)
    self.order = flex.size_t(h.atoms().extract_tmp_as_size_t())
    assert self.order.size() == keep.count(True)
    self.n_keep = self.order.size()
    h = h.deep_copy()
    h.atoms().reset_i_seq()
    h.de_deuterate()
    assert h.atoms_size() == self.n_keep
    h.atoms().reset_i_seq()
    h.atoms_reset_serial()
    self._model = mmtbx.model.manager(
      model_input       = None,
      pdb_hierarchy     = h,
      crystal_symmetry  = model.crystal_symmetry(),
      restraint_objects = model.get_restraint_objects(),
      log               = null_out())

  def get_single_model(self):
    return self._model

  def shrink(self, array):
    assert array.size() == self.size
    return array.select(self.order)

  def check(self, array):
    assert array.size() == self.size
    for i, j in self.pairs:
      assert array[i] == array[j]

  def average(self, array):
    for i, j in self.pairs:
      ai = flex.vec3_double([array[i],])
      aj = flex.vec3_double([array[j],])
      aa = (ai+aj)*0.5
      array[i] = aa[0]
      array[j] = aa[0]
    return array

  def expand(self, g_short):
    assert self.n_keep == g_short.size()
    g_all = flex.vec3_double(self.size, (0,0,0))
    g_all = g_all.set_selected(self.order, g_short)
    for i, j in self.pairs:
      assert g_all[i] != (0,0,0), i
      assert g_all[j] == (0,0,0), j
    for i, j in self.pairs:
      g_all[j] = g_all[i]
    for g in g_all:
      assert g != (0,0,0)
    return g_all

def create_restraints_manager(params, model, working_folder="ase", log=None):
  #
  # Conformers: one restraints manager per single-conformer model (built by
  # the routes below), gradients combined; scratch folder working_folder/<X>.
  #
  if(conformers.n_conformers(model.get_hierarchy()) > 1 and
     params.altlocs.handling == "conformers"):
    return conformers.from_conformers(
      model          = model,
      build          = lambda m, f: create_restraints_manager(params, m, f, log),
      working_folder = working_folder,
      log            = log)
  restraints_source = restraints.restraints(params = params, model = model,
    working_folder = working_folder)
  #
  # Expansion.
  #
  if(params.expansion):
    assert not params.cluster.clustering
    #
    # TMP development insert to investigate Rss BEGIN
    #
    if(params.debug_rss):
      from collections import OrderedDict
      import time
      from libtbx import easy_pickle
      result = OrderedDict()
      RSS = [10,9,8,7,6,5,4,3,2,1]
      for i, rss in enumerate(RSS):
        t0 = time.time()
        params.cluster.select_within_radius=rss
        FE = restraints.from_expansion(
          params            = params,
          restraints_source = restraints_source)
        _, g = FE.target_and_gradients(sites_cart = model.get_sites_cart())
        result[rss] = g
        print("Trying rss", rss, time.time()-t0)
      easy_pickle.dump(
        file_name="%s_gradients.pkl"%params.output_file_name_prefix, obj = result)
      STOP()
    #
    # TMP development insert to investigate Rss END
    #
    print("Restraints are from restraints.from_expansion")
    return restraints.from_expansion(
      params            = params,
      restraints_source = restraints_source)
  #
  # No expansion
  #
  else:
    if(params.cluster.clustering):
      fragment_manager = create_fragment_manager(params = params, model = model,
        working_folder = working_folder)
      return cluster_restraints.from_cluster(
        restraints_manager = restraints_source.restraints_manager,
        fragment_manager   = fragment_manager,
        parallel_params    = params.parallel)
    else:
      # restraints=cctbx clustering=false expansion=false
      return restraints_source.restraints_manager

def create_calculator(params,
                      restraints_manager,
                      model,
                      geometry_rmsd_manager,
                      fmodel=None,
                      hdm=None,
                      exclude_selection=None):
  # Define freeze selection
  fr = params.refine.freeze
  freeze_selection = None
  if fr.selection is not None:
    iselection = model.selection(string = fr.selection).iselection()
    freeze_selection = ~model.select_nearby_atoms(
      iselection = iselection,
      r_min      = fr.r_min,
      level      = int(fr.level),
      passes     = fr.passes)
  #
  if(params.refine.refine_sites):
    if(params.refine.mode == "refine"):
      assert model is not None
      return calculator.sites(
        fmodel                = fmodel,
        geometry_rmsd_manager = geometry_rmsd_manager,
        max_bond_rmsd         = params.refine.max_bond_rmsd,
        hdm                   = hdm,
        exclude_selection     = exclude_selection,
        restraints_manager    = restraints_manager,
        dump_gradients        = params.dump_gradients,
        max_shift             = params.refine.stpmax)
    else:
      # XXX
      # XXX hdm are not used XXX
      # XXX
      return calculator.sites_opt(
        restraints_manager = restraints_manager,
        model              = model,
        max_shift          = params.refine.stpmax,
        shift_eval         = params.refine.shift_evaluation,
        exclude_selection  = exclude_selection,
        freeze_selection   = freeze_selection,
        debug              = params.debug)

def set_qm_defaults(params, log):
  outl = ''
  if params.quantum.engine_name=='mopac':
    if params.quantum.method==Auto:
      params.quantum.method='PM7'
      outl += '  Setting QM method to PM7 (mopac)\n'
    if params.quantum.basis==Auto:
      params.quantum.basis=''
  if params.quantum.engine_name=='xtb':
    if params.quantum.method==Auto:
      params.quantum.method=' --gfn 2 --etemp 500 --acc 0.1 --gbsa h2o'
      print(' Default method for xtb is %s' % (
          params.quantum.method,
          ), file=log)
      params.quantum.basis = ''
  if params.quantum.engine_name=='aimnet2':
    print(file=log)
    print("Setting aimnet2 specific defaults:", file=log)
    if params.quantum.method==Auto:
      params.quantum.method="aimnet2-qr"
      print(f"  Default Aimnet2 model: {params.quantum.method}", file=log)
      try:
        import warnings
        warnings.filterwarnings('ignore', module='aimnet2calc')
        from aimnet2calc import AIMNet2ASE # intentional
      except ModuleNotFoundError as e:
        print(str(e), file=log)
        if not params.debug:
          raise Sorry("AIMNet2 does not seem to be installed or available.")
  else:
    if params.quantum.method==Auto:
      params.quantum.method='HF'
      outl += '  Setting QM method to HF\n'
    if params.quantum.basis==Auto:
      params.quantum.basis='6-31g'
      outl += '  Setting QM basis to 6-31g\n'
  if outl:
    print(' Setting QM defaults', file=log)
    print(outl, file=log)

  if params.quantum.engine_name=='mopac':
    if params.quantum.basis:
      print('  Because engine is %s basis set %s ignored' % (
        params.quantum.engine_name,
        params.quantum.basis,
        ), file=log)
      params.quantum.basis = ''
    if params.quantum.method=='hf': # default
      print(' Default method set as PM7 (mopac)', file=log)
      params.quantum.method='PM7'


def check_qm_sees_crystal_neighbours(model, params, map_data):
  """
  With restraints=qm the QM engine only sees the asymmetric unit unless the
  model is expanded (expansion=True) or fragments are built from the expansion
  (cluster.clustering=True). Refuse to silently ignore crystal neighbours that
  are within cluster.select_within_radius. Real-space refinement is exempt: the
  map box is not a physical lattice.
  """
  if(params.restraints != "qm"): return
  if(params.cluster.clustering or params.expansion): return
  if(map_data is not None): return
  radius = params.cluster.select_within_radius
  siiu, ops = get_siiu(
    pdb_hierarchy        = model.get_hierarchy(),
    crystal_symmetry     = model.crystal_symmetry(),
    select_within_radius = radius)
  if(len(ops) == 0): return
  raise Sorry(
    "Symmetry-related neighbours are within cluster.select_within_radius=%g A "
    "of the model, but with restraints=qm, clustering=False and expansion=False "
    "the QM engine only sees the asymmetric unit and ignores them. Set "
    "expansion=True (QM on the model and its neighbours) or clustering=True "
    "(fragments built with neighbours), or reduce select_within_radius if the "
    "neighbours are known to be negligible." % radius)

def run(model, fmodel, map_data, params, rst_file, prefix, log):
  set_qm_defaults(params,log)
  n_conf = conformers.n_conformers(model.get_hierarchy())
  if(n_conf > 1 and params.altlocs.handling == "native" and
     params.restraints == "qm"):
    raise Sorry("altlocs.handling=native needs restraints=cctbx: a QM engine cannot handle alternative conformations. Use altlocs.handling=conformers.")
  if(n_conf > 1 and params.altlocs.handling == "conformers"):
    params.refine.minimizer = "lbfgs"
    params.refine.gradient_only = True
    print("  alternative conformations require minimizer=lbfgs and gradient_only=True", file=log)
  check_qm_sees_crystal_neighbours(
    model = model, params = params, map_data = map_data)
  if(params.cluster.clustering):
    params.refine.gradient_only = True
    print(" params.gradient_only", params.refine.gradient_only, file=log)
  if (model.size() > params.max_atoms):
    print("Too many atoms. Can take forever or crash."*50, file=log)

  geometry_rmsd_manager = calculator.cctbx_geometry(model=model)

  #
  if(params.refine.dry_run): return
  #
  monitor = results.manager(
    model                 = model,
    geometry_rmsd_manager = geometry_rmsd_manager,
    log                   = log)
  monitor.update(fmodel = fmodel, model = model)
  if(params.rst_file is None):
    if(params.output_file_name_prefix is None):
      params.output_file_name_prefix = prefix
    if(os.path.exists(params.output_folder_name) is False):
      os.mkdir(params.output_folder_name)
    params.output_file_name_prefix = os.path.basename(params.output_file_name_prefix)
    params.rst_file = os.path.abspath(params.output_folder_name + "/" + \
      params.output_file_name_prefix + ".rst.pickle")
  if os.path.isfile(params.rst_file):
    os.remove(params.rst_file)
  start_fmodel = fmodel
  start_ph = None # is it used anywhere? I don't see where it is used!

  hdm = None
  use_hdm = len(exchangeable_hd_pairs(model.get_hierarchy())) > 0
  if use_hdm:
    print("Using hd_mapper", file=log)
    hdm = hd_mapper(model = model)

  exclude_selection = None
  if params.refine.exclude is not None:
    exclude_selection = model.selection(string = params.refine.exclude)

  if [exclude_selection, hdm].count(None)==0:
    raise Sorry("exclude_selection is incompatible with presence of H/D.")
  if params.cluster.clustering:
    if exclude_selection is not None:
      raise Sorry("exclude_selection is incompatible with clustering=True")

  if hdm is not None:
    print("Using hd_mapper again!", file=log)
    model_for_restraints = hdm.get_single_model()
  elif exclude_selection is not None:
    model_for_restraints = model.select(~exclude_selection)
  else:
    model_for_restraints = model

  restraints_manager = create_restraints_manager(
    params=params, model=model_for_restraints, log=log)

  # XXX This needs to be consolidated !

  if(map_data is not None and params.refine.mode == "refine"):
    model.geometry_statistics(use_hydrogens=False).show()
    show_cc(
      map_data        = map_data,
      xray_structure  = model.get_xray_structure(),
      log             = log)
    O = calculator.sites_real_space(
      model                 = model,
      geometry_rmsd_manager = geometry_rmsd_manager,
      max_bond_rmsd         = params.refine.max_bond_rmsd,
      map_data              = map_data,
      data_weight           = params.refine.data_weight,
      refine_cycles         = params.refine.number_of_refine_cycles,
      stpmax                = params.refine.stpmax,
      gradient_only         = params.refine.gradient_only,
      line_search           = params.refine.line_search,
      number_of_weight_search_cycles = params.refine.number_of_weight_search_cycles,
      restraints_manager    = restraints_manager,
      max_iterations        = params.refine.max_iterations_refine,
      log                   = log)
    model = O.run()
    of = open("real_space_refined.pdb", "w")
    print(model.model_as_pdb(output_cs=True), file=of)
    of.close()
    model.geometry_statistics(use_hydrogens=False).show()
    show_cc(
      map_data=map_data,
      xray_structure=model.get_xray_structure(),
      log=log)
    return
  else:
    calculator_manager = create_calculator(
      fmodel                = start_fmodel,
      model                 = model,
      geometry_rmsd_manager = geometry_rmsd_manager,
      params                = params,
      restraints_manager    = restraints_manager,
      hdm                   = hdm,
      exclude_selection     = exclude_selection)
    if(params.refine.mode == "refine"):
      #
      # Optimize H
      #
      if fmodel.f_obs().d_min()<1.2:
        import mmtbx.hydrogens
        print("Optimize H with d.o.f.", file=log)
        print ("  start: r_work=%6.4f r_free=%6.4f"%(fmodel.r_work(), fmodel.r_free()))
        mmtbx.hydrogens.fit_rotatable2(model=model, fmodel=fmodel)
        print ("  final: r_work=%6.4f r_free=%6.4f"%(fmodel.r_work(), fmodel.r_free()))
        monitor.update(fmodel = fmodel, model = model)
      #
      # PRE-OPTIMIZATION (CODITIONAL)
      #
      stats = monitor.get_stats()
      if(params.refine.refine_sites is True and
         stats.rama_z_score().whole.value is not None and (
         stats.rama_z_score().whole.value<-2.5 or
         stats.clash().score>10) and
         fmodel.f_obs().d_min()>3):
        if params.refine.exclude is not None:
          print("Initial geometry regularization skipped because of exclude_selection.", file=log)
        else:
          print("Initial geometry regularization:")
          print("  r_work: %6.4f r_free: %6.4f"%(
            fmodel.r_work(), fmodel.r_free()))
          calculator_sites_opt_pre = calculator.sites_opt(
            restraints_manager      = restraints_manager,
            model                   = model,
#            dump_gradients          = False,
            max_shift               = params.refine.stpmax,
            shift_eval              = params.refine.shift_evaluation,
            use_callback_after_step = True,
            convergence_threshold   = 0.3)
          driver.opt(
            params     = params,
            calculator = calculator_sites_opt_pre,
            monitor    = monitor)
          fmodel.xray_structure.set_sites_cart(
            sites_cart = monitor.model.get_sites_cart())
          fmodel.update_xray_structure(update_f_calc  = True)
          print("  r_work: %6.4f r_free: %6.4f"%(
            fmodel.r_work(), fmodel.r_free()))
        # PRE-OPTIMIZATION (CODITIONAL) END
      driver.refine(
        params                = params,
        fmodel                = fmodel,
        geometry_rmsd_manager = geometry_rmsd_manager,
        calculator            = calculator_manager,
        monitor               = monitor)
    else:
      driver.opt(
        params     = params,
        calculator = calculator_manager,
        monitor    = monitor)
    monitor.finalize(
      input_file_name_prefix  = prefix,
      output_file_name_prefix = params.output_file_name_prefix,
      output_folder_name      = params.output_folder_name)
