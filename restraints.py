from __future__ import division
from __future__ import absolute_import

import os
import ase.units as ase_units
import mmtbx.restraints
from libtbx.utils import Sorry
from .charges import charges_class
from scitbx.array_family import flex
from .plugin.ase.mopac_qr import Mopac
from .plugin.ase.pyscf_qr import Pyscf
from .plugin.ase.terachem_qr import TeraChem
from .plugin.ase.turbomole_qr import Turbomole
from .plugin.ase.orca_qr import Orca
from .plugin.ase.gaussian_qr import Gaussian
from .plugin.ase.xtb_qr import GFNxTB
from .plugin.ase.server_qr import RestAPICalculator
from .plugin.tools import qr_tools
from libtbx import group_args
import math
from qrefine.super_cell import expand
import qrefine.completion as model_completion
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from collections import OrderedDict
from libtbx import adopt_init_args
import time
from libtbx import easy_pickle

def model_from_hierarchy(pdb_hierarchy, crystal_symmetry, cif_objects=None):
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = True
  params.pdb_interpretation.restraints_library.cdl = False
  params.pdb_interpretation.sort_atoms = False
  params.pdb_interpretation.flip_symmetric_amino_acids = False
  params.pdb_interpretation.correct_hydrogens=False
  model = mmtbx.model.manager(
    model_input       = pdb_hierarchy.as_pdb_input(),
    crystal_symmetry  = crystal_symmetry,
    restraint_objects = cif_objects,
    log               = null_out())
  model.process(make_restraints=True, grm_normalization=False,
    pdb_interpretation_params = params,)
  return model

class restraints(object):
  """
  Create CCTBX or QM restraints manager.
  """
  def __init__(self, params, model):
    self.params = params
    self.model = model
    self.cif_objects      = model.get_restraint_objects()
    self.pdb_hierarchy    = model.get_hierarchy()
    self.crystal_symmetry = model.crystal_symmetry()
    self.pi_params        = model.get_current_pdb_interpretation_params()
    self.restraints_manager = None
    self.update(
      pdb_hierarchy    = model.get_hierarchy(),
      crystal_symmetry = model.crystal_symmetry())

  def source_of_restraints_qm(self):
    return self.params.restraints == "qm"

  def update(self, pdb_hierarchy, crystal_symmetry):
    #
    if self.restraints_manager is not None:
      if(not self.source_of_restraints_qm()):
        size = self.restraints_manager.geometry_restraints_manager.geometry.\
          sites_cart_used_for_pair_proxies().size()
      else:
        size = self.restraints_manager.system_size
        assert self.restraints_manager.pdb_hierarchy.atoms().size() == size
    #
    # IMPORTANT! This assumes expansion remains constant in terms of atom
    #            content (not atoms added/removed during minimization).
    #            Therefore, we do not need to re-created restraints.
    #
    if(self.restraints_manager is None or pdb_hierarchy.atoms().size()!=size):
      # This is called in expansion
      if(not self.source_of_restraints_qm()):
        model = model_from_hierarchy(
          pdb_hierarchy    = pdb_hierarchy,
          crystal_symmetry = crystal_symmetry,
          cif_objects      = self.cif_objects)
        self.restraints_manager = from_cctbx(
          restraints_manager = model.get_restraints_manager())
      else:
        assert self.source_of_restraints_qm()
        self.restraints_manager = from_qm(
          cif_objects      = self.cif_objects,
          method           = self.params.quantum.method,
          basis            = self.params.quantum.basis,
          pdb_hierarchy    = pdb_hierarchy,
          charge           = self.params.quantum.charge,
          qm_engine_name   = self.params.quantum.engine_name,
          qm_addon         = self.params.quantum.qm_addon,
          qm_addon_method  = self.params.quantum.qm_addon_method,
          memory           = self.params.quantum.memory,
          nproc            = self.params.quantum.nproc,
          url              = self.params.quantum.server_url,
          crystal_symmetry = crystal_symmetry,
          clustering       = self.params.cluster.clustering)
    return self.restraints_manager

def h_diff_sel(h1, h2):
  def getids(a):
    return a.name + a.parent().resname + a.parent().parent().resid() + \
           a.parent().parent().parent().id
  d = {}
  for a in h1.atoms():
    d[getids(a)]=a.i_seq
  #
  sel = flex.size_t()
  for a in h2.atoms():
    if not getids(a) in d:
      sel.append(a.i_seq)
  assert sel.size() == h2.atoms().size() - h1.atoms().size()
  result = ~flex.bool(h2.atoms().size(), sel)
  return result

class from_expansion(object):
  def __init__(self, params, restraints_source, pdb_hierarchy, crystal_symmetry):
    self.restraints_manager  = restraints_source.restraints_manager
    self.restraints_source   = restraints_source
    self.pdb_hierarchy       = pdb_hierarchy
    self.crystal_symmetry    = crystal_symmetry
    self.params              = params
    self.size                = self.pdb_hierarchy.atoms().size()
    #
    self.pdb_hierarchy.atoms().reset_i_seq() # XXX May be unnecessary
    #
    self.expansion = expand(
      pdb_hierarchy        = self.pdb_hierarchy.deep_copy(),
      crystal_symmetry     = self.crystal_symmetry,
      select_within_radius = self.params.cluster.select_within_radius)

    if restraints_source.source_of_restraints_qm():
      self.pdb_hierarchy_super_completed = model_completion.run(
        pdb_hierarchy          = self.expansion.ph_super_sphere.deep_copy(),
        crystal_symmetry       = self.expansion.cs_box,
        model_completion       = False,
        original_pdb_filename  = None,
        append_to_end_of_model = False, #XXX
        use_reduce             = self.params.use_reduce)
    else:
      self.pdb_hierarchy_super_completed = self.expansion.ph_super_sphere.deep_copy()

    # Selection of master copy
    selection = flex.bool(
      self.pdb_hierarchy_super_completed.atoms().size(), False)
    self.selection = selection.set_selected(
      flex.size_t(range(self.pdb_hierarchy.atoms().size())), True)
    # At this point here we are sure the model is complete. So make sure the
    # call above only changes (completes) the explansion part and leaves the
    # master copy intact!
    sites_cart_master = self.pdb_hierarchy.atoms().extract_xyz()
    sites_cart_all = self.pdb_hierarchy_super_completed.atoms().extract_xyz()
    sites_cart_all = sites_cart_all.set_selected(self.selection, sites_cart_master)
    self.pdb_hierarchy_super_completed.atoms().set_xyz(sites_cart_all)
    #
    self.pdb_hierarchy.atoms().reset_i_seq()
    self.expansion.ph_super_sphere.atoms().reset_i_seq()
    self.pdb_hierarchy_super_completed.atoms().reset_i_seq()


    self.restraints_manager = self.restraints_source.update(
        pdb_hierarchy    = self.pdb_hierarchy_super_completed,
        crystal_symmetry = self.expansion.cs_box)

    self.h_diff_sel = h_diff_sel(
      h1 = self.expansion.ph_super_sphere,
      h2 = self.pdb_hierarchy_super_completed)

  def __call__(self, selection_and_sites_cart):
    # XXX Not sure what this is and why!?
    return self.target_and_gradients(
      sites_cart = selection_and_sites_cart[1],
      selection  = selection_and_sites_cart[0],
      index      = selection_and_sites_cart[2])

  def target_and_gradients(self, sites_cart, selection=None, index=None):
    self._update(sites_cart = sites_cart)
    energy, gradients = self.restraints_manager.target_and_gradients(
        sites_cart = self.pdb_hierarchy_super_completed.atoms().extract_xyz())
    gradients = gradients.select(self.selection)
    energy=0 # XXX UNDEFINED in most cases (or may be even all cases since
             # selection above)
    #
    return energy, gradients

  def energies_sites(self, sites_cart, compute_gradients=True):
    tg = self.target_and_gradients(sites_cart=sites_cart)
    return group_args(
      target    = tg[0],
      gradients = tg[1])

  def _update(self, sites_cart):
    # Update coordinates in main hierarchy
    self.pdb_hierarchy.atoms().set_xyz(sites_cart)
    # Propagate the update into expanded hierarchy
    self.expansion.update(sites_cart = self.pdb_hierarchy.atoms().extract_xyz())
    xyz = self.pdb_hierarchy_super_completed.atoms().extract_xyz()
    xyz = xyz.set_selected(self.h_diff_sel,
      self.expansion.ph_super_sphere.atoms().extract_xyz())
    self.pdb_hierarchy_super_completed.atoms().set_xyz(xyz)

#-------------------------------------------------------------------------------

def get_cctbx_gradients(ph, cs, rm_only=False):
  model = model_from_hierarchy(
    pdb_hierarchy    = ph,
    crystal_symmetry = cs,
    cif_objects      = None)
  rm = model.get_restraints_manager()
  if(not rm_only):
    gradients = rm.energies_sites(
      sites_cart=ph.atoms().extract_xyz(), compute_gradients=True).gradients
    return group_args(model = model, gradients = gradients)
  else:
    def get_g(sites_cart):
      target = 0 # pretend it is undefined
      return target, rm.energies_sites(
        sites_cart=sites_cart, compute_gradients=True).gradients
    return get_g

class from_altlocs2(object):
  def __init__(self, model, method, params=None):
    adopt_init_args(self, locals())
    self.cs = self.model.crystal_symmetry()
    ph = self.model.get_hierarchy()
    self.d = OrderedDict()
    self.conf_ind  = ph.get_conformer_indices().conformer_indices
    self.n_altlocs = flex.max(self.conf_ind)
    for ci in set(self.conf_ind): # ci=0 is for blanc
      sel_ci = self.conf_ind == ci
      sel = sel_ci if ci == 0 else sel_ci | (self.conf_ind == 0)
      model_conformer = self.model.select(sel)
      ph_conformer = model_conformer.get_hierarchy()
      ci_ph_conformer = ph_conformer.get_conformer_indices().conformer_indices
      rm = self._setup_restraints_managers(model = model_conformer)
      self.d[ci] = group_args(
        c_selection          = sel,
        c_zero               = ci_ph_conformer == 0,
        c_one                = ci_ph_conformer == 1,
        sel_ci               = sel_ci,
        target_and_gradients = rm)
    self.sel_W_empty = \
      True if not 0 in self.d.keys() else self.d[0].c_selection.count(True) == 0

  def _setup_restraints_managers(self, model):
    restraints_source = restraints(
      params = self.params, model = model)
    return from_expansion(
        params            = self.params,
        restraints_source = restraints_source,
        pdb_hierarchy     = model.get_hierarchy(),
        crystal_symmetry  = model.crystal_symmetry()).target_and_gradients

  def target_and_gradients(self, sites_cart):
    self.model.set_sites_cart(sites_cart=sites_cart)
    g_result  = flex.vec3_double(self.conf_ind.size(), [0,0,0])
    if not self.sel_W_empty:
      g_blanks  = flex.vec3_double(self.d[0].c_selection.count(True))
    for ci, v in zip(self.d.keys(), self.d.values()):
      if ci==0: continue
      _, g_ci_blank_ = v.target_and_gradients(
        sites_cart = sites_cart.select(v.c_selection))
      g_ci = g_ci_blank_.select(v.c_one)
      g_result = g_result.set_selected(v.sel_ci, g_ci)
      if not self.sel_W_empty: g_blanks += g_ci_blank_.select(v.c_zero)
    if self.method=="subtract":
      if not self.sel_W_empty:
        W = self.d[0]
        _, g_blank = W.target_and_gradients(
          sites_cart = sites_cart.select(W.c_selection))
        result = g_result.set_selected(
          W.c_selection, g_blanks-((self.n_altlocs-1)*g_blank))
      else: result = g_result
    elif self.method=="average":
      if not self.sel_W_empty:
        result = g_result.set_selected(
          self.d[0].c_selection, g_blanks*(1/self.n_altlocs))
      else:
        result = g_result
    else: assert 0
    energy=0 # undefined!

    # XXX tmp debugging info
    N = flex.mean( flex.sqrt((result).dot()) )
    print("<|gradient|>", self.method, N)

    return energy, result

def from_cctbx_altlocs(ph, cs, method="subtract", option=2):
  assert method in ["subtract", "average"]
  g_result = flex.vec3_double(ph.atoms().size(), [0,0,0])
  conf_ind = ph.get_conformer_indices().conformer_indices
  n_altlocs = flex.max(conf_ind)
  sel_W = conf_ind == 0
  sel_W_empty = sel_W.count(True) == 0
  g_blanks = flex.vec3_double(sel_W.count(True))
  for ci in range(1, n_altlocs+1):
    sel_ci = conf_ind == ci
    ph_conformer = ph.select((conf_ind == ci) | sel_W)
    ci_ph_conformer = ph_conformer.get_conformer_indices().conformer_indices
    g_ci_blank_ = get_cctbx_gradients(ph=ph_conformer, cs=cs).gradients
    g_ci = g_ci_blank_.select(ci_ph_conformer == 1)
    g_result = g_result.set_selected(sel_ci, g_ci)
    g_blanks += g_ci_blank_.select(ci_ph_conformer == 0)
  if(method=="subtract"):
    """
    Not suitable for QM as this needs to calculate gradients using the whole model
    """
    # Option 1
    if(option==1):
      if(not sel_W_empty):
        g_blank = get_cctbx_gradients(ph=ph.select(sel_W), cs=cs).gradients
        result = g_result.set_selected(sel_W, g_blanks-((n_altlocs-1)*g_blank))
      else:
        # Both are exwctly equivalent
        #result = get_cctbx_gradients(ph=ph, cs=cs).gradients
        result = g_result
    # Option 2
    if(option==2):
      g_blank = get_cctbx_gradients(ph=ph, cs=cs).gradients.select(sel_W)
      result = g_result.set_selected(sel_W, g_blank)
    #
    # Options 1 and 2 are identical. Disabled for performance and because it
    # expectedly crashes when altloc is ' '.
    # assert approx_equal(g_result_1, g_result_2)
  elif(method=="average"):
    result = g_result.set_selected(sel_W, g_blanks*(1/n_altlocs))
    # DEBUG ph.select(sel_W).write_pdb_file("sel_W.pdb")
  else: assert 0
  return result

class from_altlocs(object):
  def __init__(self, restraints_source, pdb_hierarchy, crystal_symmetry,
               method):
    self.restraints_manager = restraints_source.restraints_manager
    self.restraints_source  = restraints_source
    self.pdb_hierarchy      = pdb_hierarchy
    self.crystal_symmetry   = crystal_symmetry
    self.method             = method

  def __call__(self, selection_and_sites_cart):
    return self.target_and_gradients(
      sites_cart = selection_and_sites_cart[1],
      selection  = selection_and_sites_cart[0],
      index      = selection_and_sites_cart[2])

  def target_and_gradients(self, sites_cart, selection=None, index=None):
    gradient = from_cctbx_altlocs(
      ph=self.pdb_hierarchy, cs=self.crystal_symmetry, method=self.method)
    energy=None # undefined!
    return energy, gradient

#-------------------------------------------------------------------------------

class from_cctbx(object):
  def __init__(self, restraints_manager, fragment_extracts=None,
              file_name="./ase/tmp_ase.pdb"):
    self.geometry_restraints_manager = restraints_manager
    self.file_name = file_name
    self.fragment_extracts = fragment_extracts
    self.total_time = 0
    self.number_of_target_and_gradients_calls = 0

  def __call__(self, selection_and_sites_cart):
    return self.target_and_gradients(
      sites_cart = selection_and_sites_cart[1],
      selection  = selection_and_sites_cart[0],
      index      = selection_and_sites_cart[2])

  def select(self, selection):
    grm = self.geometry_restraints_manager.select(selection = selection)
    return from_cctbx(restraints_manager = grm)

  def energies_sites(self, sites_cart, compute_gradients=True):
    tg = self.target_and_gradients(sites_cart=sites_cart)
    return group_args(
      target    = tg[0],
      gradients = tg[1])

  def target_and_gradients(self, sites_cart, selection=None, index=None):

    if(selection is not None): ### clustering
      super_selection = self.fragment_extracts.fragment_super_selections[index]
      grm = self.fragment_extracts.super_sphere_geometry_restraints_manager

      #es = grm.energies_sites(
      #  sites_cart=sites_cart, compute_gradients=True)
      #es.gradients = es.gradients.select(super_selection)[:selection.count(True)]
      # Is this the same?
      es = grm.select(super_selection).energies_sites(
        sites_cart=sites_cart.select(super_selection), compute_gradients=True)
      es.gradients = es.gradients[:selection.count(True)]
      es.gradients = es.gradients * flex.double(
            self.fragment_extracts.fragment_scales[index])
    else:
      es = self.geometry_restraints_manager.energies_sites(
        sites_cart=sites_cart, compute_gradients=True)
    return es.target, es.gradients

class from_qm(object):
  def __init__(self,
      fragment_extracts          = None,
      pdb_hierarchy              = None,
      charge                     = None,
      qm_engine_name             = None,
      qm_addon                   = None,
      qm_addon_method            = None,
      file_name                  = "./ase/tmp_ase.pdb",
      crystal_symmetry           = None,
      clustering                 = False,
      cif_objects                = None,
      # change to quantum phil scope !!!!
      method                     = 'rhf',
      basis                      = "sto-3g",
      memory                     = None,
      nproc                      = 1,
      url                        = None
  ):
    self.fragment_extracts  = fragment_extracts
    self.method = method
    self.basis = basis
    self.memory = memory
    self.nproc = nproc
    self.qm_addon = qm_addon
    self.qm_addon_method = qm_addon_method
    self.url = url

    self.crystal_symmetry = crystal_symmetry
    self.pdb_hierarchy = pdb_hierarchy
    self.qm_engine_name = qm_engine_name
    self.file_name = file_name
    self.working_folder = os.path.split(self.file_name)[0]+ "/"
    if(os.path.exists(self.working_folder) is not True):
      os.mkdir(self.working_folder)
    if(charge is None and clustering is False):
      charge_service = charges_class(
        pdb_hierarchy = self.pdb_hierarchy.deep_copy(), # XXX Something bad happens otherwise!
        crystal_symmetry = crystal_symmetry,
        cif_objects = cif_objects)
      self.charge = charge_service.get_total_charge()
    else: self.charge = charge
    self.clustering = clustering
    self.qm_engine = self.create_qm_engine()
    self.system_size = self.pdb_hierarchy.atoms_size()

  def create_qm_engine(self):
    if(self.qm_engine_name == "turbomole"):
      calculator = Turbomole()
    elif(self.qm_engine_name == "terachem"):
      ### if TeraChem has problem reading pdb file, update TeraChem version.
      calculator = TeraChem(gpus="4",
                            basis=self.basis,
                            dftd="no",
                            watcheindiis="yes",
                            scf="diis+a")
    elif(self.qm_engine_name == "mopac"):
      calculator = Mopac()
    elif(self.qm_engine_name == "pyscf"):
      calculator = Pyscf()
    elif(self.qm_engine_name == "orca"):
      calculator = Orca()
    elif(self.qm_engine_name == "gaussian"):
      calculator = Gaussian()
    elif(self.qm_engine_name == "torchani"):
      from .plugin.ase.torchani_qr import TorchAni
      calculator = TorchAni()
    elif(self.qm_engine_name == "aimnet2-old"):
      from .plugin.ase.aimnet2_qr_old import AIMNet2CalculatorOLD
      calculator = AIMNet2CalculatorOLD()
    elif(self.qm_engine_name == "aimnet2"):
      from .plugin.ase.aimnet2_qr import AIMNet2Calculator
      calculator = AIMNet2Calculator('aimnet2-qr')
    elif(self.qm_engine_name == "xtb"):
      calculator = GFNxTB()
    elif(self.qm_engine_name == "server"):
      calculator = RestAPICalculator(url=self.url)
    else:
      raise Sorry("qm_calculator needs to be specified.")
    #
    # set to appropriate values
    #
    for attr in ['charge',
                 'basis',
                 'method',
                 'memory',
                 'nproc',
                 ]:
      value = getattr(self, attr, None)
      func = getattr(calculator, 'set_%s' % attr, None)
      action=False
      if func is not None:
        if value is not None:
          #print '  Setting %s to %s' % (attr, value)
          func(value)
          action=True
      # XXX Avoid bare prints. Fix by propagating log channel.
      #if not action:
      #  if value and not func:
      #    print '  No function available to set %s to %s' % (attr, value)
    return calculator

  def __call__(self,fragment_selection_and_sites_cart):
    return self.target_and_gradients(
      sites_cart = fragment_selection_and_sites_cart[1],
      selection  = fragment_selection_and_sites_cart[0],
      index      = fragment_selection_and_sites_cart[2])

  def energies_sites(self, sites_cart, compute_gradients=True):
    tg = self.target_and_gradients(sites_cart=sites_cart)
    return group_args(
      target    = tg[0],
      gradients = tg[1])

  def target_and_gradients(self,sites_cart, selection=None, index=None):
    if(self.clustering):
      from .fragment import get_qm_file_name_and_pdb_hierarchy
      from .fragment import charge
      from .fragment import write_mm_charge_file
      #
      qm_pdb_file, ph = get_qm_file_name_and_pdb_hierarchy(
                          fragment_extracts=self.fragment_extracts,
                          index=index)
      #
      qm_charge = charge(fragment_extracts=self.fragment_extracts,
                                      index=index)
      charge_file =  write_mm_charge_file(fragment_extracts=self.fragment_extracts,
                                      index=index)
      gradients_scale = self.fragment_extracts.fragment_scales[index]
    else:
      self.pdb_hierarchy.atoms().set_xyz(sites_cart)
      self.pdb_hierarchy.write_pdb_file(file_name=self.file_name)
      ph = self.pdb_hierarchy## return pdb_hierarchy
      qm_pdb_file = self.file_name
      qm_charge = self.charge
      charge_file = None
      selection =flex.bool(self.system_size, True)
      gradients_scale = [1.0]*self.system_size
    define_str=''
    atoms = ase_atoms_from_pdb_hierarchy(ph, self.crystal_symmetry, self.qm_engine_name)
    unit_convert = ase_units.mol/ase_units.kcal # ~ 23.06
    self.qm_engine.set_label(qm_pdb_file[:-4])
    cwd = os.getcwd()

    #FOR DEBUGGING distance check
    # print ''
    # print '*distance check before QM calc*'
    # thr=0.6
    # for i in range(0,len(atoms)-1):
    #   for j in range(i,len(atoms)):
    #       if i==j: continue
    #       x=atoms[i].position[0]-atoms[j].position[0]
    #       y=atoms[i].position[1]-atoms[j].position[1]
    #       z=atoms[i].position[2]-atoms[j].position[2]
    #       dist=math.sqrt(x*x+y*y+z*z)
    #       if(dist<=thr):
    #         print 'WARNING: atoms ', i,j,' are closer than', thr,' A -> ',dist
    self.qm_engine.run_qr(atoms,
                          charge       = qm_charge,
                          pointcharges = charge_file,
                          coordinates  = qm_pdb_file[:-4]+".xyz",
                          define_str   = define_str, # for Turbomole
      )
    os.chdir(cwd)
    if self.qm_addon is not None:
      tool_e,tool_g= qr_tools.qm_toolbox(atoms,
                              charge=qm_charge,
                              pointcharges=charge_file,
                              label=qm_pdb_file[:-4],
                              addon=self.qm_addon,addon_method=self.qm_addon_method)
      energy = (self.qm_engine.energy_free+tool_e)*unit_convert
      ase_gradients = (tool_g-self.qm_engine.forces)*unit_convert
    else:
      energy = self.qm_engine.energy_free*unit_convert
      ase_gradients = (-1.0) * self.qm_engine.forces*unit_convert
    # remove capping and neigbouring buffer
    gradients = ase_gradients[:selection.count(True)]
    gradients =  flex.vec3_double(gradients)
    ## TODO
    ## unchange the altloc gradient, averagely scale the non-altloc gradient
    gradients = gradients*flex.double(gradients_scale)
    return energy, gradients

from ase import Atoms
def ase_atoms_from_pdb_hierarchy(ph, crystal_symmetry, qm_engine_name):

  def read_pdb_hierarchy(pdb_hierarchy):
    positions = []
    symbols = []
    for chain in pdb_hierarchy.chains():
      for residue_group in chain.residue_groups():
        for atom in residue_group.atoms():
          element = atom.element.strip()
          if (len(element) == 2):
            element = element[0] + element[1].lower()
          symbols.append(element)
          positions.append(list(atom.xyz))
    return symbols, positions
  symbols, positions = read_pdb_hierarchy(ph)
  if(qm_engine_name == "torchani"):
    unit_cell = crystal_symmetry.unit_cell().parameters()
    return Atoms(symbols=symbols, positions=positions, pbc=True, cell=unit_cell)
  else:
    #
    # XXX Ugly work-around to by-pass inability of ASE to handle D atoms.
    #
    symbols = ["H" if symbol=="D" else symbol for symbol in symbols]
    return Atoms(symbols=symbols, positions=positions)
