from __future__ import division

import os
import ase.units as ase_units
import mmtbx.restraints
from libtbx.utils import Sorry
from charges import charges_class
from scitbx.array_family import flex
from plugin.ase.mopac_qr import Mopac
from plugin.ase.pyscf_qr import Pyscf
from plugin.ase.terachem_qr import TeraChem
from plugin.ase.turbomole_qr import Turbomole
from plugin.ase.orca_qr import Orca
from plugin.ase.gaussian_qr import Gaussian
from plugin.ase.xtb_qr import GFNxTB
from plugin.tools import qr_tools
from libtbx import group_args
import math
from qrefine.super_cell import expand
import qrefine.completion as model_completion
from libtbx.utils import null_out

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
    if(not self.source_of_restraints_qm()):
      model = mmtbx.model.manager(
        model_input       = None,
        restraint_objects = self.cif_objects,
        pdb_hierarchy     = pdb_hierarchy,
        process_input     = True,
        crystal_symmetry  = crystal_symmetry,
        pdb_interpretation_params = self.pi_params,
        log               = null_out())
      model.process_input_model(make_restraints=True, grm_normalization=False)
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
        crystal_symmetry = crystal_symmetry,
        clustering       = self.params.cluster.clustering)
    return self.restraints_manager

class from_expansion(object):
  def __init__(self, restraints_source, pdb_hierarchy, crystal_symmetry):
    self.restraints_manager = restraints_source.restraints_manager
    self.restraints_source  = restraints_source
    self.pdb_hierarchy      = pdb_hierarchy
    self.crystal_symmetry   = crystal_symmetry
    self.pdb_hierarchy_super_completed = None
    self.selection = None
    self.size = self.pdb_hierarchy.atoms().size()
    self.crystal_symmetry_ss = None
    self._expand()
    self.sites_cart_previous = self.pdb_hierarchy.atoms().extract_xyz()

  def __call__(self, selection_and_sites_cart):
    return self.target_and_gradients(
      sites_cart = selection_and_sites_cart[1],
      selection  = selection_and_sites_cart[0],
      index      = selection_and_sites_cart[2])

  def target_and_gradients(self, sites_cart, selection=None, index=None):
    self._update(sites_cart = sites_cart)
    energy, gradients = self.restraints_manager.target_and_gradients(
        sites_cart = self.pdb_hierarchy_super_completed.atoms().extract_xyz())
    gradients = gradients.select(self.selection)
    return energy, gradients

  def energies_sites(self, sites_cart, compute_gradients=True):
    tg = self.target_and_gradients(sites_cart=sites_cart)
    return group_args(
      target    = tg[0],
      gradients = tg[1])

  def _update(self, sites_cart, threshold = 0.1):
    shift_max = flex.max(
      flex.sqrt((sites_cart - self.sites_cart_previous).dot()))
    if(shift_max > threshold):
      self.pdb_hierarchy.atoms().set_xyz(sites_cart)
      self._expand()
      self.sites_cart_previous = sites_cart

  def _expand(self):
    expansion = expand(
      pdb_hierarchy        = self.pdb_hierarchy,
      crystal_symmetry     = self.crystal_symmetry,
      select_within_radius = 10.0)
    pdb_hierarchy_super = expansion.ph_super_sphere
    pdb_hierarchy_super.write_pdb_file(file_name="supersphere.pdb",
      crystal_symmetry = expansion.cs_box)
    self.crystal_symmetry_ss = expansion.cs_box
    # Use same route for CCTBX and QM !
    # if(self.restraints_source.source_of_restraints_qm()):
    if 1:
      self.pdb_hierarchy_super_completed = model_completion.run(
        #pdb_hierarchy         = pdb_hierarchy_super,
        crystal_symmetry      = expansion.cs_box,
        model_completion      = False,
        pdb_filename          = "supersphere.pdb",
        original_pdb_filename = None)
    else:
      self.pdb_hierarchy_super_completed = pdb_hierarchy_super
    selection = flex.bool(
      self.pdb_hierarchy_super_completed.atoms().size(), False)
    self.selection = selection.set_selected(
      flex.size_t(xrange(self.pdb_hierarchy.atoms().size())), True)
    self.restraints_manager = self.restraints_source.update(
      pdb_hierarchy    = self.pdb_hierarchy_super_completed,
      crystal_symmetry = expansion.cs_box)

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

    # XXX This is very dirty, just to get things runing for now.

    def get_grm(ph, cs):
      import mmtbx.monomer_library.server
      import mmtbx.monomer_library.pdb_interpretation
      from mmtbx import monomer_library

      mon_lib_srv = mmtbx.monomer_library.server.server()
      ener_lib    = mmtbx.monomer_library.server.ener_lib()
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

    def get_grads2(ph, cs):
      grm = get_grm(ph=ph, cs=cs)
      es = grm.energies_sites(sites_cart=ph.atoms().extract_xyz(),
        compute_gradients=True)
      return es.gradients

    cs = self.crystal_symmetry
    ph = self.pdb_hierarchy
    g_result = flex.vec3_double(ph.atoms().size(), [0,0,0])
    conf_ind = ph.get_conformer_indices()
    n_altlocs = flex.max(conf_ind)
    sel_W = conf_ind == 0
    g_blanks = flex.vec3_double(sel_W.count(True))
    for ci in range(1, n_altlocs+1):
      sel_ci = conf_ind == ci
      sel_ci_blank = (conf_ind == ci) | sel_W
      ph_ci_blank = ph.select(sel_ci_blank)
      g_ci_blank_ = get_grads2(ph=ph_ci_blank, cs=cs)
      g_ci = g_ci_blank_.select(ph_ci_blank.get_conformer_indices() == 1)
      g_result = g_result.set_selected(sel_ci, g_ci)
      g_blanks += g_ci_blank_.select( ph_ci_blank.get_conformer_indices()==0 )
    if(self.method=="subtract"):
      g_blank = get_grads2(ph=ph.select(sel_W), cs=cs)
      g_result = g_result.set_selected(sel_W, g_blanks-((n_altlocs-1)*g_blank))
    elif(self.method=="average"):
      g_result = g_result.set_selected(sel_W, g_blanks*(1/n_altlocs))
    else: assert 0

    # Make sure gradient_only is set correctly!!!
    energy=0 # undefined!
    return energy, g_result



class from_cctbx(object):
  def __init__(self, restraints_manager, fragment_extracts=None,
              file_name="./ase/tmp_ase.pdb"):
    self.geometry_restraints_manager = restraints_manager
    self.file_name = file_name
    self.fragment_extracts = fragment_extracts

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
#      charge_service             = None,
      cif_objects                = None,
      # change to quantum phil scope !!!!
      method                     = 'rhf',
      basis                      = "sto-3g",
      memory                     = None,
      nproc                      = None,
  ):
    self.fragment_extracts  = fragment_extracts
    self.method = method
    self.basis = basis
    self.memory = memory
    self.nproc = nproc
    self.qm_addon = qm_addon
    self.qm_addon_method = qm_addon_method

    self.pdb_hierarchy = pdb_hierarchy
    self.qm_engine_name = qm_engine_name
    self.file_name = file_name
    self.working_folder = os.path.split(self.file_name)[0]+ "/"
    if(os.path.exists(self.working_folder) is not True):
      os.mkdir(self.working_folder)
    if(charge is None and clustering is False):
      #raw_records = pdb_hierarchy.as_pdb_string(crystal_symmetry=crystal_symmetry)
      #cc = charges_class(raw_records=raw_records)
      #self.charge = cc.get_total_charge()
      #@Nigel
      raw_records = pdb_hierarchy.as_pdb_string(crystal_symmetry=crystal_symmetry)
      charge_service = charges_class(raw_records = raw_records,
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
    elif(self.qm_engine_name == "ani"):
      from plugin.ase.ani_qr import Ani
      calculator = Ani()
    elif(self.qm_engine_name == "torchani"):
      from plugin.ase.torchani_qr import TorchAni
      calculator = TorchAni(method=self.method)
    elif(self.qm_engine_name == "xtb"):
      calculator = GFNxTB()
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
      from fragment import get_qm_file_name_and_pdb_hierarchy
      from fragment import charge
      from fragment import write_mm_charge_file
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
    atoms = ase_atoms_from_pdb_hierarchy(ph)
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
def ase_atoms_from_pdb_hierarchy(ph):

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
  return Atoms(symbols=symbols, positions=positions)
