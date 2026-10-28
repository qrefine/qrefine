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
from qrefine import fragment_capping
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
import time
from libtbx import easy_pickle
from qrefine.utils.hierarchy_utils import model_from_hierarchy

class restraints(object):
  """
  Create CCTBX or QM restraints manager.
  """
  def __init__(self, params, model, working_folder="ase"):
    self.params = params
    self.working_folder = working_folder
    self.model = model.deep_copy()
    self.cif_objects      = model.get_restraint_objects()
    self.pdb_hierarchy    = model.get_hierarchy()
    self.crystal_symmetry = model.crystal_symmetry()
    self.restraints_manager = None
    self.update(
      pdb_hierarchy    = model.get_hierarchy(),
      crystal_symmetry = model.crystal_symmetry())

  def source_of_restraints_qm(self):
    return self.params.restraints == "qm"

  def update(self, pdb_hierarchy, crystal_symmetry, charge=None):
    # charge: the total charge of pdb_hierarchy when the caller knows it (a
    # capped expansion); quantum.charge set by the user wins.
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
          restraints_manager = model.get_restraints_manager(),
          file_name          = os.path.join(self.working_folder, "tmp_ase.pdb"))
      else:
        assert self.source_of_restraints_qm()
        self.restraints_manager = from_qm(
          cif_objects      = self.cif_objects,
          method           = self.params.quantum.method,
          basis            = self.params.quantum.basis,
          pdb_hierarchy    = pdb_hierarchy,
          charge           = self.params.quantum.charge if
            self.params.quantum.charge is not None else charge,
          qm_engine_name   = self.params.quantum.engine_name,
          qm_addon         = self.params.quantum.qm_addon,
          qm_addon_method  = self.params.quantum.qm_addon_method,
          memory           = self.params.quantum.memory,
          nproc            = self.params.quantum.nproc,
          url              = self.params.quantum.server_url,
          crystal_symmetry = crystal_symmetry,
          clustering       = self.params.cluster.clustering,
          file_name        = os.path.join(self.working_folder, "tmp_ase.pdb"))
    return self.restraints_manager

class from_expansion(object):
  def __init__(self, params, restraints_source):
    self.model               = restraints_source.model
    self.restraints_manager  = restraints_source.restraints_manager
    self.restraints_source   = restraints_source
    self.params              = params
    self.size                = self.model.get_hierarchy().atoms().size()
    #
    self.model.get_hierarchy().atoms().reset_i_seq() # XXX May be unnecessary
    #
    self.expansion = expand(
      pdb_hierarchy        = self.model.get_hierarchy().deep_copy(),
      crystal_symmetry     = self.model.crystal_symmetry(),
      select_within_radius = self.params.cluster.select_within_radius)

    # DEBUGGING
    #self.expansion.ph_super_sphere.write_pdb_file(
    #  file_name = "ss.pdb",
    #  crystal_symmetry = self.expansion.cs_box)

    if restraints_source.source_of_restraints_qm():
      # Cap the super sphere (cut out of the super cell) as reduce2
      # capping=True would, one H per cut bond; the master is complete
      # (reduce2 qmready=True) and gets no caps. The capper re-places the
      # caps on every geometry update (_update).
      self.capper = fragment_capping.capper(
        parent_hierarchy = self.expansion.SCM.super_cell_hierarchy,
        selection        = self.expansion.SCM.keep_selection,
        crystal_symmetry = self.expansion.cs_box,
        cif_objects      = restraints_source.cif_objects)
      assert self.capper.fragment_hierarchy.is_similar_hierarchy(
        self.expansion.ph_super_sphere)
      r = self.capper.cap(
        sites_cart = self.expansion.ph_super_sphere.atoms().extract_xyz())
      self.pdb_hierarchy_super_completed = r.pdb_hierarchy
      # Selection of master copy: the first self.size fragment atoms
      self.selection = flex.bool(r.pdb_hierarchy.atoms().size(),
        r.fragment_selection.iselection()[:self.size])
      # ASSERT MASTER IS UNTOUCHED
      assert self.pdb_hierarchy_super_completed.select(
        self.selection).is_similar_hierarchy(self.model.get_hierarchy())
    else:
      self.capper = None
      self.pdb_hierarchy_super_completed = \
        self.expansion.ph_super_sphere.deep_copy()
      # Selection of master copy
      self.selection = flex.bool(
        self.pdb_hierarchy_super_completed.atoms().size(),
        flex.size_t(range(self.size)))
    #
    self.model.get_hierarchy().atoms().reset_i_seq()
    self.expansion.ph_super_sphere.atoms().reset_i_seq()
    self.pdb_hierarchy_super_completed.atoms().reset_i_seq()

    self.restraints_manager = self.restraints_source.update(
        pdb_hierarchy    = self.pdb_hierarchy_super_completed,
        crystal_symmetry = self.expansion.cs_box,
        charge           = self.capper.charge if self.capper else None)

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
    self.model.get_hierarchy().atoms().set_xyz(sites_cart)
    # Propagate the update into expanded hierarchy
    self.expansion.update(
      sites_cart = self.model.get_hierarchy().atoms().extract_xyz())
    xyz = self.expansion.ph_super_sphere.atoms().extract_xyz()
    if self.capper is not None:
      # re-place the caps on the new geometry (same atoms every time)
      xyz = self.capper.cap(sites_cart=xyz).pdb_hierarchy.atoms().extract_xyz()
    assert xyz.size() == self.pdb_hierarchy_super_completed.atoms().size()
    self.pdb_hierarchy_super_completed.atoms().set_xyz(xyz)

#-------------------------------------------------------------------------------

def get_cctbx_gradients(ph, cs):
  model = model_from_hierarchy(
    pdb_hierarchy    = ph,
    crystal_symmetry = cs,
    cif_objects      = None)
  rm = model.get_restraints_manager()
  gradients = rm.energies_sites(
    sites_cart=ph.atoms().extract_xyz(), compute_gradients=True).gradients
  return group_args(model = model, gradients = gradients)

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
    os.makedirs(self.working_folder, exist_ok=True)
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
      model = 'aimnet2-qr' if self.method == 'rhf' else self.method  # replace class default method
      calculator = AIMNet2Calculator(model)
    elif(self.qm_engine_name == "xtb"):
      calculator = GFNxTB()
    elif(self.qm_engine_name == "server"):
      calculator = RestAPICalculator(url=self.url)
    elif(self.qm_engine_name == "fairchem"):
      from .plugin.ase.fairchem import FAIRChemCalculator2
      calculator = FAIRChemCalculator2()
    else:
      raise Sorry("qm_calculator needs to be specified.")
    #
    # set to appropriate values
    #
    # aimnet2-qr plugin doesnt use these, however set_charge exists in the AIMNET2
    #  base calculator but expects an Atom object (it tries to update tensors)
    # set_charge is done inside our AIMNet2Calculator
    # This is needed for the offical https://github.com/isayevlab/AIMNet2 not Romans version
    if not self.qm_engine_name == "aimnet2":
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
      qm_pdb_file, ph, fragment_selection = get_qm_file_name_and_pdb_hierarchy(
                          fragment_extracts=self.fragment_extracts,
                          index=index)
      #
      qm_charge = charge(fragment_extracts=self.fragment_extracts,
                                      index=index)
      charge_file =  write_mm_charge_file(fragment_extracts=self.fragment_extracts,
                                      index=index)
    else:
      self.pdb_hierarchy.atoms().set_xyz(sites_cart)
      self.pdb_hierarchy.write_pdb_file(file_name=self.file_name)
      ph = self.pdb_hierarchy## return pdb_hierarchy
      qm_pdb_file = self.file_name
      qm_charge = self.charge
      charge_file = None
      selection = flex.bool(self.system_size, True)
      fragment_selection = selection
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
    # drop the capping H (fragment_selection), then the buffer atoms from
    # symmetry copies: the fragment's master atoms come first (the super
    # sphere starts with the master copy) and selection counts them
    gradients = flex.vec3_double(ase_gradients).select(fragment_selection)
    gradients = gradients[:selection.count(True)]
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
