from __future__ import division

import os
import ase.units as ase_units
import mmtbx.restraints
from libtbx.utils import Sorry
from charges import get_total_charge_from_pdb
from scitbx.array_family import flex
from clustering import betweenness_centrality_clustering
from plugin.ase.mopac_qr import Mopac
from plugin.ase.pyscf_qr import Pyscf
from plugin.ase.terachem_qr import TeraChem
from plugin.ase.turbomole_qr import Turbomole

class from_cctbx(object):
  def __init__(self, processed_pdb_file, has_hd, fragment_extracts=None,
              file_name="./ase/tmp_ase.pdb"):
    geometry = processed_pdb_file.geometry_restraints_manager(
      show_energies                = False,
      assume_hydrogens_all_missing = not has_hd,
      plain_pairs_radius           = 5.0)
    self.geometry_restraints_manager = mmtbx.restraints.manager(
       geometry = geometry, normalization = False)
    self.file_name = file_name
    self.fragment_extracts = fragment_extracts

  def __call__(self,selection_and_sites_cart):
    return self.target_and_gradients(
      sites_cart = selection_and_sites_cart[1],
      selection  = selection_and_sites_cart[0],
      index      = selection_and_sites_cart[2])

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
      file_name                  = "./ase/tmp_ase.pdb",
      #clustering_method          = betweenness_centrality_clustering,
      #maxnum_residues_in_cluster = 20,
      #charge_embedding           = False,
      crystal_symmetry           = None,
      clustering                 = False,
      #shared_disk                = True,
      basis                      = "sto-3g"):
    self.fragment_extracts  = fragment_extracts
    self.basis = basis
    self.pdb_hierarchy = pdb_hierarchy
    self.qm_engine_name = qm_engine_name
    self.file_name = file_name
    self.working_folder = os.path.split(self.file_name)[0]+ "/"
    if(os.path.exists(self.working_folder) is not True):
      os.mkdir(self.working_folder)
    if(charge is None and clustering is False):
      raw_records = pdb_hierarchy.as_pdb_string(crystal_symmetry=crystal_symmetry)
      self.charge = get_total_charge_from_pdb(raw_records=raw_records)
    else: self.charge = charge
    self.clustering = clustering
    #self.shared_disk = shared_disk
    self.qm_engine = self.create_qm_engine()
    self.system_size = self.pdb_hierarchy.atoms_size()

  def create_qm_engine(self):
    if(self.qm_engine_name == "turbomole"):
      calculator = Turbomole()
    elif(self.qm_engine_name == "terachem"):
      ### if TeraChem has problem reading pdb file, update TeraChem version.
      calculator = TeraChem( gpus="4",basis=self.basis,dftd="yes",
                      watcheindiis="yes", scf="diis+a")#
    elif(self.qm_engine_name == "mopac"):
      calculator = Mopac()
    elif(self.qm_engine_name == "pyscf"):
      calculator = Pyscf()
    else:
      raise Sorry("qm_calculator needs to be specified.")
    return calculator

  def __call__(self,fragment_selection_and_sites_cart):
    print "here __call__"
    return self.target_and_gradients(
      sites_cart = fragment_selection_and_sites_cart[1],
      selection  = fragment_selection_and_sites_cart[0],
      index      = fragment_selection_and_sites_cart[2])

  def target_and_gradients(self,sites_cart, selection=None, index=None):
    if(self.clustering):
      from fragment import get_qm_file_name_and_pdb_hierarchy
      from fragment import charge
      from fragment import write_mm_charge_file
      #
      print index
      qm_pdb_file, ph = get_qm_file_name_and_pdb_hierarchy(
                          fragment_extracts=self.fragment_extracts,
                          index=index,
                          original_pdb_filename=os.path.abspath(self.file_name))
      #
      qm_charge = charge(fragment_extracts=self.fragment_extracts,
                                      index=index)
      charge_file =  write_mm_charge_file(fragment_extracts=self.fragment_extracts,
                                      index=index)
      print charge_file
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
    define_str = '\n\na coord\n*\nno\nb all def-SV(P)\n*\neht\n\n'\
    + str(qm_charge)\
    + '\n\nscf\niter\n300\n\ncc\nmemory\n4000\n*\ndft\non\nfunc\nb-p\n*\nri\non\nm\n1000\n*\n* '
    atoms = ase_atoms_from_pdb_hierarchy(ph)
    self.qm_engine.label = qm_pdb_file[:-4]
    self.qm_engine.run_qr(atoms,charge=qm_charge, pointcharges=charge_file,
          coordinates=qm_pdb_file[:-4]+".xyz", define_str=define_str)
    unit_convert = ase_units.mol/ase_units.kcal
    energy = self.qm_engine.energy_free*unit_convert
    ase_gradients = (-1.0) * self.qm_engine.forces*unit_convert
    gradients = ase_gradients[:selection.count(True)]#remove capping and neibouring buffer
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
