"""
Total formal charge of a model (mmtbx.ligands.electrons) and point charges
of the atoms around a QM fragment (monomer library partial charges) for
charge embedding.
"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
from libtbx.utils import Sorry
from mmtbx.monomer_library import server
from iotbx.pdb import amino_acid_codes as aac

get_class = iotbx.pdb.common_residue_names_get_class

default_ion_charges = {
  "PT" : 2,
  "CL" : -1,
  "CD" : 2,
  "ZN" : 2,
  'MG' : 2,
  'CU' : 2,
  }

def get_mon_lib_server(cif_objects=None):
  """
  Monomer library server that knows the restraints of cif_objects (a list of
  (file_name, cif_object) pairs as mmtbx.model.manager.get_restraint_objects
  returns them).
  """
  mon_lib_server = server.server()
  if cif_objects:
    for cif_object in cif_objects:
      mon_lib_server.process_cif_object(cif_object[1])
  return mon_lib_server

def get_cif_objects(ligand_cif_file_names=None):
  cif_objects = []
  if ligand_cif_file_names:
    for fn in ligand_cif_file_names:
      cif_object = iotbx.cif.reader(fn, strict=False).model()
      cif_objects.append((fn, cif_object))
  return cif_objects

class charges_class:
  """
  pdb_hierarchy and crystal_symmetry, or pdb_filename (the hierarchy and
  symmetry are read from the file; crystal_symmetry overrides the file's).
  """
  def __init__(self,
               pdb_hierarchy=None,
               crystal_symmetry=None,
               pdb_filename=None,
               ligand_cif_file_names=None,
               cif_objects=None,
               verbose=False,
               ):
    self.verbose = verbose
    assert not (ligand_cif_file_names and cif_objects)
    if ligand_cif_file_names:
      cif_objects = get_cif_objects(ligand_cif_file_names)
    self.cif_objects = cif_objects
    assert [pdb_hierarchy, pdb_filename].count(None) == 1
    if pdb_filename is not None:
      pdb_inp = iotbx.pdb.input(file_name=pdb_filename)
      pdb_hierarchy = pdb_inp.construct_hierarchy()
      if crystal_symmetry is None:
        crystal_symmetry = pdb_inp.crystal_symmetry()
    self.update_pdb_hierarchy(
      pdb_hierarchy=pdb_hierarchy,
      crystal_symmetry=crystal_symmetry)
    assert self.crystal_symmetry, 'There is no CRYST1 record in the input file'
    self.hetero_charges = default_ion_charges
    self.mon_lib_server = get_mon_lib_server(cif_objects)

  def update_pdb_hierarchy(self, pdb_hierarchy, crystal_symmetry):
    self.pdb_hierarchy = pdb_hierarchy
    self.crystal_symmetry = crystal_symmetry
    # merge atoms from clustering
    self.pdb_hierarchy.merge_atoms_at_end_to_residues()

  def get_total_charge(self, list_charges=False, verbose=False):
    from mmtbx.ligands import electrons
    return electrons.run(pdb_hierarchy=self.pdb_hierarchy,
                         crystal_symmetry=self.crystal_symmetry,
                         cif_objects=self.cif_objects,
                         return_formal_charges=list_charges,
                         verbose=verbose or self.verbose)

  def write_pdb_hierarchy_qxyz_file(self,
                                    file_name="qxyz_cctbx.dat",
                                    exclude_water=True,
                                    charge_scaling_positions=None,
                                    scale=0,
                                    ):
    self.write_charge_and_coordinates_from_hierarchy(
      file_name=file_name,
      qxyz_order='qxyz',
      exclude_water=exclude_water,
      charge_scaling_positions=charge_scaling_positions,
      scale=scale)

  def write_pdb_hierarchy_xyzq_file(self,
                                    file_name="xyzq_cctbx.dat",
                                    exclude_water=True,
                                    charge_scaling_positions=None,
                                    scale=0,
                                    ):
    self.write_charge_and_coordinates_from_hierarchy(
      file_name=file_name,
      qxyz_order='xyzq',
      exclude_water=exclude_water,
      charge_scaling_positions=charge_scaling_positions,
      scale=scale)

  def write_charge_and_coordinates_from_hierarchy(self,
                                                  file_name,
                                                  qxyz_order='qxyz',
                                                  exclude_water=True,
                                                  charge_scaling_positions=None,
                                                  scale=0,
                                                  assert_no_alt_loc=False,
                                                  ):
    if qxyz_order not in ('qxyz', 'xyzq'):
      raise Sorry('invalid qxyz_order parameter "%s"' % qxyz_order)
    qxyz = []
    for rg in self.pdb_hierarchy.residue_groups():
      if exclude_water:
        if get_class(rg.atom_groups()[0].resname) == 'common_water': continue
      if assert_no_alt_loc and len(rg.atom_groups()) > 1:
        raise Sorry("Contains alt. locs.")
      qxyz += get_partial_point_charges(rg, self.mon_lib_server,
                                        hetero_charges=self.hetero_charges)
    if len(qxyz) == 0: return
    scale_partial_point_charges(qxyz, charge_scaling_positions, scale=scale)
    outl = ""
    with open(file_name, "w") as qxyz_file:
      if qxyz_order == 'qxyz': # terachem
        qxyz_file.write(str(self.pdb_hierarchy.atoms_size()) + "  \n")
        qxyz_file.write("  \n")
      for i, item in enumerate(qxyz):
        if item[0] is None: ## point charge may be set as zero
          outl += ' %s has zero/None partial charge\n' % \
            self.pdb_hierarchy.atoms()[i].quote()
        if qxyz_order == 'qxyz':
          item_list = item + ["  \n"]
        else:
          item_list = item[1:] + item[0:1] + ["  \n"]
        qxyz_file.write("  ".join(str(elm) for elm in item_list))
    if outl:
      print('WARNINGS')
      print(outl)
      raise Sorry('point charges are not set.')

def _get_restraints_from_resname(resname, mon_lib_server):
  input_resname = resname
  restraints = mon_lib_server.get_comp_comp_id_direct(resname)
  if restraints is None:
    resname = aac.three_letter_l_given_three_letter_d.get(resname, None)
    if resname is not None:
      restraints = mon_lib_server.get_comp_comp_id_direct(resname)
  assert restraints, 'no restraints for "%s" found' % input_resname
  return restraints

def get_partial_point_charges(rg,
                              mon_lib_server,
                              hetero_charges=None):
  """
  [charge, x, y, z] of every atom of the residue group rg from the monomer
  library partial charges (ions: the charge of the atom record, else
  hetero_charges by residue name).
  """
  # need to handle HT1, HT2 and HXT
  v2_to_3 = {' HA3':' HA1',
             ' HB3':' HB1',
             ' HG3':' HG1',
             'HG13':'HG11',
             ' HD3':' HD1',
             ' HE3':' HE1',
             ' HT1':' HA1', #TODO check!
             ' HT2':' HA2', #TODO check!
             }
  misc = {' OXT' : ' O  ',
          }
  tmp = []
  for ag in rg.atom_groups():
    restraints = _get_restraints_from_resname(ag.resname, mon_lib_server)
    atom_dict = restraints.atom_dict()
    for atom in ag.atoms():
      # ions
      if get_class(ag.resname)=="common_element":
        assert len(ag.atoms())==1
        if not atom.charge.strip():
          if hetero_charges:
            key = atom.parent().resname
            charge = hetero_charges.get(key.strip(), None)
            if charge:
              tmp.append([charge]+list(atom.xyz))
            else:
              raise Sorry('no charge found in the model file or hetero_charges for "%s"' % atom.quote())
          else:
            raise Sorry('no charge found in the model file for "%s"' % atom.quote())
        else:
          tmp.append([atom.charge_as_int()]+list(atom.xyz))
        continue
      # other atoms
      cif = atom_dict.get(atom.name.strip(), None)
      if cif is None:
        if atom.name in [" H1 ", " H2 ", " H3 "]: # needs calculating...
          tmp.append([0.26]+list(atom.xyz))
          continue
        elif atom.name in [' HXT']: # zero charges for fake atoms
          tmp.append([0.00]+list(atom.xyz))
          continue
        if atom.name in v2_to_3:
          cif = atom_dict.get(v2_to_3[atom.name].strip())
        elif atom.name in misc:
          cif = atom_dict.get(misc[atom.name].strip())
        elif atom.name.find("'")>-1:
          name = atom.name.replace("'", "*")
          cif = atom_dict.get(name.strip(), None)
      assert cif, "%s" % atom_dict
      tmp.append([cif.partial_charge]+list(atom.xyz))
  return tmp

def scale_partial_point_charges(qxyz,
                                charge_scaling_positions=None,
                                scale=0):
  """
  Scale (in place) the charges of qxyz whose position is one of
  charge_scaling_positions (within 1.e-3 A) by scale.
  """
  def partial_charge_in_charge_scaling_positions(partial_charge,
                                                 charge_scaling_positions):
    for xyz in charge_scaling_positions:
      same_point = ( abs(xyz[0] - partial_charge[1]) < 1.0E-3 and
                     abs(xyz[1] - partial_charge[2]) < 1.0E-3 and
                     abs(xyz[2] - partial_charge[3]) < 1.0E-3 )
      if same_point:
        return True
    return False
  #####
  if charge_scaling_positions is not None:
    for item in qxyz:
      if partial_charge_in_charge_scaling_positions(item,
                                                    charge_scaling_positions):
        item[0] =  item[0]*scale
