import os
import sys
from libtbx.utils import Sorry
import iotbx
from mmtbx.chemical_components import get_cif_dictionary
from mmtbx.monomer_library import server

from utils import hierarchy_utils
from iotbx.pdb import amino_acid_codes as aac

get_class = iotbx.pdb.common_residue_names_get_class

default_ion_charges = {
  "PT" : 2,
  "CL" : -1,
  "CD" : 2,
  "ZN" : 2,
  'MG' : 2,
  }
allowable_amino_acid_charges = {
  "ARG" : 1,
  "CYS" : 0,
  'GLU' : -1,
  'ASP' : -1,
  }
charge_per_aa_polymer = {}
hydrogens_per_aa_polymer = {}
non_hydrogens_per_aa_polymer = {}
non_standard_amino_acids = { #"SAR" : None,
                            }

def get_mon_lib_server(ligand_cif_file_names=None):
  mon_lib_server = server.server()
  if ligand_cif_file_names:
    for fn in ligand_cif_file_names:
      #iotbx.cif.reader(input_string=ligand_cif,
      #                 cif_object=cif_object,
      #                 strict=False,
      #)
      mon_lib_server.process_cif(file_name=fn)
  return mon_lib_server

class chemical_component_class(dict):
  def get_total_charge(self):
    total = 0
    for atom in self.get("_chem_comp_atom", []):
      total += getattr(atom, 'charge', 0)
    return total

  def get_hydrogens(self):
    hs = []
    for atom in self.get("_chem_comp_atom", []):
      if getattr(atom, 'type_symbol') in ["H", "D"]:
        hs.append(atom)
    return hs

  def get_non_hydrogens(self):
    hs = []
    for atom in self.get("_chem_comp_atom", []):
      if not getattr(atom, 'type_symbol') in ["H", "D"]:
        hs.append(atom)
    return hs

class charges_class:
  def __init__(self,
               pdb_filename=None,
               raw_records=None,
               pdb_inp=None,
               ligand_cif_file_names=None,
               #list_charges=False,
               verbose=False,
               ):
    ppf = hierarchy_utils.get_processed_pdb(pdb_filename=pdb_filename,
                                            raw_records=raw_records,
                                            pdb_inp=pdb_inp,
                                           )
    self.pdb_inp = ppf.all_chain_proxies.pdb_inp
    self.pdb_hierarchy = ppf.all_chain_proxies.pdb_hierarchy
    self.crystal_symmetry = self.pdb_inp.crystal_symmetry_from_cryst1()
    assert self.crystal_symmetry, 'There is no CRYST1 record in the input file'

    self.hetero_charges = get_hetero_charges(self.pdb_inp,
                                             self.pdb_hierarchy,
    )
    if not self.hetero_charges:
      # some defaults
      self.hetero_charges = default_ion_charges
    self.inter_residue_bonds = get_inter_residue_bonds(ppf)
    if verbose:
      for key, item in inter_residue_bonds.items():
        if type(key)!=type(0) and len(key)==2: print key, item
    # merge atoms from clustering
    self.pdb_hierarchy = hierarchy_utils.merge_atoms_at_end_to_residues(
      self.pdb_hierarchy,
      )
    # needs hetero_charges?
    self.mon_lib_server = get_mon_lib_server(
      ligand_cif_file_names=ligand_cif_file_names)

  def get_total_charge(self,
                       list_charges=False,
                       check=False,
                       verbose=False,
                      ):
    total_charge = calculate_pdb_hierarchy_charge(
      self.pdb_hierarchy,
      hetero_charges=self.hetero_charges,
      inter_residue_bonds=self.inter_residue_bonds,
      list_charges=list_charges,
      check=check,
      verbose=verbose,
    )
    return total_charge

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
      charge_scaling_positions=None,
      scale=0,
    )

  def write_charge_and_coordinates_from_hierarchy(self,
                                                  file_name,
                                                  qxyz_order='qxyz',
                                                  exclude_water=True,
                                                  charge_scaling_positions=None,
                                                  scale=0,
                                                  assert_no_alt_loc=False,
                                                  ):
    qxyz = None
    for residue in hierarchy_utils.generate_residue_groups(
      self.pdb_hierarchy,
      assert_no_alt_loc=assert_no_alt_loc,
      exclude_water=exclude_water,
      ):
      if qxyz is None:
        qxyz = get_partial_point_charges(residue,
                                         self.mon_lib_server,
                                         hetero_charges=self.hetero_charges)
      else:
        qxyz = qxyz + get_partial_point_charges(residue,
                                                self.mon_lib_server,
                                                hetero_charges=self.hetero_charges)
    if qxyz is None: return
    scale_partial_point_charges(qxyz,charge_scaling_positions, scale=0)
    qxyz_file = open(file_name,"w+")
    if qxyz_order=='qxyz': # tetrachem?
      qxyz_file.write(str(self.pdb_hierarchy.atoms_size())+ "  \n")
      qxyz_file.write("  \n")
    elif qxyz_order=='xyzq':
      pass
    else:
      raise Sorry('invalid qxyz_order parameter "%s"' % qxyz_order)
    outl = ""
    for i, item in enumerate(qxyz):
      if item[0]==0 or item[0] is None:
        outl += ' %s has zero/None partial charge\n' % hierarchy.atoms()[i].quote()
      if qxyz_order=='qxyz':
        item_list = item + ["  \n"]
      elif qxyz_order=='xyzq':
        item_list = item[1:]+item[0:1] + ["  \n"]
      else:
        raise Sorry('invalid qxyz_order parameter "%s"' % qxyz_order)
      item_string = "  ".join(str(elm) for elm in item_list)
      qxyz_file.write(item_string)
    qxyz_file.close()

    if outl:
      print 'WARNINGS'
      print outl
      raise Sorry('point charges are not set.')

def get_aa_charge(code):
  # get from cache first
  # then look in the chemical components
  # not sure what to do about novel ligands...
  tmp = charge_per_aa_polymer.get(code, None)
  if tmp: return tmp
  l_peptide = acc.three_letter_l_given_three_letter_d.get(code, None)
  cc = chemical_component_class()
  if l_peptide:
    cc.update(get_cif_dictionary(l_peptide))
  else:
    cc.update(get_cif_dictionary(code))
  tmp = cc.get_total_charge()
  charge_per_aa_polymer[code] = tmp
  return tmp

def get_aa_polymer_hydrogens(code):
  tmp = hydrogens_per_aa_polymer.get(code, None)
  cc = chemical_component_class()
  cc.update(get_cif_dictionary(code))
  tmp = cc.get_hydrogens()
  hydrogens_per_aa_polymer[code] = tmp
  return tmp

def get_aa_polymer_non_hydrogens(code):
  tmp = non_hydrogens_per_aa_polymer.get(code, None)
  cc = chemical_component_class()
  cc.update(get_cif_dictionary(code))
  tmp = cc.get_non_hydrogens()
  non_hydrogens_per_aa_polymer[code] = tmp
  return tmp

def validate_non_hydrogens_atoms(rg):
  assert len(rg.atom_groups())==1
  code = rg.atom_groups()[0].resname
  tmp = get_aa_polymer_non_hydrogens(code)
  for r in tmp:
    if r.atom_id.strip() in ['OXT']: continue
    for atom in rg.atoms():
      if atom.name.strip()==r.atom_id.strip():
        break
    else:
      assert 0, r.atom_id
  return True

def validate_hydrogens_atoms(rg):
  # rewrite using sets!!!
  assert len(rg.atom_groups())==1
  code = rg.atom_groups()[0].resname
  tmp = get_aa_polymer_hydrogens(code)
  for r in tmp:
    if r.atom_id.strip() in ['H',  # not in terminal
                             'H2', # not in polymer
                             'HXT',# obsolete?
                             'HD2',# side-chain acid ASP
                             'HD1', 'HE2', # HIS
                             'HG', # CYS
                            ]: continue
    for atom in rg.atoms():
      if atom.name.strip()==r.atom_id.strip():
        break
    else:
      assert 0, r.atom_id
  return True

def validate_all_atoms(rg):
  return True
  validate_non_hydrogens_atoms(rg)
  validate_hydrogens_atoms(rg)
  return True

def calculate_residue_charge(rg,
                             assert_contains_hydrogens=True,
                             assert_no_alt_loc=True,
                             hetero_charges=None,
                             inter_residue_bonds=None,
                             verbose=False,
                             ):
  if verbose: print '-'*80
  def _terminal(names, check):
    for name in check:
      if name not in names:
        break
    else:
      return True
    return False
  def n_terminal(residue_name, atom_names):
    if residue_name in ["PRO"]:
      check_names = [[' H2 ',' H3 '],
                     [' H 1', ' H 2'], # CHARMM...
                     [' HN1', ' HN2'], # BABEL...
                     ]
    else:
      check_names = [[' H1 ',' H2 ',' H3 '],
                     ]
    for check_name in check_names:
      rc = _terminal(atom_names, check_name)
      if rc: break
    return rc
  def n_capping(residue_name, atom_names):
    if residue_name in ["PRO"]:
      check_names = [[' H2 ']]
    else:
      check_names = [[' H1 ',' H2 '],
                     [' H  ',' H2 '], # from finalise
                     ]
    for check_name in check_names:
      rc = _terminal(atom_names, check_name)
      if rc: break
    return rc
  def nh2_terminal(atom_names):
    return _terminal(atom_names, [' HT1', ' HT2'])
  def nh3_terminal(atom_names):
    return _terminal(atom_names, [' HT1', ' HT2', ' HT3'])
  def c_terminal(atom_names):
    rc = _terminal(atom_names, [' OXT'])
    if not rc: rc = _terminal(atom_names, [' OT1', ' OT2'])
    return rc
  def c_capping(atom_names):
    rc = _terminal(atom_names, [' HC '])
    return rc
  def covalent_bond(i_seqs, inter_residue_bonds):
    for i_seq in i_seqs:
      if i_seq in inter_residue_bonds:
        return True
    return False
  ############
  max_charge=1
  if assert_no_alt_loc:
    if len(rg.atom_groups())>1:
      raise Sorry("alt locs in %s" % hierarchy_utils.display_residue_group(rg))
  # ions
  # needs to be centralised!!!
  resname = rg.atom_groups()[0].resname
  if get_class(resname)=="common_element":
    atom = rg.atoms()[0]
    if not atom.charge.strip():
      if hetero_charges:
        charge = hetero_charges.get( atom.parent().resname.strip(), None)
        if charge is None:
          raise Sorry('no charge found in the model file or hetero_charges for "%s"' % atom.quote())
        else:
          return charge, charge
      else:
        raise Sorry('no charge found in the model file for "%s"' % atom.quote())
    else:
      return atom.charge_as_int(), atom.charge_as_int()
  # others
  hs=0
  atom_names = []
  atom_i_seqs = []
  for atom in rg.atoms():
    if verbose: print '...',atom.quote()
    if atom.element_is_hydrogen(): hs+=1
    atom_names.append(atom.name)
    atom_i_seqs.append(atom.i_seq)
  if verbose: print get_class(resname)
  if assert_contains_hydrogens:
    if hs==0:
      hydrogens = get_aa_polymer_hydrogens(resname)
      if len(hydrogens)!=0:
        if verbose:
          for atom in rg.atoms(): print 'H',atom.quote()
        raise Sorry("no hydrogens: %s" % hierarchy_utils.display_residue_group(rg))
  ag = rg.atom_groups()[0]
  charge = get_aa_charge(ag.resname)
  rc = get_aa_charge(ag.resname)
  if ag.resname in ['GLU', 'ASP']:
    rc=-1 # reporting only
  annot = ''
  if verbose:
    print '%s\nstarting charge: %s' % ('*'*80, charge)
  if ( get_class(ag.resname) in ["common_amino_acid", "modified_amino_acid"] or
       ag.resname in acc.three_letter_l_given_three_letter_d
       ):
    if verbose:
      print ag.id_str()
      print 'number of hydrogens',len(get_aa_polymer_hydrogens(ag.resname))
    poly_hs = len(get_aa_polymer_hydrogens(ag.resname))-2
    diff_hs = hs-poly_hs
    if verbose: print 'charge: %s poly_hs: %s diff_hs: %s' % (charge,
                                                              poly_hs,
                                                              diff_hs,
                                                            )
    if verbose: print atom_names
    if n_terminal(ag.resname, atom_names):
      diff_hs-=1
      max_charge+=1
      if verbose:
        print 'n_terminal'
        print 'charge: %s poly_hs: %s diff_hs: %s' % (charge,
                                                      poly_hs,
                                                      diff_hs,
        )
      annot += 'N-term. '
    elif nh3_terminal(atom_names):
      diff_hs-=1
      max_charge+=1
      if verbose: print 'nh3_terminal True'
      annot += 'NH3-term. '
    elif nh2_terminal(atom_names):
      diff_hs-=1
      max_charge+=1
      if verbose: print 'nh2_terminal True'
      annot += 'NH2-term. '
    elif n_capping(ag.resname, atom_names):
      diff_hs-=1
      if verbose: print 'n_capping True'
      annot += 'N-capp. '
    else:
      if verbose: print 'no N term'
    if c_terminal(atom_names):
      diff_hs-=1
      max_charge+=1
      if verbose:
        print 'c_terminal'
        print 'charge: %s poly_hs: %s diff_hs: %s' % (charge,
                                                      poly_hs,
                                                      diff_hs,
                                                    )
      annot += 'C-term. '
    elif c_capping(atom_names):
      diff_hs-=1
      #max_charge+=1
      if verbose:
        print 'c_capping'
        print 'charge: %s poly_hs: %s diff_hs: %s' % (charge,
                                                      poly_hs,
                                                      diff_hs,
                                                    )
      annot += 'C-capp. '
    else:
      if verbose: print 'no C term'
    if covalent_bond(atom_i_seqs, inter_residue_bonds):
      diff_hs+=1
      if verbose:
        print 'covalent_bond',atom_i_seqs, inter_residue_bonds
      annot += 'Coval. '
    if hierarchy_utils.is_n_terminal_atom_group(ag):
      diff_hs-=1
    if verbose:
      print 'residue: %s charge: %s poly_hs: %2s diff_hs: %2s total: %2s %s' % (
        ag.resname,
        charge,
        poly_hs,
        diff_hs,
        charge+diff_hs,
        annot,
      )
    charge+=diff_hs
    if charge: verbose=0
    if verbose:
      print '  %s charge: %-2s poly_hs: %s diff_hs: %-2s' % (ag.id_str(),
                                                             charge,
                                                             poly_hs,
                                                             diff_hs,
                                                           )
    assert abs(charge)<=max_charge, 'residue %s charge %s is greater than %s' % (
      rg.atoms()[0].quote(),
      charge,
      max_charge,
    )
    if resname in allowable_amino_acid_charges:
      assert allowable_amino_acid_charges[resname]-1 <= charge <= allowable_amino_acid_charges[resname]+1, 'resname %s charge %s range %s %s' % (
        resname,
        charge,
        allowable_amino_acid_charges[resname]-1,
        allowable_amino_acid_charges[resname]+1,
        )
  else:
    restraints = _get_restraints_from_resname(ag.resname)
    ag_names = set()
    for atom in ag.atoms():
      ag_names.add(atom.name.strip())
    atom_dict = restraints.atom_dict()
    cif_names = set()
    total = 0
    for name, atom in atom_dict.items():
      cif_names.add(name)
      total += atom.partial_charge # should use formal charge!!!
    assert len(cif_names)==len(cif_names.intersection(ag_names))
    assert len(ag_names)==len(cif_names.intersection(ag_names))
    assert abs(total-int(total))<0.01, 'sum of parial charges fo %s not accurate %f' % (ag.name, total)
    charge = int(total)
    annot = 'non-polymer'
  return charge, rc, annot

def _get_restraints_from_resname(resname, mon_lib_server):
  input_resname = resname
  restraints = mon_lib_server.get_comp_comp_id_direct(resname)
  if restraints is None:
    resname = aac.three_letter_l_given_three_letter_d.get(resname, None)
    if resname is not None:
      restraints = mon_lib_server.get_comp_comp_id_direct(resname)
  if restraints is None:
    assert restraints, 'no restraints for "%s" found' % input_resname
  return restraints

def get_partial_point_charges(rg,
                              mon_lib_server,
                              hetero_charges=None):
  """
  This function relies only on the residue group and monomer library server
  """
  #assert 0
  v2_to_3 = {' HA3':' HA1',
             ' HB3':' HB1',
             ' HG3':' HG1',
             'HG13':'HG11',
             ' HD3':' HD1',
             ' HE3':' HE1',
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
            #print 'using hetero_charges for :%s' % key
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

def calculate_pdb_hierarchy_charge(hierarchy,
                                   hetero_charges=None,
                                   inter_residue_bonds=None,
                                   assert_no_alt_loc=True,
                                   list_charges=False,
                                   check=None,
                                   assert_correct_chain_terminii=True,
                                   verbose=False,
                                   ):
  charge = 0
  charges = []
  annotations = []
  if inter_residue_bonds is None: inter_residue_bonds=[]
  if assert_no_alt_loc:
    # see if we can squash into a single conf.
    hierarchy = hierarchy_utils.attempt_to_squash_alt_loc(hierarchy)
    if hierarchy is None: raise Sorry('too many alt locs to squash')
  residue_types = []
  for residue in hierarchy_utils.generate_residue_groups(
      hierarchy,
      assert_no_alt_loc=assert_no_alt_loc,
      exclude_water=True,
      ):
    validate_all_atoms(residue)
    assert len(residue.atom_groups())==1
    residue_types.append(get_class(residue.atom_groups()[0].resname))
    tmp, rc, annot = calculate_residue_charge(
      residue,
      hetero_charges=hetero_charges,
      inter_residue_bonds=inter_residue_bonds,
      assert_no_alt_loc=assert_no_alt_loc,
      verbose=verbose,
    )
    annotations.append(annot)
    if check:
      print residue
      key = 'PRO%s' % residue.parent().id
      key += '.%s' % residue.resseq.strip()
      key += '.%s' % residue.atom_groups()[0].resname
      key = key.replace('HIS', 'HSD')
      print key
      print ' CHARMM %f Phenix %f' % (check[key], tmp)
      if 1:
        print inter_residue_bonds
      assert abs(check[key]-tmp)<0.001
      assert 0
    if list_charges:
      outl = residue.id_str()
      for ag in residue.atom_groups():
        outl += '%s ' % ag.resname
      #print 'CHARGE %s %2d (%2d) %s' % (outl, tmp, rc, annot)
      charges.append([outl, tmp, annot])
    charge += tmp
    if verbose: # or display_residue_charges:
      if tmp:
        outl = '-'*80
        outl += '\nNON-ZERO CHARGE current %2d base %2d total %2d' % (tmp,rc,charge)
        for i, atom in enumerate(residue.atoms()):
          if i==0:
            outl += ' "%s"' % (atom.parent().id_str())
            if tmp!=rc: outl += ' DIFF %2d' % (tmp-rc)
          assert abs(tmp-rc)<=2, outl
          outl += "\n%s" % atom.quote()
        outl += '\n%s' % ('-'*80)
        print outl
  # check annotations
  if residue_types and assert_correct_chain_terminii:
    assert filter(None, annotations), 'No terminal or capping hydrogens found'
  if list_charges:
    #print 'CHARGE',charge
    charges.append(['Total', charge])
    return charges
  return charge

def write_pdb_hierarchy_xyzq_file(hierarchy,
                                  file_name="xyzq_cctbx.dat",
                                  hetero_charges=None,
                                  exclude_water=True,
                                  charge_scaling_positions=None,
                                  scale=0,
                                  ):
  write_charge_and_coordinates_from_hierarchy(hierarchy,
                                              file_name=file_name,
                                              qxyz_order='xyzq',
                                              hetero_charges=hetero_charges,
                                              exclude_water=exclude_water,
                                              charge_scaling_positions=None,
                                              scale=0,
                                              )
def scale_partial_point_charges(qxyz,
                                charge_scaling_positions=None,
                                scale=0):
  def partial_charge_in_charge_scaling_positions(partial_charge,
                                                 charge_scaling_positions):
    scaling = False
    for xyz in charge_scaling_positions:
      same_point = ( abs(xyz[0] - partial_charge[1]) < 1.0E-3 and
                     abs(xyz[1] - partial_charge[2]) < 1.0E-3 and
                     abs(xyz[2] - partial_charge[3]) < 1.0E-3 )
      if same_point:
        scaling = True
        break
    return scaling
  #####
  if charge_scaling_positions != None:
    for item in qxyz:
      if partial_charge_in_charge_scaling_positions(item,
                                                    charge_scaling_positions):
        item[0] =  item[0]*scale

def get_hetero_charges_FORMUL(pdb_inp):
  # get the hetero charges from the FORMUL record
  hetero_charges = {}
  for line in pdb_inp.heterogen_section():
    if line.find("FORMUL")==-1: continue
    tmp = line.split()
    hetero_charges.setdefault(tmp[2], 0)
    charge = tmp[-1].replace(")", '')
    sign = None
    if charge.find('-')>-1: sign=-1
    if charge.find('+')>-1: sign=1
    if sign:
      charge=charge.replace("+", '').replace('-','')
      charge = int(charge)
      hetero_charges[tmp[2]]=charge*sign
  return hetero_charges

def get_charge_from_restraints(resname):
  pass

def get_hetero_charges_DB(pdb_hierarchy):
  for atom_group in pdb_hierarchy.atom_groups():
    restraints = _get_restraints_from_resname(atom_group.resname)
    if restraints.is_peptide():
      print 'PEPTIDE',atom_group.resname
      continue
    print atom_group.resname
    print restraints
    print dir(restraints)
    restraints.show()
    get_charge_from_restraints(atom_group.resname)
  assert 0

def get_hetero_charges(pdb_inp, pdb_hierarchy):
  if 1:
    return get_hetero_charges_FORMUL(pdb_inp)
  else:
    return get_hetero_charges_DB(pdb_hierarchy)

def get_inter_residue_bonds(ppf, verbose=False):
  # must use this before changing the hierarchy
  # used for S-S bonds etc...
  inter_residue_bonds = {}
  grm = ppf.geometry_restraints_manager()
  if not hasattr(grm, 'get_all_bond_proxies'): return inter_residue_bonds
  atoms = ppf.all_chain_proxies.pdb_hierarchy.atoms()
  for bond in grm.get_all_bond_proxies():
    if not hasattr(bond, 'get_proxies_with_origin_id'): continue
    for p in bond.get_proxies_with_origin_id():
      assert p.origin_id==0
      r1 = atoms[p.i_seqs[0]]
      r2 = atoms[p.i_seqs[1]]
      # exclude peptide links
      # but maybe should include all for completeness
      if r1.name.strip()=='C' and r2.name.strip()=='N': continue
      #outl = 'bonding %s %s' % ( r1.quote(),r2.quote())
      r1=r1.parent().parent()
      r2=r2.parent().parent()
      if r1.id_str()!=r2.id_str():
        inter_residue_bonds[p.i_seqs] = True
        if verbose:
          print 'bonding',atoms[p.i_seqs[0]].quote(), atoms[p.i_seqs[1]].quote()
        for i in range(2):
          inter_residue_bonds.setdefault(p.i_seqs[i], [])
          inter_residue_bonds[p.i_seqs[i]].append(inter_residue_bonds[p.i_seqs])
        #assert not verbose
  return inter_residue_bonds

def run(pdb_filename,
        ligand_file_names=None,
        list_charges=False,
        assert_correct_chain_terminii=True,
        verbose=False):
  assert 0
  data = {}
  if os.path.exists(pdb_filename.replace('.pdb', '.psf')):
    f=file(pdb_filename.replace('.pdb', '.psf'), 'rb')
    lines = f.readlines()
    f.close()
    for line in lines:
      tmp = line.split()
      if not tmp: continue
      if tmp[1].find('PRO')==-1: continue
      key = '.'.join(tmp[1:4])
      data.setdefault(key, 0)
      data[key]+=float(tmp[6])

  return total_charge

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  list_charges=False
  if len(args)>1 and args[1]:
    list_charges=True
  total_charge = run(*tuple(args), list_charges=list_charges)
  print "total_charge",total_charge
