import math
import sys
from string import letters

import iotbx
from mmtbx.monomer_library import server
from scitbx import matrix
from scitbx.math import dihedral_angle

from iotbx.pdb import amino_acid_codes as aac

mon_lib_server = server.server()
get_class = iotbx.pdb.common_residue_names_get_class

from utils import hierarchy_utils

def d_squared(xyz1, xyz2):
  d2 = 0
  for i in range(3):
    d2 += (xyz2[i]-xyz1[i])**2
  return d2

def get_bond_vector(a1,a2,unit=False):
  vector = []
  l = 0
  for i in range(3):
    vector.append(a1.xyz[i]-a2.xyz[i])
    l+=vector[i]**2
  if unit:
    l=math.sqrt(l)
    for i in range(3):
      vector[i] /= l
  return tuple(vector)

def construct_xyz(ba, bv,
                  aa, av,
                  da, dv,
                  period=3,
                  ):
  assert ba is not None
  assert aa is not None
  assert da is not None
  rn = matrix.col(ba.xyz)
  rca = matrix.col(aa.xyz)
  rc = matrix.col(da.xyz)
  rcca = rc -rca

  e0 = (rn - rca).normalize()
  e1 = (rcca - (rcca.dot(e0))*e0).normalize()
  e2 = e0.cross(e1)

  pi = math.pi
  alpha = math.radians(av)
  phi = math.radians(dv)

  rh_list = []
  for n in range(0, period):
    rh = rn + bv * (math.sin(alpha)*(math.cos(phi + n*2*pi/period)*e1 +
                                     math.sin(phi + n*2*pi/period)*e2) -
                    math.cos(alpha)*e0)
    rh_list.append(rh)
  return rh_list

def _add_atom_to_chain(atom, ag):
  rg = _add_atom_to_residue_group(atom, ag)
  chain = ag.parent().parent()
  tc = iotbx.pdb.hierarchy.chain()
  tc.id = chain.id
  tc.append_residue_group(rg)
  return tc

def _add_atom_to_residue_group(atom, ag):
  tag = iotbx.pdb.hierarchy.atom_group()
  tag.resname = ag.resname
  tag.append_atom(atom)
  rg = iotbx.pdb.hierarchy.residue_group()
  rg.resseq = ag.parent().resseq
  rg.append_atom_group(tag)
  for i, c in enumerate(letters):
    if c==ag.parent().parent().id:
      break
  atom.tmp = i
  return rg

def get_atoms_by_names(ag, l=None, all_or_nothing=True):
  assert l
  rc = []
  for name in l:
    atom = ag.get_atom(name)
    rc.append(atom)
  if len(l)!=len(filter(None, rc)): return None
  return rc

def add_n_terminal_hydrogens_to_atom_group(ag,
                                           use_capping_hydrogens=False,
                                           append_to_end_of_model=False,
                                           retain_original_hydrogens=True,
                                           n_ca_c=None,
                                          ):
  rc=[]
  if n_ca_c is not None:
    n, ca, c = n_ca_c
  else:
    n = ag.get_atom("N")
    if n is None: return 'no N'
    ca = ag.get_atom("CA")
    if ca is None: return 'no CA'
    c = ag.get_atom("C")
    if c is None: return 'no C'
  atom = ag.get_atom('H')
  dihedral=120.
  if atom:
    dihedral = dihedral_angle(sites=[atom.xyz,
                                     n.xyz,
                                     ca.xyz,
                                     c.xyz,
                                   ],
                              deg=True)
  if retain_original_hydrogens: pass
  else:
    if ag.get_atom("H"): # maybe needs to be smarter or actually work
      ag.remove_atom(ag.get_atom('H'))
  #if use_capping_hydrogens and 0:
  #  for i, atom in enumerate(ag.atoms()):
  #    if atom.name == ' H3 ':
  #      ag.remove_atom(i)
  #      break
  # add H1
  rh3 = construct_xyz(n, 0.9,
                      ca, 109.5,
                      c, dihedral,
                     )
  # this could be smarter
  possible = ['H', 'H1', 'H2', 'H3', 'HT1', 'HT2']
  h_count = 0
  for h in possible:
    if ag.get_atom(h): h_count+=1
  number_of_hydrogens=3
  if use_capping_hydrogens:
    number_of_hydrogens-=1
    #if ag.atoms()[0].parent().resname=='PRO':
    #  number_of_hydrogens=-1
    #  # should name the hydrogens correctly
  if h_count>=number_of_hydrogens: return []
  for i in range(0, number_of_hydrogens):
    name = " H%d " % (i+1)
    if retain_original_hydrogens:
      if i==0 and ag.get_atom('H'): continue
    if ag.get_atom(name.strip()): continue
    if ag.resname=='PRO':
      if i==0:
        continue
    atom = iotbx.pdb.hierarchy.atom()
    atom.name = name
    atom.element = "H"
    atom.xyz = rh3[i]
    atom.occ = n.occ
    atom.b = n.b
    atom.segid = ' '*4
    if append_to_end_of_model and i+1==number_of_hydrogens:
      rg = _add_atom_to_chain(atom, ag)
      rc.append(rg)
    else:
      ag.append_atom(atom)
  return rc

def add_n_terminal_hydrogens_to_residue_group(residue_group,
                                              use_capping_hydrogens=False,
                                              append_to_end_of_model=False,
                                             ):
  rc=[]
  for ag, (n, ca, c) in generate_atom_group_atom_names(residue_group,
                                                       ['N', 'CA', 'C'],
                                                       ):
    tmp = add_n_terminal_hydrogens_to_atom_group(
      ag,
      use_capping_hydrogens=use_capping_hydrogens,
      append_to_end_of_model=append_to_end_of_model,
      n_ca_c=[n,ca,c],
    )
    assert type(tmp)!=type(''), 'not string "%s" %s' % (tmp, type(tmp))
    rc += tmp
  return rc

def add_n_terminal_hydrogens(hierarchy,
                             #residue_selection=None,
                             add_to_chain_breaks=False,
                            ):
  assert 0
  # add N terminal hydrogens because Reduce only does it to resseq=1
  # needs to be alt.loc. aware for non-quantum-refine
  for chain_i, chain in enumerate(hierarchy.chains()):
    for res_i, residue_group in enumerate(chain.residue_groups()):
      if len(residue_group.atom_groups())>1: continue
      atom_group = residue_group.atom_groups()[0]
      if get_class(atom_group.resname) not in ["common_amino_acid",
                                               "modified_amino_acid",
                                             ]:
        continue
      if res_i==0: # need better switch
        add_n_terminal_hydrogens_to_atom_group(atom_group)
  hierarchy.atoms_reset_serial()
  hierarchy.atoms().reset_i_seq()
  return hierarchy

def add_c_terminal_oxygens_to_atom_group(ag,
                                         use_capping_hydrogens=False,
                                         append_to_end_of_model=False,
                                         c_ca_n=None,
                                        ):
  #
  # do we need ANISOU
  #
  rc = []
  atom_name=' OXT'
  atom_element = 'O'
  bond_length=1.231
  if use_capping_hydrogens:
    if ag.get_atom(atom_name.strip()): return []
    atom_name=" HC "
    atom_element="H"
    bond_length=1.
  if ag.get_atom(atom_name.strip()): return []
  if c_ca_n is not None:
    c, ca, n = c_ca_n
  else:
    c = ag.get_atom("C")
    if c is None: return
    ca = ag.get_atom("CA")
    if ca is None: return
    n = ag.get_atom("N")
    if n is None: return
  atom = ag.get_atom('O')
  dihedral = dihedral_angle(sites=[atom.xyz,
                                   c.xyz,
                                   ca.xyz,
                                   n.xyz,
                                 ],
                            deg=True)
  ro2 = construct_xyz(c, bond_length,
                      ca, 120.,
                      n, dihedral,
                      period=2,
                     )
  oxys = [' O  ', atom_name]
  for i in range(0,2):
    name = oxys[i]
    atom = ag.get_atom(name.strip())
    if atom:
      pass #atom.xyz = ro2[i]
    else:
      atom = iotbx.pdb.hierarchy.atom()
      atom.name = name
      atom.element = atom_element
      atom.occ = c.occ
      atom.b = c.b
      atom.segid = ' '*4
      atom.xyz = ro2[i]
      if append_to_end_of_model:
        chain = _add_atom_to_chain(atom, ag)
        rc.append(chain)
      else:
        # add the atom to the hierarchy
        ag.append_atom(atom)
  return rc

def add_c_terminal_oxygens_to_residue_group(residue_group,
                                            use_capping_hydrogens=False,
                                            append_to_end_of_model=False,
                                          ):
  rc=[]
  for ag, (c, ca, n) in generate_atom_group_atom_names(residue_group,
                                                       ['C', 'CA', 'N'],
                                                       ):
    tmp = add_c_terminal_oxygens_to_atom_group(
      ag,
      use_capping_hydrogens=use_capping_hydrogens,
      append_to_end_of_model=append_to_end_of_model,
      c_ca_n = [c, ca, n],
    )
    rc += tmp
  return rc

def add_c_terminal_oxygens(hierarchy,
                          ):
  assert 0
  for chain_i, chain in enumerate(hierarchy.chains()):
    for res_i, residue_group in enumerate(chain.residue_groups()):
      if len(residue_group.atom_groups())>1: continue
      atom_group = residue_group.atom_groups()[0]
      if get_class(atom_group.resname) not in ["common_amino_acid",
                                               "modified_amino_acid",
                                             ]:
        continue
      if capping_hydrogens:
        assert 0
      if res_i==len(chain.residue_groups())-1: # need better switch
        add_c_terminal_oxygens_to_atom_group(atom_group)
  hierarchy.atoms_reset_serial()
  hierarchy.atoms().reset_i_seq()
  return hierarchy

def _add_hydrogens_to_atom_group_using_bad(ag,
                                           atom_name,
                                           atom_element,
                                           bond_atom,
                                           angle_atom,
                                           dihedral_atom,
                                           bond_length,
                                           angle,
                                           dihedral,
                                           append_to_end_of_model=False,
                                           ):
  rc = []
  if ag.get_atom(atom_name.strip()): return []
  if type(bond_atom)==type(''):
    ba = ag.get_atom(bond_atom)
    #print bond_atom,ba.quote()
    if ba is None: return
  else: ba = bond_atom
  if type(angle_atom)==type(''):
    aa = ag.get_atom(angle_atom)
    #print angle_atom,aa.quote()
    if aa is None: return
  else: aa = angle_atom
  if type(dihedral_atom)==type(''):
    da = ag.get_atom(dihedral_atom)
    #print dihedral_atom, da.quote()
    if da is None: return
  else: da = dihedral_atom
  ro2 = construct_xyz(ba, bond_length,
                      aa, angle,
                      da, dihedral,
                      period=1,
                     )
  atom = iotbx.pdb.hierarchy.atom()
  atom.name = atom_name
  atom.element = atom_element
  atom.occ = ba.occ
  atom.b = ba.b
  # altloc???
  atom.hetero = ba.hetero
  atom.segid = ' '*4
  atom.xyz = ro2[0]
  if append_to_end_of_model:
    chain = _add_atom_to_chain(atom, ag)
    rc.append(chain)
  else:
    ag.append_atom(atom)
  return rc

def add_cys_hg_to_atom_group(ag,
                             append_to_end_of_model=False,
                            ):
  #
  # do we need ANISOU
  #
  rc = _add_hydrogens_to_atom_group_using_bad(
    ag,
    ' HG ',
    'H',
    'SG',
    'CB',
    'CA',
    1.,
    120.,
    160.,
    append_to_end_of_model=append_to_end_of_model,
   )
  return rc

def add_cys_hg_to_residue_group(rg,
                                append_to_end_of_model=False,
                               ):
  rc=[]
  for ag in rg.atom_groups():
    if ag.resname not in ['CYS']: continue
    rc += add_cys_hg_to_atom_group(
      ag,
      append_to_end_of_model=append_to_end_of_model,
    )
  return rc

def conditional_add_cys_hg_to_atom_group(geometry_restraints_manager,
                                         rg,
                                         ):
  sgs = []
  for atom in rg.atoms():
    if atom.name.strip()=='SG' and atom.parent().resname=='CYS':
      sgs.append(atom.i_seq)
  assert len(sgs) in [0, 1]
  sg_bonds = []
  if sgs:
    for bond in geometry_restraints_manager.get_all_bond_proxies():
      if not hasattr(bond, 'get_proxies_with_origin_id'): continue
      for p in bond.get_proxies_with_origin_id():
        assert p.origin_id==0
        if sgs[0] in p.i_seqs:
          sg_bonds.append(p.i_seqs)
  if len(sg_bonds)==1:
    add_cys_hg_to_residue_group(rg)

def iterate_over_threes(hierarchy,
                        geometry_restraints_manager,
                        use_capping_hydrogens=False,
                        append_to_end_of_model=False,
                        verbose=False,
                        ):
  atoms = hierarchy.atoms()
  ###
  def get_residue_group(residue):
    for atom in residue.atoms():
      atom = atoms[atom.i_seq]
      break
    return atom.parent().parent()
  ###
  additional_hydrogens=hierarchy_utils.smart_add_atoms()
  for three in hierarchy_utils.generate_protein_fragments(
    hierarchy,
    geometry_restraints_manager,
    backbone_only=False,
    use_capping_hydrogens=use_capping_hydrogens,
  ):
    if verbose: print three
    if not len(three): continue
    ptr=0
    assert three.are_linked()
    if use_capping_hydrogens:
      for i in range(len(three)):
        rg = get_residue_group(three[i])
        conditional_add_cys_hg_to_atom_group(geometry_restraints_manager,
                                              rg,
                                              )
    # check if N-term residue - FVA
    n_term_done = False
    if three[0].resname in ['FVA',
                            ]:
      n_term_done = True
      ptr+=1
      assert ptr==1, 'ptr (%d) is not 1' % ptr
    if three.start and not n_term_done:
      ptr+=1
      assert ptr==1, 'ptr (%d) is not 1' % ptr
      rg = get_residue_group(three[0])
      rc = add_n_terminal_hydrogens_to_residue_group(
        rg,
        use_capping_hydrogens=use_capping_hydrogens,
        append_to_end_of_model=append_to_end_of_model,
      )
      if rc: additional_hydrogens.append(rc)
      #hierarchy.reset_i_seq_if_necessary()
    c_term_done = False
    if three[-1].resname in ['ETA',
                             ]:
      c_term_done = True
      ptr-=1
      assert ptr==0, 'ptr (%d) is not 0' % ptr
    if three.end and not c_term_done:
      ptr-=1
      assert ptr==0, 'ptr (%d) is not 0' % ptr
      rg = get_residue_group(three[-1])
      rc = add_c_terminal_oxygens_to_residue_group(
        rg,
        use_capping_hydrogens=use_capping_hydrogens,
        append_to_end_of_model=append_to_end_of_model,
      )
      if rc: additional_hydrogens.append(rc)
      #hierarchy.reset_i_seq_if_necessary()
    else:
      pass
  return additional_hydrogens

def iterate_using_original(hierarchy,
                           geometry_restraints_manager,
                           original_hierarchy,
                           use_capping_hydrogens=False,
                           append_to_end_of_model=False,
                           ):
  slots=[]
  start=12
  assert len(original_hierarchy.models())==1
  assert len(hierarchy.models())==1
  for chain in original_hierarchy.chains():
    for org in chain.residue_groups():
      protein = True
      for atom_group in org.atom_groups():
        if(get_class(atom_group.resname) not in ["common_amino_acid",
                                                 "modified_amino_acid",
                                               ] and
           atom_group.resname not in aac.three_letter_l_given_three_letter_d):
          protein=False
          break
      if not protein:
        slots.append(False)
        continue
      org_atom1_quote = org.atoms()[0].quote()[start:]
      for rg in hierarchy.residue_groups():
        if rg.atoms()[0].quote()[start:]==org_atom1_quote:
          slots.append(rg)
          break
      else:
        slots.append(0)
    slots.append(None)

  ptr=0
  additional_hydrogens=[]
  for i in range(len(slots)):
    start=False
    end=False
    if slots[i]:
      if i==0: start=True
      elif not slots[i-1]: start=True
      if i==len(slots)-1: end=True
      elif not slots[i+1]: end=True
      # does not work for chain ends
    else: continue
    rg = slots[i]
    conditional_add_cys_hg_to_atom_group(geometry_restraints_manager,
                                         rg)
    if start:
      ptr+=1
      assert ptr==1
      if hierarchy_utils.is_n_terminal_residue(rg):
        rc = None
      else:
        rc = add_n_terminal_hydrogens_to_residue_group(
          rg,
          use_capping_hydrogens=use_capping_hydrogens,
          append_to_end_of_model=append_to_end_of_model,
        )
      if rc: additional_hydrogens.append(rc)
    if end:
      ptr-=1
      assert ptr==0
      rc = add_c_terminal_oxygens_to_residue_group(
        rg,
        use_capping_hydrogens=use_capping_hydrogens,
        append_to_end_of_model=append_to_end_of_model,
      )
      if rc: additional_hydrogens.append(rc)
    else:
      pass
  return additional_hydrogens

def add_terminal_hydrogens(
    hierarchy,
    geometry_restraints_manager,
    add_to_chain_breaks=False,
    use_capping_hydrogens=False,  # instead of terminal H
    append_to_end_of_model=False, # useful for Q|R
    #use_capping_only_on_chain_breaks=False,
    original_hierarchy=None,
    occupancy=1.,
    verbose=False,
    ):
  # add N terminal hydrogens because Reduce only does it to resseq=1
  # needs to be alt.loc. aware for non-quantum-refine
  if original_hierarchy:
    additional_hydrogens = iterate_using_original(
      hierarchy,
      geometry_restraints_manager,
      original_hierarchy,
      use_capping_hydrogens=use_capping_hydrogens,
      append_to_end_of_model=append_to_end_of_model,
      )
  else:
    additional_hydrogens=iterate_over_threes(
      hierarchy,
      geometry_restraints_manager,
      use_capping_hydrogens=use_capping_hydrogens,
      append_to_end_of_model=append_to_end_of_model,
    )

  if append_to_end_of_model and additional_hydrogens:
    tmp = []
    for group in additional_hydrogens:
      for chain in group:
        tmp.append(chain)
    _add_atoms_from_chains_to_end_of_hierarchy(hierarchy, tmp)

def _add_atoms_from_chains_to_end_of_hierarchy(hierarchy, chains):
  lookup = {}
  for chain in chains:
    lookup.setdefault(chain.id, [])
    lookup[chain.id].append(chain)
  model = hierarchy.models()[0]
  for i, chain_group in sorted(lookup.items()):
    tc = iotbx.pdb.hierarchy.chain()
    tc.id = i
    for chain in chain_group:
      for rg in chain.residue_groups():
        tc.append_residue_group(rg.detached_copy())
    model.append_chain(tc)

def _add_atoms_from_residue_groups_to_end_of_hierarchy(hierarchy, rgs):
  assert 0
  chains = {}
  for rg in rgs:
    for atom in rg.atoms():
      # this is a bad idea
      cid = atom.tmp
      if cid not in chains:
        chains[cid] = iotbx.pdb.hierarchy.chain()
        chains[cid].id = letters[cid]
      chains[cid].append_residue_group(rg)
  model = hierarchy.models()[0]
  for i, chain in sorted(chains.items()):
    model.append_chain(chain)

def remove_acid_side_chain_hydrogens(hierarchy):
  removes = {"GLU" : "HE2",
             "ASP" : "HD2",
             }
  for ag in hierarchy.atom_groups():
    r = removes.get(ag.resname, None)
    if r is None: continue
    atom = ag.get_atom(r)
    if atom:
      ag.remove_atom(atom)
  hierarchy.atoms_reset_serial()
  hierarchy.atoms().reset_i_seq()
  return hierarchy

def _eta_peptide_h(hierarchy,
                   geometry_restraints_manager,
                   verbose=False,
                   ):
  atoms = hierarchy.atoms()
  ###
  def get_residue_group(residue):
    for atom in residue.atoms():
      atom = atoms[atom.i_seq]
      break
    return atom.parent().parent()
  ###
  for three in hierarchy_utils.generate_protein_fragments(
    hierarchy,
    geometry_restraints_manager,
    backbone_only=False,
    #use_capping_hydrogens=use_capping_hydrogens,
    ):
    if len(three)==1: continue
    if three[-1].resname!='ETA': continue
    print three
    eta = get_residue_group(three[-1])
    print dir(eta)
    previous = get_residue_group(three[-2])
    print previous
    print dir(previous)
    for ag in previous.atom_groups(): # smarter?
      previous_c = ag.get_atom('C')
      previous_o = ag.get_atom('O')
    for ag in eta.atom_groups(): # needs to be conformers...
      atom_name = ' H  '
      if ag.get_atom(atom_name):
        assert 0
      else:
        for atom in ag.atoms(): print atom.format_atom_record()
        rc = _add_hydrogens_to_atom_group_using_bad(
          ag,
          atom_name,
          'H',
          'N',
          previous_c, #'CA', 
          previous_o, #'CB',
          1.,
          120.,
          180.,
          #append_to_end_of_model=append_to_end_of_model,
        )
        assert rc is not None
        print '-'*80
        for atom in ag.atoms(): print atom.format_atom_record()
  #      hierarchy.show()
  #assert 0

def generate_atom_group_atom_names(rg, names):
  '''
  Generate all alt. loc. groups of names
  '''
  atom_groups = rg.atom_groups()
  atom_altlocs = {}
  for ag in atom_groups:
    for atom in ag.atoms():
      atom_altlocs.setdefault(atom.parent().altloc, [])
      atom_altlocs[atom.parent().altloc].append(atom)
  keys = atom_altlocs.keys()
  if len(keys)>1 and '' in keys:
    for key in keys:
      if key=='': continue
      for atom in atom_altlocs['']:
        atom_altlocs[key].append(atom)
    del atom_altlocs['']
  for key, item in atom_altlocs.items():
    atoms=[]
    for name in names:
      for atom in item:
        if atom.name.strip()==name.strip():
          atoms.append(atom)
          break
      else:
        assert 0
    yield atoms[0].parent(), atoms

def _h_h2_on_N(hierarchy,
               geometry_restraints_manager,
               verbose=False,
               ):
  atoms = hierarchy.atoms()
  ###
  def get_residue_group(residue):
    for atom in residue.atoms():
      atom = atoms[atom.i_seq]
      break
    return atom.parent().parent()
  ###
  def get_atom_from_residue_group(residue, label):
    h = None
    for ag in residue.atom_groups():
      h = ag.get_atom(label)
      if h: break
    return h
  ###
  n_done = []
  for three in hierarchy_utils.generate_protein_fragments(
    hierarchy,
    geometry_restraints_manager,
    backbone_only=False,
    #use_capping_hydrogens=use_capping_hydrogens,
    ):
    if len(three)==1: continue
    for i, residue in enumerate(three):
      if not i: continue
      residue = get_residue_group(residue)
      h = get_atom_from_residue_group(residue, 'H')
      if h is None:
        for ag, (n, ca, c) in generate_atom_group_atom_names(residue,
                                                             ['N', 'CA', 'C'],
                                                            ):
          if ag.resname in ['PRO']: continue
          if n in n_done: continue
          n_done.append(n)
          dihedral = 0
          rh3 = construct_xyz(n, 0.9,
                              ca, 109.5,
                              c, dihedral,
                            )
          atom = iotbx.pdb.hierarchy.atom()
          atom.name = ' H  '
          atom.element = "H"
          atom.xyz = rh3[0]
          atom.occ = n.occ
          atom.b = n.b
          atom.segid = ' '*4
          ag.append_atom(atom)
          assert ag.resname!='PRO'

def special_case_hydrogens(hierarchy,
                           geometry_restraints_manager,
                           verbose=False,
                          ):
  for special_case in [
      #_eta_peptide_h,
      _h_h2_on_N,
                       ]:
    rc = special_case(hierarchy,
                      geometry_restraints_manager,
                      verbose=verbose,
                      )
  hierarchy.sort_atoms_in_place()

def complete_pdb_hierarchy(hierarchy,
                           geometry_restraints_manager,
                           use_capping_hydrogens=False,
                           append_to_end_of_model=False,
                           pdb_filename=None,
                           pdb_inp=None,
                           original_pdb_filename=None,
                           verbose=False,
                           debug=False,
                          ):
  for ag in hierarchy.atom_groups():
    if get_class(ag.resname) in ['common_rna_dna']:
      raise Sorry('')
  from mmtbx.building import extend_sidechains
  params=None
  original_hierarchy = None
  if use_capping_hydrogens:
    params = hierarchy_utils.get_pdb_interpretation_params()
    params.link_distance_cutoff=1.8 # avoid linking across a single missing AA
    if original_pdb_filename:
      original_pdb_inp = iotbx.pdb.input(original_pdb_filename)
      original_hierarchy = original_pdb_inp.construct_hierarchy()
  if debug:
    output = hierarchy_utils.write_hierarchy(pdb_filename,
                                             pdb_inp,
                                             hierarchy,
                                             'temp1',
                                           )
  #
  # assume model is heavy-atom complete
  #
  if not use_capping_hydrogens:
    if debug:
      ppf = hierarchy_utils.get_processed_pdb(pdb_filename=output)
    else:
      raw_records = hierarchy_utils.get_raw_records(pdb_inp, hierarchy)
      ppf = hierarchy_utils.get_processed_pdb(raw_records=raw_records,
                                              params=params,
                                            )
      sites_cart = hierarchy.atoms().extract_xyz()
      ppf.all_chain_proxies.pdb_hierarchy.atoms().set_xyz(sites_cart)
    n_changed = extend_sidechains.extend_protein_model(
      ppf.all_chain_proxies.pdb_hierarchy,
      mon_lib_server,
      add_hydrogens=False,
    )
    if debug:
      print 'number of side chains changed',n_changed
      output = hierarchy_utils.write_hierarchy(pdb_filename,
                                               pdb_inp,
                                               ppf.all_chain_proxies.pdb_hierarchy,
                                               'temp2',
                                             )
  #
  # need to use Reduce/ReadySet! to add hydrogens
  #
  if not use_capping_hydrogens:
    output = hierarchy_utils.write_hierarchy(
      pdb_filename,
      pdb_inp,
      ppf.all_chain_proxies.pdb_hierarchy,
      'readyset_input',
    )
    hierarchy = hierarchy_utils.add_hydrogens_using_ReadySet(output)
  #
  # remove side chain acid hydrogens - maybe not required since recent changes
  #
  if debug:
    ppf = hierarchy_utils.get_processed_pdb(pdb_filename=output,
                                            params=params,
                                          )
  else:
    raw_records = hierarchy_utils.get_raw_records(pdb_inp, hierarchy)
    ppf = hierarchy_utils.get_processed_pdb(raw_records=raw_records,
                                            params=params,
                                          )
    sites_cart = hierarchy.atoms().extract_xyz()
    ppf.all_chain_proxies.pdb_hierarchy.atoms().set_xyz(sites_cart)
  remove_acid_side_chain_hydrogens(ppf.all_chain_proxies.pdb_hierarchy)
  #
  # add hydrogens in special cases
  #  eg ETA
  #  eg N - H, H2
  #
  if debug:
    ppf = hierarchy_utils.get_processed_pdb(pdb_filename=output,
                                            params=params,
                                          )
  else:
    hierarchy = ppf.all_chain_proxies.pdb_hierarchy
    raw_records = hierarchy_utils.get_raw_records(pdb_inp, hierarchy)
    ppf = hierarchy_utils.get_processed_pdb(raw_records=raw_records,
                                            params=params,
                                          )
    sites_cart = hierarchy.atoms().extract_xyz()
    ppf.all_chain_proxies.pdb_hierarchy.atoms().set_xyz(sites_cart)
  special_case_hydrogens(ppf.all_chain_proxies.pdb_hierarchy,
                         ppf.geometry_restraints_manager(),
                         #use_capping_hydrogens=use_capping_hydrogens,
                         #append_to_end_of_model=append_to_end_of_model,
                         #original_hierarchy=original_hierarchy,
                         verbose=verbose,
                       )
  #
  # add terminals atoms including hydrogens and OXT - more docs here...
  #
  if debug:
    output = hierarchy_utils.write_hierarchy(
      pdb_filename,
      pdb_inp,
      ppf.all_chain_proxies.pdb_hierarchy,
      'temp3',
    )
    ppf = hierarchy_utils.get_processed_pdb(pdb_filename=output,
                                            params=params,
                                           )
  else:
    hierarchy = ppf.all_chain_proxies.pdb_hierarchy
    raw_records = hierarchy_utils.get_raw_records(pdb_inp, hierarchy)
    ppf = hierarchy_utils.get_processed_pdb(raw_records=raw_records,
                                            params=params,
                                           )
    sites_cart = hierarchy.atoms().extract_xyz()
    ppf.all_chain_proxies.pdb_hierarchy.atoms().set_xyz(sites_cart)
  add_terminal_hydrogens(ppf.all_chain_proxies.pdb_hierarchy,
                         ppf.geometry_restraints_manager(),
                         use_capping_hydrogens=use_capping_hydrogens,
                         append_to_end_of_model=append_to_end_of_model,
                         original_hierarchy=original_hierarchy,
                         verbose=verbose,
                        ) # in place
  ppf.all_chain_proxies.pdb_hierarchy.atoms().set_chemical_element_simple_if_necessary()
  ppf.all_chain_proxies.pdb_hierarchy.sort_atoms_in_place()
  #display_hierarchy_atoms(ppf.all_chain_proxies.pdb_hierarchy)
  #ppf.all_chain_proxies.pdb_hierarchy.atoms_reset_serial()
  #ppf.all_chain_proxies.pdb_hierarchy.atoms().reset_i_seq()
  return ppf

def run(pdb_filename=None,
        pdb_hierarchy=None,
        crystal_symmetry=None,
        model_completion=True,
        original_pdb_filename=None,
        append_to_end_of_model=True,
        ):
  #
  # function as be used in two main modes
  #   1. completing a model with hydrogens in a protein-like manner
  #   2. completing a cluster with hydrogens in a QM-sensible manner
  #
  if pdb_hierarchy:
    assert crystal_symmetry
    assert pdb_filename is None

  if model_completion:
    use_capping_hydrogens=False
    fname = 'complete'
  else:
    use_capping_hydrogens=True
    fname = 'capping'
    #assert 0 # model has H
  params=None
  if use_capping_hydrogens:
    params = hierarchy_utils.get_pdb_interpretation_params()
    params.link_distance_cutoff=1.8
  if pdb_hierarchy:
    raw_records = hierarchy_utils.get_raw_records(
      pdb_inp=None,
      pdb_hierarchy=pdb_hierarchy,
      crystal_symmetry=crystal_symmetry,
    )
    ppf = hierarchy_utils.get_processed_pdb(raw_records=raw_records)
    sites_cart = pdb_hierarchy.atoms().extract_xyz()
    ppf.all_chain_proxies.pdb_hierarchy.atoms().set_xyz(sites_cart)
  else:
    ppf = hierarchy_utils.get_processed_pdb(pdb_filename=pdb_filename,
                                            params=params,
                                          )
  ppf = complete_pdb_hierarchy(
    ppf.all_chain_proxies.pdb_hierarchy,
    ppf.geometry_restraints_manager(),
    use_capping_hydrogens=use_capping_hydrogens,
    append_to_end_of_model=append_to_end_of_model, # needed for clustering 
                                                   # code and Molprobity
    pdb_filename=pdb_filename,   # used just for naming of debug output
    pdb_inp=ppf.all_chain_proxies.pdb_inp, # used in get_raw_records. why?
    original_pdb_filename=original_pdb_filename, # used to define breaks in
                                                 # main chain for capping
    verbose=False,
  )
  if pdb_filename:
    output = hierarchy_utils.write_hierarchy(
      pdb_filename,
      ppf.all_chain_proxies.pdb_inp,
      ppf.all_chain_proxies.pdb_hierarchy,
      fname,
    )
  return ppf.all_chain_proxies.pdb_hierarchy

def display_hierarchy_atoms(hierarchy, n=5):
  #print '-'*80
  for i, atom in enumerate(hierarchy.atoms()):
    if i>n: break

if __name__=="__main__":
  def _fake_phil_parse(arg):
    def _boolean(s):
      if s.lower() in ['1', 'true']: return True
      elif s.lower() in ['0', 'false']: return False
      else: return s
    rc = {arg.split('=')[0] : _boolean(arg.split('=')[1])}
    return rc
  args = sys.argv[1:]
  del sys.argv[1:]
  kwds={}
  remove=[]
  for i, arg in enumerate(args):
    if arg.find('=')>-1:
      kwds.update(_fake_phil_parse(arg))
      remove.append(i)
  remove.reverse()
  for r in remove: del args[r]
  if 'test_from_clustering' in args:
    args.remove('test_from_clustering')
    ppf = hierarchy_utils.get_processed_pdb(args[0])
    sites_cart = ppf.all_chain_proxies.pdb_hierarchy.atoms().extract_xyz()
    sites_cart[0]=(4.123456789, 7.7, 1.5)
    ppf.all_chain_proxies.pdb_hierarchy.atoms().set_xyz(sites_cart)
    kwds['pdb_hierarchy'] = ppf.all_chain_proxies.pdb_hierarchy
    kwds['crystal_symmetry'] = ppf.all_chain_proxies.pdb_inp.crystal_symmetry()
    display_hierarchy_atoms(kwds['pdb_hierarchy'])
    rc = run(None, **kwds)
    #print '='*80
    display_hierarchy_atoms(rc)
    assert 0, 'FINISHED TESTING'
  #print 'run',args,kwds
  run(*tuple(args), **kwds)
