from __future__ import division
from __future__ import print_function
import iotbx.pdb
from cctbx.array_family import flex
import mmtbx.model
import sys
import time
from functools import cmp_to_key
assert 0
atom_database = {'H' : {'valence' : 1},
                 'D' : {'valence' : 1},
                 #
                 'C' : {'valence' : 4},
                 'N' : {'valence' : 3, 'lone pairs' : 1},
                 'O' : {'valence' : 2},
                 'F' : {'valence' : 1},
                 #
                 'P' : {'valence' : 3, 'lone pairs' : 1},
                 'S' : {'valence' : 2, 'lone pairs' : 2},
                 'Cl': {'valence' : 1},
                 'Ca' : {'valence' : -2, 'metal' : True},
                 #
                 'Se' : {'valence' : 2},
                 'Cu' : {'valence' : -2,
                         #'charge' : 2,
                         'metal' : True},
                 'Zn' : {'valence' : -2,
                         #'charge' : 2,
                         'metal' : True},
                 }
# via attila
# transition metals need to have a mutliplicity set
atom_database['Cu'] = {'valence': 1, 'lone pairs': 1}

def distance2(xyz1, xyz2):
  sum = 0
  for i in range(3): sum+=(xyz2[i]-xyz1[i])**2
  return sum

class atom_property(dict):
  def __init__(self):
    for element, data in atom_database.items():
      self[element] = data

  def __repr__(self):
    outl = 'atom properties\n'
    for element, data in self.items():
      outl += '  %-2s : %s\n' % (element, data)
    return outl

  def get_valence(self, element, effective=True):
    assert effective
    return self.get(element.strip(), {}).get('valence', None)

  def get_lone_pairs(self, element):
    return self.get(element.strip(), {}).get('lone pairs', 0)

  def is_metal(self, element):
    return self.get(element.strip().capitalize(), {}).get('metal', False)

class electron_distribution(dict):
  def __init__(self,
               hierarchy,
               grm,
               alternative_location_id=None,
               alternative_location_index=None,
               verbose=False,
               ):
    alternative_location_id='A'
    self.properties = atom_property()
    self.hierarchy = hierarchy
    self.atoms = self.hierarchy.atoms()
    self.grm = grm
    self.verbose=verbose
    if filter(None, hierarchy.get_conformer_indices()):
      assert (alternative_location_id is not None or
              alternative_location_index is not None)
    for atom in hierarchy.atoms():
      e = self.properties.get_valence(atom.element)
      assert e is not None, ' element %s not found' % atom.element
      self[atom.i_seq] = e
    self.form_bonds()

  def __repr__(self):
    show_all = False
    show_unpaired = True
    show_empty_bonds = True
    atoms = self.hierarchy.atoms()
    outl = 'elec. dist.\n'
    for key, electrons in self.items():
      if type(key)==type(tuple([])):
        if(show_empty_bonds and electrons==0) or show_all:
          outl += '  %s-%s : %d\n' % (atoms[key[0]].quote(),
                                      atoms[key[1]].quote(),
                                      electrons,
          )
      else:
        assert abs(electrons)<10
        if(show_unpaired and electrons) or show_all:
          outl += '  %s  : %d\n' % (atoms[key].quote(), electrons)
    return outl

  def _generate_atoms(self):
    for key, electrons in self.items():
      if type(key)==type(tuple([])): continue
      yield key

  def _generate_bonds(self):
    for key, electrons in self.items():
      if type(key)==type(tuple([])):
        yield key

  def __setitem__(self, i_seq, electrons):
    if electrons<-1:
      if self.properties.get_lone_pairs(self.atoms[i_seq].element):
        electrons+=2
    dict.__setitem__(self, i_seq, electrons)

  def _add_electron_to_bond(self, i_seqs):
    if 0:
      atoms = self.hierarchy.atoms()
      print(i_seqs, atoms[i_seqs[0]].quote(), atoms[i_seqs[1]].quote())
      print(self)
    self[i_seqs]+=1
    self[i_seqs[0]]-=1
    self[i_seqs[1]]-=1
    if 0: print(self)

  def form_bonds(self, extend_based_on_proximity=False, verbose=False):
    if self.verbose or verbose: verbose=1
    #verbose=1
    atoms = self.hierarchy.atoms()
    def _is_max_bond_valence(i_seq):
      max_valence = self.properties.get_valence(atoms[i_seq].element)
      lp = self.properties.get_lone_pairs(atoms[i_seq].element)
      if lp: max_valence += lp*2
      for bond in self._generate_bonds():
        if i_seq in bond:
          max_valence -= self[bond]
          if max_valence==0: break
      return max_valence==0
    def _can_denote_electron_to_covalent_bond(i_seq, j_seq, verbose=False):
      if verbose:
        print('processing %s %s' % (atoms[i_seq].quote(), atoms[j_seq].quote()))
      if self[i_seq]>0 and self[j_seq]>0:
        if verbose:
            print('bonding %s %s' % (atoms[i_seq].quote(), atoms[j_seq].quote()))
        return True
      elif self[i_seq]==0 and self[j_seq]==0:
        return False
      atom1 = atoms[i_seq]
      if atom1.element_is_hydrogen() and self[i_seq]==0: return False
      atom2 = atoms[j_seq]
      if atom2.element_is_hydrogen() and self[j_seq]==0: return False
      if _is_max_bond_valence(i_seq) or _is_max_bond_valence(j_seq):
        return False
      assert i_seq==atom1.i_seq
      assert j_seq==atom2.i_seq
      if atom1.element_is_hydrogen():
        hydrogen = atom1
        other = atom2
      elif atom2.element_is_hydrogen():
        hydrogen = atom2
        other = atom1
      else:
        if self.properties.get_lone_pairs(atom1.element):
          lone_pair = atom1
          other = atom2
        elif self.properties.get_lone_pairs(atom2.element):
          lone_pair = atom2
          other = atom1
        else:
          return False
        if verbose:
          print('other-lp   %s-%s' % (other.quote(), lone_pair.quote()))
        if self[other.i_seq]>0:
          return True
        return False
      if self.properties.get_lone_pairs(other.element):
        #self[other.i_seq]+=2
        if verbose: print('hydrogen-X lone pair TRUE')
        return True
      return None
    ###
    def is_metal(atom1, atom2):
      is_metal_count = [0,1][self.properties.is_metal(atom1.element)]
      is_metal_count+= [0,1][self.properties.is_metal(atom2.element)]
      if is_metal_count==2: assert 0
      return is_metal_count
    ###
    def generate_bonds_from_simple(simple,
                                   sort_on_lone_pairs=False,
                                   ):
      def _get_sum_lone_pairs(bp):
        i_seq, j_seq = bp.i_seqs
        lp1 = self.properties.get_lone_pairs(atoms[i_seq].element)
        lp2 = self.properties.get_lone_pairs(atoms[j_seq].element)
        return lp1+lp2
      def _sort_lone_pairs(bp1, bp2):
        slp1 = _get_sum_lone_pairs(bp1)
        slp2 = _get_sum_lone_pairs(bp2)
        if slp2>slp1: return -1
        return 1
      if sort_on_lone_pairs:
        l = []
        for bp in simple:
          l.append(bp)
        # l.sort(_sort_lone_pairs)
        l.sort(key=cmp_to_key(_sort_lone_pairs))
        for bp in l:
          yield bp
      else:
        assert 0
    ###
    xrs = self.hierarchy.extract_xray_structure(
      crystal_symmetry=self.grm.crystal_symmetry)
    simple, asu = self.grm.get_all_bond_proxies(sites_cart=xrs.sites_cart())
    # need to filter out H-bonds
    # look for metal coordination
    # this needs to be intergrated with GRM to get correct metal coordination
    metal_coordination = []
    tmp = {}
    for bp in simple:
      assert bp.i_seqs not in self
      i_seq, j_seq = bp.i_seqs
      assert i_seq in self
      assert j_seq in self
      atom1 = atoms[i_seq]
      atom2 = atoms[j_seq]
      if is_metal(atom1, atom2):
        tmp[distance2(atom1.xyz, atom2.xyz)] = (atom1, atom2)
    # look for single (non-metal) bonds
    for bp in simple:
      if is_metal(atoms[bp.i_seqs[0]], atoms[bp.i_seqs[1]]): continue
      assert bp.i_seqs not in self
      i_seq, j_seq = bp.i_seqs
      assert i_seq in self
      assert j_seq in self
      mc = None
      if i_seq in metal_coordination:
        mc = atoms[i_seq]
        other = atoms[j_seq]
      elif j_seq in metal_coordination:
        mc = atoms[j_seq]
        other = atoms[i_seq]
      if mc:
        if other.element_is_hydrogen():
          continue
      self[bp.i_seqs]=0
      if _can_denote_electron_to_covalent_bond(i_seq, j_seq):
        self._add_electron_to_bond(bp.i_seqs)
        if verbose: print('single: %s-%s\n%s' % (atoms[i_seq].quote(), atoms[j_seq].quote(),self))
    # look for double bonds
    for bp in generate_bonds_from_simple(simple,
                                         sort_on_lone_pairs=True,
                                         ):
      if bp.i_seqs not in self: continue
      i_seq, j_seq = bp.i_seqs
      assert i_seq in self
      assert j_seq in self
      while self[i_seq]>0 and self[j_seq]>0:
        self._add_electron_to_bond(bp.i_seqs)
        if verbose: print('double',self)
        if verbose: print('bonding 2',atoms[i_seq].quote(), atoms[j_seq].quote())
    # look for hyper-valance bonds
    for bp in simple:
      if bp.i_seqs not in self: continue
      if verbose: print('hyper',self)
      i_seq, j_seq = bp.i_seqs
      assert i_seq in self
      assert j_seq in self
      while _can_denote_electron_to_covalent_bond(i_seq,
                                                  j_seq,
                                                  verbose=verbose):
        self._add_electron_to_bond(bp.i_seqs)
    # remove HG on sulfur bridge
    self.check_sulfur_bridge()

  def check_sulfur_bridge(self, verbose=False):
    atoms = self.hierarchy.atoms()
    for i_seq in self._generate_atoms():
      for j_seq in self._generate_atoms():
        if j_seq==i_seq: break
        atom1 = atoms[i_seq]
        atom2 = atoms[j_seq]
        if self[i_seq]<0 and self[j_seq]<0:
          bond0 = bond1 = bond2 = None
          for key in self:
            if type(key)==type(tuple([])):
              if i_seq in key and j_seq in key:
                bond0 = key
              elif i_seq in key and not bond1:
                other1=list(key)
                other1.remove(i_seq)
                if atoms[other1[0]].element_is_hydrogen():
                  bond1 = key
              elif j_seq in key and not bond2:
                other2=list(key)
                other2.remove(j_seq)
                if atoms[other2[0]].element_is_hydrogen():
                  bond2 = key
          if bond0 and bond1 and bond2:
            if verbose:
              print('-'*80)
              print(bond0, bond1, bond2)
              print(atoms[bond0[0]].quote())
              print(atoms[bond0[1]].quote())
              print(atoms[bond1[0]].quote())
              print(atoms[bond1[1]].quote())
              print(atoms[bond2[0]].quote())
              print(atoms[bond2[1]].quote())
            self[bond1]-=1
            self[bond2]-=1
            self[bond1[0]]+=1
            self[bond1[1]]+=1
            self[bond2[0]]+=1
            self[bond2[1]]+=1

  def extend_based_on_proximity(self):
    # use available electrons and proximity
    # needs more care or does not need bond proxies
    if extend_based_on_proximity and 0:
      rc = self.get_possible_covalent_bonds()
      for i_seq, j_seq in rc:
        if verbose:
          print('  forming bond between %s %s' % (atoms[i_seq].quote(),
                                                  atoms[j_seq].quote()))
        assert (i_seq, j_seq) not in self
        self[(i_seq, j_seq)] = 1
        self[i_seq]-=1
        self[j_seq]-=1

  def get_possible_covalent_bonds(self):
    rc = []
    atoms = self.hierarchy.atoms()
    for i_seq in self._generate_atoms():
      if self[i_seq]<1: continue
      for j_seq in self._generate_atoms():
        if j_seq==i_seq: break
        if self[j_seq]<1: continue
        atom1 = atoms[i_seq]
        atom2 = atoms[j_seq]
        # exclude H-H
        if atom1.element_is_hydrogen() and atom2.element_is_hydrogen(): continue
        # terminal atoms on a single amino acid C..N
        if not (atom1.element_is_hydrogen() or atom2.element_is_hydrogen()):
          continue
        d2 = distance2(atoms[i_seq].xyz ,atoms[j_seq].xyz)
        if atom1.element_is_hydrogen() or atom2.element_is_hydrogen():
          if d2<1.5:
            rc.append([i_seq, j_seq])
          continue
        assert d2>9, ' %s-%s is %0.1f' % (atoms[i_seq].quote(),
                                          atoms[j_seq].quote(),
                                          d2,
                                          )
    return rc

  def validate_atomic_formal_charges(self, verbose=False):
    data = {'*'   : {'N'  : [-1,0,1],
                     'OXT': [1,0],
                     },
            'LYS' : {'NZ' : [-1]},
            'GLU' : {'OE2': [1]},
            'ASP' : {'OD2': [1]},
            }
    rc = []
    for i_seq in self._generate_atoms():
      atom = self.atoms[i_seq]
      residue_data = data.get(atom.parent().resname, {})
      residue_data.update(data['*'])
      if not residue_data:
        if self[i_seq]:
          rc.append(i_seq)
        assert 0
      else:
        if self[i_seq] in residue_data.get(atom.name.strip(), [0]): continue
      residue_data = data['*'] # needs to be only AA?
      if self[i_seq] in residue_data.get(atom.name.strip(), [0]): continue
      rc.append(i_seq)
    if verbose:
      for i_seq in rc:
        print(i_seq, self.atoms[i_seq].quote())
    return rc

  def get_total_charge(self):
    total=0
    for key, electrons in self.items():
      if type(key)==type(tuple([])): continue
      total+=electrons
    return total*-1

  def get_charged_atoms(self):
    rc = []
    atoms = self.hierarchy.atoms()
    for key, electrons in self.items():
      if type(key)==type(tuple([])): continue
      if electrons:
        rc.append([ atoms[key],electrons])
    return rc

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
