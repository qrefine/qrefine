from __future__ import division
import iotbx.pdb
from cctbx.array_family import flex
import mmtbx.model
import sys
import time

atom_database = {'H' : {'valence' : 1},
                 #
                 'C' : {'valence' : 4},
                 'N' : {'valence' : 3, 'lone pairs' : 1},
                 'O' : {'valence' : 2},
                 #
                 'P' : {'valence' : 3, 'lone pairs' : 1},
                 'S' : {'valence' : 2, 'lone pairs' : 2},
                 #'CL': {'valence' : 1},
                 }

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
    return self.get(element.strip(), {}).get('lone pairs', None)

class electron_distribution(dict):
  def __init__(self,
               hierarchy,
               grm,
               alternative_location_id=None,
               alternative_location_index=None,
               ):
    alternative_location_id='A'
    self.properties = atom_property()
    self.hierarchy = hierarchy
    self.atoms = self.hierarchy.atoms()
    self.grm = grm
    if filter(None, hierarchy.get_conformer_indices()):
      assert (alternative_location_id is not None or
              alternative_location_index is not None)
    for atom in hierarchy.atoms():
      e = self.properties.get_valence(atom.element)
      assert e is not None, ' element %s not found' % atom.element
      self[atom.i_seq]   = e
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

  def form_bonds(self, extend_based_on_proximity=False, verbose=False):
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
        print 'processing %s %s' % (atoms[i_seq].quote(), atoms[j_seq].quote())
      if self[i_seq]>0 and self[j_seq]>0:
        if verbose:
            print 'bonding %s %s' % (atoms[i_seq].quote(), atoms[j_seq].quote())
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
          print 'other-lp   %s-%s' % (other.quote(), lone_pair.quote())
        if self[other.i_seq]>0:
          return True
        return False
      if self.properties.get_lone_pairs(other.element):
        #self[other.i_seq]+=2
        if verbose: print 'hydrogen-X lone pair TRUE' 
        return True
      return None
    ###
    xrs = self.hierarchy.extract_xray_structure()
    simple, asu = self.grm.get_all_bond_proxies(sites_cart=xrs.sites_cart())
    # need to filter out H-bonds and metal coordination ...
    for bp in simple:
      if verbose: print 'single',self
      assert bp.i_seqs not in self
      self[bp.i_seqs]=0
      i_seq, j_seq = bp.i_seqs
      assert i_seq in self
      assert j_seq in self
      if _can_denote_electron_to_covalent_bond(i_seq, j_seq):
        self[bp.i_seqs]+=1
        self[i_seq]-=1
        self[j_seq]-=1
    # look for double bonds
    for bp in simple:
      if verbose: print 'double',self
      i_seq, j_seq = bp.i_seqs
      assert i_seq in self
      assert j_seq in self
      while self[i_seq]>0 and self[j_seq]>0:
        self[bp.i_seqs]+=1
        self[i_seq]-=1
        self[j_seq]-=1
    # look for hyper-valance bonds
    for bp in simple:
      if verbose: print 'hyper',self
      i_seq, j_seq = bp.i_seqs
      assert i_seq in self
      assert j_seq in self
      while _can_denote_electron_to_covalent_bond(i_seq,
                                                  j_seq,
                                                  verbose=verbose):
        self[bp.i_seqs]+=1
        self[i_seq]-=1
        self[j_seq]-=1
    # use available electrons and proximity
    if extend_based_on_proximity or 1:
      rc = self.get_possible_covalent_bonds()
      for i_seq, j_seq in rc:
        if verbose:
          print '  forming bond between %s %s' % (atoms[i_seq].quote(),
                                                  atoms[j_seq].quote())
        assert (i_seq, j_seq) not in self
        self[(i_seq, j_seq)] = 1
        self[i_seq]-=1
        self[j_seq]-=1

  def get_possible_covalent_bonds(self):
    def distance2(xyz1, xyz2):
      sum = 0
      for i in range(3): sum+=(xyz2[i]-xyz1[i])**2
      return sum
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
        print i_seq, self.atoms[i_seq].quote()
    return rc

  def get_total_charge(self):
    total=0
    for key, electrons in self.items():
      if type(key)==type(tuple([])): continue
      total+=electrons
    return total*-1

def run(pdb_filename=None,
        raw_records=None,
        ):
  if pdb_filename:
    # Read file into pdb_input class
    inp = iotbx.pdb.input(file_name=pdb_filename)
  elif raw_records:
    inp = iotbx.pdb.input(lines=raw_records, source_info='lines from PDB')
  else:
    assert 0

  # create a model manager
  model = mmtbx.model.manager(
      model_input = inp)
  # get xray structure
  xrs = model.get_xray_structure()
  grm = model.get_restraints_manager()
  t0=time.time()
  atom_valences = electron_distribution(model.get_hierarchy(), # needs to be altloc free
                                        model.get_restraints_manager().geometry,
                                        )
  print atom_valences
  total_charge = atom_valences.get_total_charge()
  print 'total_charge',total_charge
  print 'time %0.1f' % (time.time()-t0)
  rc = atom_valences.validate_atomic_formal_charges()
  return total_charge

  # get number of atoms in the input model
  n_atoms = model.get_number_of_atoms()

  # extract atom coordinates
  old_sites_cart = model.get_sites_cart()
  # generate random additions
  random_addition = flex.vec3_double(
    flex.random_double(size=n_atoms*3)-0.5)
  # actually add them to old coordinates
  new_xyz = old_sites_cart + random_addition

  # Update coordinates in model manager
  model.set_sites_cart(sites_cart=new_xyz)


  # reset B-factors (min=1, max=20)
  # generate array of new B-factors
  new_b = flex.random_double(size=n_atoms, factor=19) + 1
  # set them in xray structure
  xrs.set_b_iso(values=new_b)
  # update model manager with this xray structure
  model.set_xray_structure(xrs)
  # output result in PDB format to the screen
  print model.model_as_pdb()
  print "END"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
