
"""
Backbone atom positions at the edges of a sub-hierarchy (the point charges
around a QM fragment): the charges of the atoms whose bond to the fragment
was cut are scaled (charge embedding).
"""
import iotbx.pdb

AMINO_ACID_CLASSES = ("common_amino_acid", "modified_amino_acid",
                      "d_amino_acid")

def is_amino_acid(residue):
  get_class = iotbx.pdb.common_residue_names_get_class
  for ag in residue.atom_groups():
    if(get_class(name=ag.resname) in AMINO_ACID_CLASSES): return True
  return False

def backbone_nitrogen(residue):
  if(is_amino_acid(residue) is True):
    for atom in residue.atoms():
      if(atom.name.strip() == "N"):
        return list(atom.xyz)

def backbone_nitrogen_H(residue):
  if(is_amino_acid(residue) is True):
    for atom in residue.atoms():
      if(atom.name.strip() == "H"):
        return list(atom.xyz)

def backbone_carbon(residue):
  if(is_amino_acid(residue) is True):
    for atom in residue.atoms():
      if(atom.name.strip() == "C"):
        return list(atom.xyz)

def backbone_carbon_O(residue):
  if(is_amino_acid(residue) is True):
    for atom in residue.atoms():
      if(atom.name.strip() == "O"):
        return list(atom.xyz)

def is_nterminal_residue(chain_id, residue_id, pdb_hierarchy):
  nterminal = False
  for chain in pdb_hierarchy.chains():
    if(chain.id == chain_id and chain.residue_groups()[0].resseq_as_int() == residue_id):
      nterminal = True
      break
  return nterminal

def is_cterminal_residue(chain_id, residue_id, pdb_hierarchy):
  cterminal = False
  for chain in pdb_hierarchy.chains():
    if(chain.id == chain_id and chain.residue_groups()[-1].resseq_as_int() == residue_id):
      cterminal = True
      break
  return cterminal

def get_edge_atom_positions(pdb_hierarchy, sub_pdb_hierarchy, charge_embed=False):
  positions = []
  for chain in sub_pdb_hierarchy.chains():
    residues = chain.residue_groups()
    for i in range(len(residues)):
      resid = residues[i].resseq_as_int()
      residue = residues[i]
      if(len(residues) == 1):
        if(not is_nterminal_residue(chain.id, resid, pdb_hierarchy)):
          positions.append(backbone_nitrogen(residue))
          if(charge_embed):
            positions.append(backbone_nitrogen_H(residue))
        if(not is_cterminal_residue(chain.id, resid, pdb_hierarchy)):
          positions.append(backbone_carbon(residue))
          if(charge_embed):
            positions.append(backbone_carbon_O(residue))
      if(len(residues) > 1):
        if(i == 0):
          if(not is_nterminal_residue(chain.id, resid, pdb_hierarchy)):
            positions.append(backbone_nitrogen(residue))
            if(charge_embed):
              positions.append(backbone_nitrogen_H(residue))
          if(resid != residues[i + 1].resseq_as_int() - 1):
            positions.append(backbone_carbon(residue))
            if(charge_embed):
              positions.append(backbone_carbon_O(residue))
        elif(i == len(residues) - 1):
          if not is_cterminal_residue(chain.id, resid, pdb_hierarchy):
            positions.append(backbone_carbon(residue))
            if(charge_embed):
              positions.append(backbone_carbon_O(residue))
          if(resid != residues[i - 1].resseq_as_int() + 1):
            positions.append(backbone_nitrogen(residue))
            if(charge_embed):
              positions.append(backbone_nitrogen_H(residue))
        elif(i != 0 and i != len(residues) - 1):
          if(resid != residues[i + 1].resseq_as_int() - 1):
            positions.append(backbone_carbon(residue))
            if(charge_embed):
              positions.append(backbone_carbon_O(residue))
          if(resid != residues[i - 1].resseq_as_int() + 1):
            positions.append(backbone_nitrogen(residue))
            if(charge_embed):
              positions.append(backbone_nitrogen_H(residue))
  return [x for x in positions if x is not None]
