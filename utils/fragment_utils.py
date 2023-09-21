
def get_backbone_connections(pdb_hierarchy):
  backbone_connections = []
  count = 0
  for chain in pdb_hierarchy.chains():
    residues = chain.residue_groups()
    if(not is_amino_acid(residues[0])):
      count = count + len(residues)
    else:
      for i in range(len(residues)):
        residue = residues[i]
        resid = residue.resseq_as_int()
        count = count + 1
        if(len(residues) > 1 and i < len(residues) - 1):
          if(resid == residues[i + 1].resseq_as_int() - 1):
            backbone_connections.append([count, count + 1])
  return backbone_connections

def is_amino_acid(residue):
  atom_C, atom_N, atom_O, atom_CA, amino_acid = (False,) * 5
  for atom in residue.atoms():
    if(atom.name.strip() == "N"):
      atom_N = True
    if(atom.name.strip() == "CA"):
      atom_CA = True
    if(atom.name.strip() == "O"):
      atom_O = True
    if(atom.name.strip() == "C"):
      atom_C = True
  if(atom_C is True and atom_N is True and atom_O is True and atom_CA is True):
    amino_acid = True
  return amino_acid

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

  ## just for amino acids
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
  positions = filter(lambda x: x is not None, positions)
  return positions
