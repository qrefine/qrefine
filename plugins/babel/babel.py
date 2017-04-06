##
# this is from https://github.com/cstein/blogsamples
##

import openbabel

# quench tedious output of non-conforming pdb-format
openbabel.obErrorLog.SetOutputLevel(-1)

def OBMolFromFilename(filename):
  obmol = openbabel.OBMol()
  obconv = openbabel.OBConversion()
  obconv.SetInFormat("pdb")
  obconv.ReadFile(obmol, filename)
  return obmol

def OBMolToFilename(obmol, filename):
  obconv = openbabel.OBConversion()
  obconv.SetOutFormat("pdb")
  obconv.WriteFile(obmol, filename)

def OBAtomFromIndex(mol, index):
  return mol.GetAtom(index)

def OBSmartMatches(mol, pattern):
  """This function matches a SMARTS pattern to a molecule
  """
  obpat = openbabel.OBSmartsPattern()
  obpat.Init(pattern)
  obpat.Match(mol)
  return [m for m in obpat.GetUMapList()]

def OBMolMinimize(mol):
  """Minimize a molecule
  """
  ff = openbabel.OBForceField.FindForceField("MMFF94")
  ff.Setup(mol)
  ff.ConjugateGradients(100, 1.0e-5)
  return mol

def OBStructureFromSmiles(smilesstring, filename=None):
  mol = openbabel.OBMol()
  obConversion = openbabel.OBConversion()
  obConversion.SetInAndOutFormats("smi", "pdb")
  obConversion.ReadString(mol, smilesstring)
  mol.AddHydrogens() #False, True, 7.4)
  builder = openbabel.OBBuilder()
  builder.Build(mol)
  mol = OBMolMinimize(mol)
  if filename is None: return mol
  # save structures in subfolder molecules
  obConversion.WriteFile(mol, "molecules/%s.pdb" % filename)

def atom_at_edge_positions(edge_positions,atom):
  edge_atom = False
  for p in edge_positions:
    if abs(atom.GetX()-p[0]) < 1.0E-3 and abs(atom.GetY()-p[1]) < 1.0E-3  and abs(atom.GetZ()-p[2]) < 1.0E-3:
      edge_atom = True
  return edge_atom

def add_complete_hydrogens(file_in,file_out,edge_positions=None):
  mol = OBMolFromFilename(file_in)
  for residue in openbabel.OBResidueIter(mol):
    rname = residue.GetName()
    chain = residue.GetChain()
    for atom in openbabel.OBResidueAtomIter(residue):
      if atom.IsCarbon() or atom.IsNitrogen():
        idx = atom.GetIdx()
        imval = atom.GetImplicitValence()
        reval = atom.GetValence()
        #if imval != reval:
        if atom.IsNitrogen() and atom.IsInRing() and rname=="HIS":continue
        if edge_positions != None and not atom_at_edge_positions(edge_positions,atom):continue
        atom.SetImplicitValence(3)
        mol.AddHydrogens(atom)
  OBMolToFilename(mol, file_out)

def get_unsaturated_atom_positions(file_in,edge_positions=None):
  positions = []
  mol = OBMolFromFilename(file_in)
  for residue in openbabel.OBResidueIter(mol):
    rname = residue.GetName()
    chain = residue.GetChain()
    for atom in openbabel.OBResidueAtomIter(residue):
      if atom.IsCarbon() or atom.IsNitrogen():
        idx = atom.GetIdx()
        imval = atom.GetImplicitValence()
        reval = atom.GetValence()
       # if imval != reval:
        if atom.IsNitrogen() and atom.IsInRing() and rname=="HIS":continue
        if edge_positions != None and not atom_at_edge_positions(edge_positions,atom):continue
        xyz = []
        xyz.append(atom.GetX())
        xyz.append(atom.GetY())
        xyz.append(atom.GetZ())
        positions.append(xyz)
  return positions