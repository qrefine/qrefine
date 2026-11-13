from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os, glob
import iotbx.pdb
import iotbx.cif
import mmtbx.model
from iotbx.pdb import hierarchy
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.utils import null_out
from mmtbx.hydrogens.tst_xh_mapper import pdb_str
from qrefine.tests.unit import run_qrefine_in, run_fmodel
from qrefine.tests.unit import fast_args, fast_args_no_weight_search

# Identifiers that do not fit the PDB format: four-character residue name
# ETHA (ethanol; restraints: chem_data geostd EOH renamed), plus
# two-character chain IDs (these alone would still fit).
chain_renames = {"A": "AA", "S": "SW"}
ligand_chain  = "LG"
ligand_resname = "ETHA"


def etha_cif(file_name):
  """
  Copy of chem_data/geostd/e/data_EOH.cif with every EOH replaced by ETHA,
  written to file_name. Returns [(atom name, element, (x,y,z))] of the
  ideal coordinates in _chem_comp_atom.
  """
  import libtbx.load_env
  src = libtbx.env.find_in_repositories(
    relative_path="chem_data/geostd/e/data_EOH.cif", test=os.path.isfile)
  assert src is not None
  text = open(src).read().replace("EOH", ligand_resname)
  assert text.find("EOH") == -1
  with open(file_name, "w") as fo:
    fo.write(text)
  block = iotbx.cif.reader(input_string=text).model()[
    "comp_%s" % ligand_resname]
  loop = block.get_loop("_chem_comp_atom")
  result = []
  for name, element, x, y, z in zip(loop["_chem_comp_atom.atom_id"],
      loop["_chem_comp_atom.type_symbol"], loop["_chem_comp_atom.x"],
      loop["_chem_comp_atom.y"], loop["_chem_comp_atom.z"]):
    result.append((name, element, (float(x), float(y), float(z))))
  assert len(result) == 9, result
  return result

def make_model(file_name, ligand_atoms):
  """
  The tst_xh_mapper fixture with chains renamed (A->AA, S->SW) and an
  ethanol ETHA (chain LG) whose O is 6 A from SER OG (altloc A) along -x,
  written with write_pdb_or_mmcif_file. Returns (the written file name,
  the number of atoms).
  """
  m = mmtbx.model.manager(
    model_input = iotbx.pdb.input(source_info=None, lines=pdb_str),
    log = null_out())
  h = m.get_hierarchy()
  for chain in h.chains():
    chain.id = chain_renames[chain.id]
  og = [a for a in h.atoms() if a.name.strip()=="OG" and
    a.parent().altloc=="A"]
  assert len(og) == 1
  xyz_o = [xyz for name, e, xyz in ligand_atoms if name == "O"][0]
  shift = matrix.col(og[0].xyz) + matrix.col((-6,0,0)) - matrix.col(xyz_o)
  ag = hierarchy.atom_group(resname=ligand_resname, altloc="")
  for name, element, xyz in ligand_atoms:
    a = hierarchy.atom()
    a.set_name(" %-3s" % name)
    a.set_element(element.rjust(2))
    a.set_xyz(tuple(matrix.col(xyz) + shift))
    a.set_occ(1.0)
    a.set_b(10.0)
    a.set_hetero(True)
    ag.append_atom(a)
  rg = hierarchy.residue_group(resseq="   1", icode=" ")
  rg.append_atom_group(ag)
  ch = hierarchy.chain(id=ligand_chain)
  ch.append_residue_group(rg)
  h.models()[0].append_chain(ch)
  h.atoms_reset_serial()
  # The point of the test: this model cannot be written in PDB format.
  assert not h.fits_in_pdb_format()
  # Nearest protein/water atom to the ligand.
  sel_lig = flex.bool([a.parent().resname==ligand_resname for a in h.atoms()])
  xyz = h.atoms().extract_xyz()
  d_min = min([(matrix.col(xyz[i])-matrix.col(xyz[j])).length()
    for i in sel_lig.iselection() for j in (~sel_lig).iselection()])
  assert d_min > 3.5, d_min
  result = h.write_pdb_or_mmcif_file(target_filename=file_name,
    target_format="pdb", crystal_symmetry=m.crystal_symmetry())
  assert result == os.path.splitext(file_name)[0]+".cif", result
  return result, h.atoms_size()

def identifiers(model):
  h = model.get_hierarchy()
  return set([c.id for c in h.chains()]), set(h.overall_counts().resnames)

def check_output(out, n_atoms):
  """Output is mmCIF with the identifiers intact and all the atoms."""
  assert out.pdb is not None
  assert out.pdb.endswith("_refined.cif"), out.pdb
  assert glob.glob(os.path.join(out.folder, "pdb", "*.pdb")) == []
  chains, resnames = identifiers(out.model)
  assert chains == set(list(chain_renames.values())+[ligand_chain]), chains
  assert ligand_resname in resnames, resnames
  assert out.model.get_hierarchy().atoms_size() == n_atoms

def n_atoms_in(file_name):
  """Number of atoms in a model file (PDB or mmCIF)."""
  return iotbx.pdb.input(file_name=file_name).atoms().size()

def n_atoms_in_mop(file_name):
  """Number of atoms in a MOPAC input (lines after the 3-line header)."""
  lines = open(file_name).read().splitlines()[3:]
  return len([l for l in lines if l.strip() and not l.startswith("Tv")])

def models_under(folder):
  """Model files under folder (recursively), split by extension."""
  pdb, cif = [], []
  for root, dirs, files in os.walk(folder):
    for f in files:
      if(f.endswith(".pdb")): pdb.append(os.path.join(root, f))
      if(f.endswith(".cif")): cif.append(os.path.join(root, f))
  return sorted(pdb), sorted(cif)

def has_ligand(file_name):
  h = iotbx.pdb.input(file_name=file_name).construct_hierarchy()
  return ligand_resname in h.overall_counts().resnames

def check_scratch(out, label):
  """
  Every model file under ase/ is non-empty and is mmCIF exactly when it
  holds ETHA (fragments of protein and water alone still fit PDB format);
  each MOPAC input <label>.mop has as many atoms as the model next to it
  (the hierarchy the engine received). Returns the base names of the
  mmCIF files.
  """
  pdb, cif = models_under(os.path.join(out.folder, "ase"))
  assert len(cif) > 0
  for f in pdb + cif:
    assert n_atoms_in(f) > 0, f
    assert has_ligand(f) == f.endswith(".cif"), f
  mop = sorted(glob.glob(os.path.join(out.folder, "ase", "**", label+".mop"),
    recursive=True))
  assert len(mop) > 0, label
  n_cif_mop = 0
  for f in mop:
    model_file = [f[:-4]+e for e in [".pdb", ".cif"]
      if os.path.isfile(f[:-4]+e)]
    assert len(model_file) == 1, (f, model_file)
    model_file = model_file[0]
    n_cif_mop += model_file.endswith(".cif")
    n_mop, n_model = n_atoms_in_mop(f), n_atoms_in(model_file)
    print(model_file, "atoms in .mop:", n_mop, "in the model:", n_model)
    assert n_mop == n_model > 0, (f, n_mop, n_model)
  assert n_cif_mop > 0
  return [os.path.basename(f) for f in cif]

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  qr.refine on a model whose identifiers do not fit the PDB format (residue
  name ETHA, two-character chain IDs) with cctbx restraints, with MOPAC
  without clustering (expansion=True) and with clustering: the refined
  model and the cycle models are pdb/*.cif with the chain IDs and ETHA
  intact and all the atoms; with MOPAC the scratch, expansion and fragment
  models under ase/ are non-empty, mmCIF exactly when they hold ETHA, and
  each MOPAC input has the atoms of the model next to it.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  #
  ligand_atoms = etha_cif("etha.cif")
  model_file, n_atoms = make_model("model.pdb", ligand_atoms)
  r = run_fmodel(args=[model_file, "etha.cif", "high_res=2.5", "r_free=0.05"],
    prefix="data")
  inputs = [model_file, r.mtz, "etha.cif"]
  # cctbx
  out = run_qrefine_in(args=inputs + fast_args + ["restraints=cctbx",
    "expansion=False", "clustering=False"], prefix=prefix, folder="cctbx")
  check_output(out, n_atoms)
  # MOPAC, no clustering: the scratch model ase/.../tmp_ase.cif
  out = run_qrefine_in(args=inputs + fast_args_no_weight_search + [
    "restraints=qm", "quantum.engine_name=mopac", "clustering=False",
    "expansion=True"], prefix=prefix, folder="mopac")
  check_output(out, n_atoms)
  names = check_scratch(out, label="tmp_ase")
  assert "tmp_ase.cif" in names, names
  # MOPAC, clustering: fragment, capped fragment, expansion files
  out = run_qrefine_in(args=inputs + fast_args_no_weight_search + [
    "restraints=qm", "quantum.engine_name=mopac", "clustering=True",
    "expansion=False", "maxnum_residues_in_cluster=2"], prefix=prefix,
    folder="mopac_cluster")
  check_output(out, n_atoms)
  names = check_scratch(out, label="*_capping")
  assert "expansion.cif" in names, names
  assert [n for n in names if n.endswith("_capping.cif")], names
  assert [n for n in names if n[:-4].isdigit()], names

if(__name__ == "__main__"):
  run()
  print("OK")
