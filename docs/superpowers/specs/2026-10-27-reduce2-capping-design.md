# Replace qrefine capping with mmtbx.hydrogens.capping; retire finalise.py and completion.py

Date: 2026-10-27. Status: approved design (review of 2026-10-27 incorporated).

## Goal

qr.refine caps the QM fragments it cuts out of the model (and the symmetry
expansion it builds) with hydrogens. Today that capping comes from
`qrefine/completion.py` through `completion.run(model_completion=False)`.
Replace it with the capping machinery of `mmtbx.reduce2`
(`mmtbx/hydrogens/capping.py`, behind `reduce2 capping=True`), refactored so
that the same code serves reduce2 and qr.refine; retire `finalise.py` and
`completion.py`; turn `qr.finalise` into a thin wrapper that runs
`mmtbx.reduce2 qmready=True`. The 38 tests of `qr.test` keep passing, and the
capping paths of qr.refine get direct tests, which they do not have now.

The old qrefine capping is an idea guide for what worked, not a reference:
where reduce2's result or a cleaner/faster structure differs from it, follow
reduce2. In particular the "input atoms first, caps appended at the end"
layout and the index slicing built on it are not kept.

## What exists today (facts established by reading and probing the code)

### qrefine side

- `completion.run(model_completion=True)` (model completion) is already
  retired: `finalise.run` raises Sorry "Use mmtbx.reduce2 with qmready=True
  flag". Only the capping mode is live.
- Capping mode, on a hierarchy that already carries H from the master:
  one H `H2` on a backbone N whose peptide partner is absent (NH2), one H
  `HC` on a backbone C whose next N is absent (aldehyde), `HG` on a free CYS
  SG; only `common_amino_acid` residues; truncated side chains, ligands,
  nucleic acids and glycans are never capped. Caps are appended at the end
  of the hierarchy as extra chains. The hierarchy is re-interpreted three
  times per call. `original_pdb_filename` is read and never used.
- Call sites: `fragment.py` three times (fragment charges at setup,
  re-capping at every gradient evaluation in
  `get_qm_file_name_and_pdb_hierarchy`, debug output in
  `write_cluster_and_fragments_pdbs`); `restraints.py` once
  (`from_expansion`, with `selection=~master`); `finalise.py`
  (qr.finalise action=capping).
- Fragment charge: `charges_class(capped hierarchy).get_total_charge()`,
  which interprets the capped hierarchy and runs `mmtbx.ligands.electrons`.
- Index assumptions on the "input atoms first" layout:
  `restraints.from_qm.target_and_gradients` takes
  `ase_gradients[:selection.count(True)]`; `restraints.from_expansion`
  takes the master as the first `model.atoms().size()` atoms of the
  completed super sphere and asserts it is unchanged; `h_diff_sel`
  identifies added atoms by (name, resname, resid, chain id). The caps of
  `from_expansion` are placed once and never moved afterwards.
- `use_reduce` (qr.refine PHIL) only reaches the retired model-completion
  mode: dead.
- `tests/regression` (MongoDB-based harness, `os.system` calls) is not run
  anywhere; its `datasets/finalise` and `tests/unit/babel_pdbs` are Open
  Babel dumps of clusters capped by the old code (standard H names, `H2`
  caps, B = 0), used only by `broken_to_fix/tst_48.py` and
  `test_reg_01_finalise.py`. Nothing in them exercises anything reduce2's
  own qmready tests do not.

### mmtbx side

- `reduce2 qmready=True` prepares a whole model for QM (H everywhere,
  neutron X-H distances, OXT and NH3+ at chain ends and breaks, clash
  check, total charge in the output name); `capping=True` runs
  `mmtbx.hydrogens.capping` as the last step. In `Program.run()` that is an
  inline sequence: `find_dangling_atoms`, Sorry on unknown residues,
  `extra_sigma_bonds` into `_CalculateTotalCharge` (electrons.run with
  `extra_sigma_bonds`, open-shell warnings), `_CapDanglingAtoms`
  (`cap_dangling_atoms`, a new unprocessed manager on the capped
  hierarchy, the Capping report section).
- `find_dangling_atoms(model)` needs a restraints manager; it examines
  non-H atoms of residues with a dictionary that miss a **heavy**
  dictionary partner. It does not examine atoms that only miss an H
  partner: protonation is reduce2's job, and the CCD's protonation differs
  from reduce2's (ASP, HIS, OXT), so a "missing H" rule would be wrong.
- `cap_dangling_atoms(model, result)` appends each H to its parent's atom
  group (so caps sit inside residues), names them (free dictionary H names,
  then H + parent name for O/N/S, then H + missing name, then HC1...),
  runs a symmetry-aware clash check. It reuses `result.bond_table` when
  present and works on an unprocessed model manager (verified: identical
  caps, 0.31 s against 0.83 s for pdb_interpretation of a 47-atom
  fragment). The capped hierarchy cannot be interpreted again.
- Only `programs/reduce2.py` and the tests call these functions.

### Probe results (on `tests/unit/data_files/m2_complete_box_large.pdb`, residues 3-5 cut from 1-7)

| case | old caps | new caps (module as is) | old charge | new charge |
|---|---|---|---|---|
| peptide cut both ends | `H2` on N3, `HC` on C5 | `HOXT` on C5 only, 0.08 A from old `HC` | 0 | -1 (N3 uncapped) |
| same, GLN 4 truncated at CD | as above (CD not capped) | `HNE2` on CD4 (index 30 of 47), `HOXT` on C5 | -1 | -1 (N3 uncapped) |

`reduce2 qmready=True capping=True` on the same fragment gives charge 0
because it builds NH3+ and OXT at the cut ends (zwitterion).

## Design decisions

### D1. Chemistry at a cut: neutral link-atom caps, defined by cut bonds (approved)

Every bond between an atom inside the fragment and an atom outside it in the
parent model becomes one H on the inside atom: NH2 and aldehyde H at a
peptide cut as today, `HG` on a cut disulfide, the natural H on a cut
glycosidic or isopeptide link, `HO3'` on a cut nucleic-acid O3'. No
protonation is guessed: the master was protonated by reduce2 and only the
atoms that lost a bond change. A cut through a phosphodiester bond leaves
the P side hypervalent and stops with the module's Sorry (approved; the old
code silently left it uncapped).

Rejected: the full `reduce2 qmready=True capping=True` per fragment (re-places
and re-optimizes every H, zwitterions at the cuts, a full reduce2 run per
fragment per gradient evaluation).

### D2. Refactor the capping sequence out of reduce2 so both callers share it

In `mmtbx/hydrogens/capping.py`:

- `find_dangling_atoms(model, guess_for_unknowns=False, log=None, cut_bonds=None)`.
  `cut_bonds` is `{i_seq: n}`: bonds of the model's atoms to atoms that are
  not in the model. An atom with no dictionary record gets one with
  `n_fill = n`, `coordination = bonded neighbours + n`, `missing_names = []`,
  `h_names` = its free dictionary H names (caps are named `H2` on a
  backbone N, `H` on a PRO N, `HG` on a cut SG, `HO3'` on a cut O3'). An
  atom that already has a dictionary record keeps it: its `n_fill` already
  counts the absent partner (cut backbone C: valence 4, used CA + O double
  = 3, fill 1). H and metals in `cut_bonds` are ignored. Downstream
  (`extra_sigma_bonds`, `cap_dangling_atoms`, hypervalent and internal-gap
  checks) is unchanged. Module docstring gets rule 9, "cut bonds".
- New `cap_model(model, guess_for_unknowns=False, cut_bonds=None,
  use_neutron_distances=True, calculate_charge=True, dangling=None, log=None)`:
  the sequence reduce2 runs inline today. `dangling` may be given (a cached
  `find_dangling_atoms` result for a model with the same atoms, placed on
  the model's current sites), else it is computed. Sorry on
  `unknown_residues` when not guessing (the text reduce2 prints now).
  Charge: `electrons.run(model=model, extra_sigma_bonds=...,
  return_formal_charges=True)` before anything is placed, with the
  open-shell warnings; `None` when `calculate_charge=False` or there is
  nothing dangling and the caller does not want it. Then
  `cap_dangling_atoms` on the model's hierarchy (in place, as now).
  Returns `group_args(dangling, capped, total_charge, warnings)` with
  `capped` the `cap_dangling_atoms` result (added, n_added, n_atoms,
  n_residues, lines) and `warnings` the open-shell lines.
- `reduce2.Program.run()` calls `cap_model` in place of its inline sequence
  and keeps its own report formatting, output naming and the unprocessed
  output manager (`_CalculateTotalCharge` and `_CapDanglingAtoms` shrink to
  the reduce2-specific parts). reduce2's behaviour and output do not
  change: `tst_reduce2_capping.py`, `tst_reduce2_qmready_charge.py`,
  `tst_reduce2_qmready_clash.py` rerun unchanged.

### D3. qrefine side: `qrefine/fragment_capping.py` replaces completion.py

```python
class capper(object):
  """Caps the QM fragment parent_hierarchy.select(selection) the way
  reduce2 capping=True does, plus one H per bond the cut severed.

  __init__(parent_hierarchy, selection, crystal_symmetry, cif_objects=None, log=None):
    fragment = parent.select(selection) (fragment order kept); cut bonds by
    distance across the selection boundary (heavy, non-metal pairs closer
    than the sum of covalent radii + 0.4 A, cctbx.eltbx.
    distance_based_connectivity on the fragment plus the parent atoms
    within 3 A of it); model_from_hierarchy(fragment, cs, cif_objects)
    (processed once); capping.find_dangling_atoms(model, cut_bonds=...)
    cached; capping.cap_model(..., dangling=cached) for the charge and the
    first caps.
  Attributes: charge (int), n_caps, lines (report), dangling.
  cap(sites_cart=None) -> group_args(pdb_hierarchy, fragment_selection):
    re-places the caps on the fragment at sites_cart (the fragment's
    current coordinates, in fragment order; None: the setup coordinates):
    an unprocessed manager on a deep copy of the fragment hierarchy at
    those sites, capping.cap_model(dangling=cached, calculate_charge=False).
    pdb_hierarchy: the capped fragment, caps inside their residues as
    reduce2 leaves them; fragment_selection: flex.bool over its atoms,
    True for the fragment atoms (False for the caps), so the fragment atom
    i is at fragment_selection.iselection()[i].
  Raises Sorry (at __init__ or cap): residue without a dictionary (supply
  its CIF), internal gap, hypervalent P or S, cap clash, fragment atom
  without a bonded neighbour."""
```

One distance-based cut detection serves both uses, and the parent is always
the super cell (`SCM.super_cell_hierarchy`): for the expansion with
`SCM.keep_selection`, for a fragment with its selection on the sphere
translated through `keep_selection` (`fragments.make_capper`). The sphere
cannot be the fragments' parent: a symmetry-copy residue at the sphere's
residue-complete edge has its peptide partner outside the sphere, which is
a cut bond, not an unsatisfied valence (found in the final review on 1bdw;
`tst_38.test_sphere_boundary_fragment`). Heavy atom pairs across the
boundary closer than ~1.9 A (C/N/O) or ~2.5 A (S-S) are covalent in any
usable model; an unselected atom of a residue group with selected atoms is
an alternate conformer, not a partner. A fragment may carry one altloc
label (qrefine's per-conformer hierarchies do); two conformers are refused.

`model_from_hierarchy` moves from `restraints.py` to
`utils/hierarchy_utils.py` (restraints.py imports it from there).

### D4. Layout: caps inside residues, consumers use the mask (approved refactor)

- `restraints.from_qm.target_and_gradients` (clustering): takes the capped
  hierarchy and `fragment_selection` from the capper and does
  `gradients = ase_gradients.select(fragment_selection)`; the written QM
  PDB is the capped hierarchy as is.
- `fragment.py`: `fragments` builds one `capper` per fragment at setup
  (`get_fragment_hierarchies_and_charges`: charge from `capper.charge`,
  `fragment_capped_initial` from `capper.cap()`), stores the cappers in
  `fragment_extracts` (`fragment_cappers`), and
  `get_qm_file_name_and_pdb_hierarchy` calls
  `capper.cap(sites_cart=current fragment sites)` instead of re-interpreting
  the fragment; it returns the file name, the hierarchy and the selection.
  `write_cluster_and_fragments_pdbs` uses the same cappers.
  `original_pdb_filename` plumbing goes (`expansion_file` stays: it is
  written for debugging).
- `restraints.from_expansion`: one `capper` on the sphere (parent super
  cell, selection `keep_selection`); `pdb_hierarchy_super_completed` is
  `capper.cap().pdb_hierarchy`; `self.selection` (the master in the
  completed hierarchy) becomes the mask of the first `model.atoms().size()`
  fragment atoms through `fragment_selection.iselection()`; the master
  assertion (`is_similar_hierarchy`) uses it. `_update` re-caps on the new
  geometry (`capper.cap(sites_cart=sphere sites)`) and copies all sites
  into `pdb_hierarchy_super_completed` (same atom content every time, so
  the `from_qm` restraints manager built on it stays valid); `h_diff_sel`
  is no longer needed there. The old `selection=~master` goes: a complete
  master has no dangling atoms (qr.refine requires a model prepared with
  `reduce2 qmready=True`; a model capped by `reduce2 capping=True` cannot
  be qr.refine input at all because its cap names stop
  pdb_interpretation).
- Charge: `from_expansion` passes `capper.charge` to
  `restraints.update(..., charge=)`, which hands it to `from_qm(charge=)`
  unless the user set `quantum.charge`. `from_qm` with `clustering=False`
  and `charge=None` still computes the charge of the (uncapped) master with
  `charges_class`. No capped hierarchy is interpreted anywhere.

### D5. qr.finalise: pass-through wrapper (approved)

`command_line/finalise.py` keeps the `qr.finalise` dispatcher and becomes:
print a note ("qr.finalise now runs mmtbx.reduce2 qmready=True; the former
options action, calculate_charge, reduce, options.neutron,
hydrogen_atom_occupancies, append_to_end_of_model, skip_validation,
stop_if_poor_model are gone; the charge is always calculated and put into
the output name; run mmtbx.reduce2 directly for other options"), then
`iotbx.cli_parser.run_program(reduce2.Program, args=args + ['qmready=True'])`.
Other arguments pass through (`capping=True`, `approach=remove` work).
Output names are reduce2's (`<prefix>H_charge<q>.pdb`). Old option names
are not translated; reduce2's own "not recognized" Sorry plus the note
tell the user what changed.

### D6. Removals (approved, including the regression dataset)

- `qrefine/finalise.py`, `qrefine/completion.py`.
- qr.refine PHIL `use_reduce` (definition, the aimnet2 auto-set, the log
  line; grep docs).
- `utils/hierarchy_utils.py` functions used only by the removed modules
  (verify by grep: `write_hierarchy`, `add_hydrogens_using_reduce`,
  `add_hydrogens_using_ReadySet`, `smart_add_atoms`,
  `get_pdb_interpretation_params`, `remove_atom_from_chain`,
  `is_n_terminal_residue`).
- `tests/unit/batch_run_finalise.py`, `tests/unit/broken_to_fix/tst_48.py`,
  `tests/unit/broken_to_fix/tst_49.py`, `tests/unit/babel_pdbs/`,
  `tests/regression/test_reg_01_finalise.py` with its registration in
  `regression_tests.py` and `qr_reg_data_finalise` in `test_reg_00_base.py`,
  `tests/regression/datasets/finalise/`.

### D7. Behaviour changes (accepted)

- Cap names: `HOXT` instead of `HC` on a cut carbonyl C (`H2` on N
  unchanged). `charges.py` `c_capping` gets `HOXT` next to `HC` for its
  non-electrons path.
- Cap geometry from capping.py's builder (0.08 A for the C cap in the
  probe); the clustered-vs-whole gradient tests (`tst_35`, `tst_36`, wdelta
  tolerance 10) confirm.
- Truncated side chains, cut disulfides, cut glycans and ligands with a
  dictionary are now capped (were silently left radical or charged).
- New Sorry conditions: residue without a dictionary in a fragment,
  internal gap (TYR without CE1), phosphodiester cut, cap clash.
- Caps in `from_expansion` follow the geometry (were stale).
- Capped fragment files have caps inside residues, not in extra chains.

## Testing

mmtbx (`mmtbx/hydrogens/tst_capping.py`, in `mmtbx/run_tests.py`):
`cut_bonds` on an N (`H2`, tetrahedral), a PRO N (`H`), a cut SG (`HG`), an
atom with a dictionary record (C with cut N: n_fill stays 1), H and metal
entries (ignored), `extra_sigma_bonds` counts them, `cap_model` returns the
charge of the capped molecule and the same caps as the two-step calls, a
cached `dangling` on an unprocessed model gives identical caps, and
`cut_bonds=None` is unchanged. Rerun `tst_reduce2_capping.py`,
`tst_reduce2_qmready_charge.py`, `tst_reduce2_qmready_clash.py`.

qrefine, new `tests/unit/tst_38.py` (no QM engine, seconds):
- peptide cut both ends of `m2_complete_box_large.pdb` 3-5: `H2` on N3 and
  `HOXT` on C5, `fragment_selection` picks exactly the fragment atoms in
  order, charge 0, X-H lengths and angles sane;
- same with GLN 4 truncated: `HNE2` too, charge 0;
- `cap(sites_cart=shifted)` re-places caps on the moved fragment with the
  same names and count, without re-interpretation;
- a PRO cut and a cut disulfide from `2lvr.pdb` (confirm it has one, else
  another test file);
- the whole model selected: no caps, charge equal to `qr.charges`;
- a residue without a dictionary: Sorry.

Rewritten: `tst_18` (cap a fragment cut from its expansion model, check
caps and selection); `tst_20` (`qr.finalise` on the Zn/CYS model: reduce2
output without `HG` on the Zn-bound CYS and with ARG `HE`; if reduce2's
qmready refuses that input, adapt the input).

Existing end-to-end checks: `restraints=qm` tests with clustering or
expansion (`tst_04`, `tst_12`, `tst_24`-`tst_26`, `tst_29`, `tst_32`-`tst_36`),
then the full `qr.test`.

Docs: `rst_files/qr/finalise.txt`, the workflow step in `rst_files/qr.txt`,
`rst_files/qr/tests.txt` rows for `tst_18`, `tst_20`, `tst_38`, `CHANGES.md`.

## Implementation phases

1. mmtbx: `cut_bonds` and `cap_model` in capping.py, tests; reduce2 `run()`
   switched to `cap_model`; reduce2 capping tests rerun.
2. qrefine: `fragment_capping.capper` + `tst_38.py` (TDD);
   `model_from_hierarchy` moved.
3. `fragment.py` and `restraints.py` switched to the capper and the mask
   (`update(charge=)`, re-capping in `_update`); QM cluster/expansion tests.
4. `qr.finalise` wrapper; `tst_18`, `tst_20` rewritten.
5. Removals (D6), docs, CHANGES.md.
6. Full `qr.test` and the mmtbx tests touched.
