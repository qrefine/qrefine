# Capping H as first-class atoms in cctbx: design

Date: 2026-11-04. Status: approved 2026-11-04 (D1-D7 settled). Repositories: cctbx_project
(iotbx, mmtbx), qrefine (tests only), phenix_regression (one test).

## 1. Goal

A model written by `mmtbx.reduce2 qmready=true capping=true` must load in
every Phenix program that builds restraints, qr.refine and phenix.refine
first. Capping H must be identifiable at any later time, in any program,
without information outside the model file, so that they can be stripped
before deposition. One mechanism, decided once.

Today such a model stops at `pdb_interpretation`: the caps carry names no
dictionary has, so they get no nonbonded energy type and the default
`stop_for_unknowns` ends the run ("Number of atoms with unknown nonbonded
energy type symbols"). reduce2 avoids this only because capping is its last
step and the manager it returns is unprocessed.

## 2. Findings

- A cap is an ordinary H in every respect but its name: one H bonded to one
  heavy atom, with the parent's occupancy, B, segid and hetero flag
  (`capping.cap_dangling_atoms`).
- Caps can carry dictionary names. `cap_names` prefers the parent's free
  dictionary H names (HB2, HB3 on a truncated CB), then 'H' + parent name
  for O, N, S, then 'H' + missing heavy atom name, then HC1, HC2. It avoids
  only names of atoms present in the conformer. So a C3' whose O3' is
  missing gets the cap HO3', a dictionary name bonded to the absent O3':
  pdb_interpretation accepts the atom, makes no bond, and the cap floats.
- `pdb_interpretation` classifies every atom of a residue against its
  dictionary (`monomer_mapping.expected_atoms` / `unexpected_atoms`) as part
  of normal processing. Candidates for caps are a by-product.
- An unknown-type rescue hook exists after proxies are built and before the
  unknown-type stop (`pdb_interpretation.py` ~6227, the experimental
  pH-dependent and enol-peptide code). The capping rescue goes at that point,
  before those two, and does not depend on them.
- Every proxy carries an origin id (`cctbx/geometry_restraints/
  auto_linking_types.py`, `linking_class`). Geometry statistics and
  validation filter on the "covalent geometry" id, 0
  (`manager.get_covalent_bond_proxies`, `validation/restraints.py`), and
  print the other ids separately. Proxies with a new "capping" id are left
  out of the reported rmsd automatically.
- Consumers that must still see caps as bonded H read all bond proxies
  regardless of origin id: the charge calculation (`electrons.py`,
  `get_all_bond_proxies`), riding H connectivity
  (`hydrogens/connectivity.py`, `get_all_bond_proxies`), and capping's own
  `bond_table` (all ids but metal coordination). So the charge of a loaded
  capped model is right, caps can ride, and re-capping a capped model finds
  its parents saturated and adds nothing.
- `iotbx.pdb.atom_selection.cache` holds the selection keywords (`water`,
  `protein`, `single_atom_residue`, ...), each a lazily computed slot. It has
  no dictionary access, but `iotbx/pdb/hierarchy.py` already imports mmtbx
  lazily inside a method, so a keyword may do the same. `phenix.pdbtools`
  builds this cache directly from the hierarchy, not through the model.
- qrefine's fragment caps (`fragment_capping.py`) are placed on hierarchies
  that are never interpreted afterwards. Nothing here touches them.

## 3. Definition

A **capping H** is an H or D atom such that

1. its residue has a restraint dictionary;
2. the dictionary cannot attach it to any present atom: its name is not in
   the dictionary, or it is but the dictionary bond partner of that name is
   absent from the atom's conformer;
3. the nearest heavy atom of its own residue within bonding distance (sum of
   covalent radii + 0.4 A, as capping's guess mode) exists and is
   altloc-compatible; this is the parent;
4. the parent has a bonded neighbour other than these H, and a dangling
   valence that accommodates it: the parent misses
   at least one heavy dictionary partner, and the number of such H on the
   parent does not exceed the parent's n_fill (rule 2 of `capping.py`, with
   the candidates themselves left out of `used`).

This is the test "if that bond were dangling, would capping put H there, and
how many". Counts and chemistry are compared, not coordinates: caps move in
refinement and the free rotation of a methyl cap depends on its environment.

Consequences:

- HB2 and HB3 on a truncated CB are dictionary atoms attached to a present
  atom. They are not caps, stay through deposition, and are correct there.
- An H with a wrong name on a complete residue fails condition 4 (its
  parent misses no heavy partner) and remains an unknown atom: the program
  stops as today.
- An H with a wrong name on a truncated residue, in excess of what the
  valence allows, fails condition 4 (count above n_fill) and stops as
  today. One within the allowance is a cap, which is what it chemically is.
- A residue without any dictionary has no caps by this definition. Such a
  model cannot be processed anyway; once the user supplies its CIF the
  definition applies with that CIF's names. For `capping_guess_for_unknowns`
  output the reduce2 description file stays the record.

## 4. Recognition: `capping.find_capping_h`

One function in `mmtbx/hydrogens/capping.py` implements section 3 and is the
only implementation of it:

```
find_capping_h(hierarchy, mon_lib_srv, bond_table=None, log=None)
  -> group_args(selection,  # flex.bool over hierarchy.atoms()
                records,    # per cap: i_seq, parent_i_seq, coordination of
                            # the parent after filling, steric number
                rejected)   # (atom label, reason) for candidates failing condition 3 or 4
```

- Candidates (condition 2) come from the hierarchy and the dictionary alone:
  per residue name, the dictionary atom names and the H-to-partner map,
  cached per name (`_dictionary_partners` already does this with the CCD
  first, then the restraint dictionary). One set lookup per H atom.
- Conditions 3 and 4 are tested only for residues that have candidates. Bonded
  neighbours of a parent come from `bond_table` when given (interpretation
  mode: exact, inter-residue links included). Without it (selection keyword
  on an unprocessed hierarchy) they are found by distance within the residue
  group and its chain neighbours, where an atom of the parent's own residue
  counts only if its name is a dictionary partner of the parent (atoms of
  the neighbouring residues count by distance alone): pdb_interpretation
  bonds nothing else inside a residue, so the keyword path stays at least
  as permissive as the bond-table path (ruling 2026-11-04, the NAD_594_HD
  case of a ligand whose name collides with an unrelated CCD entry); a link
  to a remote residue (disulfide, covalent ligand) can then be missed, which
  raises n_fill by one and makes acceptance more permissive, never stricter,
  with one exception: a parent whose only bonded neighbour is such a link
  has no neighbour in this mode, and its H are rejected (no bonded
  neighbour) where the bond-table path accepts them. Documented in the
  docstring.
- `find_dangling_atoms` is refactored so that its valence budget for one
  atom (V, present partners, missing names, n_fill) is a helper that
  `find_capping_h` calls; no second copy of the rule.

Cost: one pass over H atoms with set lookups, plus work proportional to the
number of candidate residues (a handful in practice). No pair search over the
model. The mon_lib server is the caller's; the keyword path uses the shared
server of the process (see section 7).

## 5. pdb_interpretation

At the rescue hook (after proxies, before the unknown-type stop),
regardless of `stop_for_unknowns`, when the pdb_interpretation parameter
`recognize_capping_h` is True (the default):

1. Call `find_capping_h` with the bond table built from the bond proxies
   just made.
2. For every cap, add:
   - a bond proxy cap-parent, ideal `capping.xh_length(element,
     coordination, use_neutron_distances)` from the interpretation
     parameters, sigma 0.02 A;
   - angle proxies cap-parent-X for every bonded neighbour X of the parent
     (other caps included; neighbours in other residues included), ideal by
     the parent's steric number as in `cap_dangling_atoms` (109.5, 120 or
     180 degrees), sigma 3 degrees;
   - the nonbonded energy type of a generic H (the energy library's `H`);
   - origin id `capping` on the bond and angle proxies (section 6).
   No dihedral, planarity or chirality proxies: caps are H on a dangling
   atom and the QM energy, not these restraints, decides their geometry.
3. Log one line per cap in the geometry restraints log ("Capping H: N
   recognized") and one per rejected candidate with its reason. Rejected
   candidates stay unknown atoms and reach the existing stop.
4. Keep the cap i_seqs on `all_chain_proxies` for the log and for
   `geometry_restraints_manager` bookkeeping; nothing else stores the
   answer, every consumer asks the keyword or the function.

Riding H: the riding manager parameterizes a cap from its bond and angle
proxies like any H with one heavy neighbour; where it refuses (as it already
does for dictionary H whose neighbour is missing) the cap refines as a free
atom, which is today's behaviour for such H. Verified by test, not changed.

- Limitation: the hook runs before `construct_geometry_restraints_manager`
  adds disulfides, automatic links (covalent ligands, glycans, metals), bond
  edits and symmetry bonds. A parent's partner bonded by such a later link is
  not seen: the cap gets no angle to it, and the parent's valence count
  misses it, so n_fill is over-estimated by one (permissive: an extra cap
  can be accepted, a correct one is never rejected).

- `recognize_capping_h = False` in reduce2's own H placement
  (`reduce_hydrogen.get_reduce_pdb_interpretation_params`). Placement adds
  every dictionary H, also those whose parent heavy atom is absent, and
  after `process()` removes the ones `isolated_atoms_selection()` finds
  ("when heavy atom is missing, H needs not to be placed"). Such a
  placeholder sits near the parent's neighbour and meets the cap definition
  (the backbone H of a SER without N on its CA, the H8 of a DG without C8 on
  N9), so with recognition on it got a bond, was no longer isolated and
  survived with occupancy 0 (tst_add_hydrogen_12 tst_04). During placement
  every H is provisional, so nothing there is a cap. reduce2 sets it back to
  True in `_ReinterpretModel` (so `_RebuildModel`: the qmready copy for
  regularization, the clash check and capping), which processes placed H
  only, and a finished model may carry caps. With it False a capped model's
  caps are unknown atoms again (Sorry under the default stop_for_unknowns).

## 6. Origin id `capping`

Add to `auto_linking_types.py`:

```
['capping',
 'Capping H placed by mmtbx.reduce2 capping=True on atoms with a dangling '
 'valence (truncated side chains, incomplete ligands); for QM use',
 '',
 ['Capping H', 'Capping H'],
 [0, 1]]
```

Effects, all by existing code: geometry statistics and validation report
covalent geometry only, so bond and angle rmsd are those of the model without
caps; the .geo file lists cap restraints under their own header; the charge
calculation, riding connectivity and capping's bond table still see the
bonds.

## 7. Selection keyword `capping_h`

Add the keyword `capping_h` to `iotbx.pdb.atom_selection.cache`: a slot,
`sel_capping_h`, and a parser branch beside `single_atom_residue`. The slot
is computed lazily by `mmtbx.hydrogens.capping.find_capping_h` with a
process-wide mon_lib server cached at module level in `capping.py`
(`shared_mon_lib_srv()`, the pattern of qrefine's `get_mon_lib_srv`; cctbx
has no such accessor today), imported inside the method as `hierarchy.py`
does. Without mmtbx on the path the keyword raises a Sorry naming it.

Then every program gets it for free:

```
phenix.pdbtools model.pdb remove=capping_h
model.selection("capping_h")
hierarchy.atom_selection_cache().selection("not capping_h")
```

The name follows the multi-word keywords (`single_atom_residue`,
`small_molecule`). Alternatives considered: `cappingh` (harder to read),
`capping` (reads as an operation). Decision D1.

## 8. cap_names

`cap_names` avoids every atom name of the residue's dictionary except the
parent's own free H names, instead of only the names present in the
conformer. New files never carry a cap under a dictionary name whose partner
is absent (the HO3' case) or whose partner is present but whose H is absent
(HG on a SER that has OG but lost its HG). Files written before this change
are still read correctly in the first case by condition 2 of the definition; in
the second case the dictionary attaches the cap to the present partner and
it shows up as a bond outlier of several Angstrom, not as a cap. Accepted as
a limitation of old files (D6).

## 9. Occupancy of caps: zero

Caps are not part of the deposited model and must not contribute to Fcalc.
`cap_dangling_atoms` (`_new_h`, formerly `set_occ(parent.occ)`) gives
occupancy 0 to every H it places whose name is not a dictionary atom name of
the residue (HCG, HC1: the caps of section 3). The dictionary-named H it
places (HB2, HB3 on a truncated CB) are ordinary dictionary H, not caps by
section 3, and take the parent's occupancy as any reduce2 H does (controller
ruling 2026-11-04). B, segid and hetero flag stay the parent's. Occupancy
never enters restraints or QM targets (restraints are occupancy-agnostic in
this code base), so a cap at occupancy 0 is a full atom for geometry and QM
and absent from scattering.

Code that would undo it, checked 2026-11-04:

- `mmtbx.model.manager.reset_occupancy_for_hydrogens_simple` (called by
  reduce2's placement) and `reset_occupancies_for_hydrogens` (called by
  phenix automation) copy the parent's occupancy to every bonded H. Both
  skip capping H (`self.selection("capping_h")`), so a cap stays 0 through
  a re-run of reduce2 on a capped model or an automation pipeline.
- `mmtbx.refinement.occupancies.occupancy_selections` already leaves
  occupancy-0 atoms out of individual occupancy refinement
  (`(occupancies != 1.) & (occupancies != 0.)`), but a cap inside an
  alternate conformation would be swept into its conformer's constrained
  group and pulled to the group's value. The function removes capping H
  from every selection it returns (one `remove_selections` call).
- phenix.refine's own macro cycle does not reset H occupancies (grep of
  `phenix/refinement/*.py`); qr.refine has no occupancy code at all.
- Code that derives a conformer's occupancy must take it from non-H atoms,
  which `reset_occupancy_for_hydrogens_simple` already does; nothing else
  found does otherwise.

Occupancy 0 is a property reduce2 sets, not part of the definition in
section 3: a user who edits occupancies does not change what is a cap.

## 10. qr.refine and phenix.refine

Nothing changes in either program. qr.refine loads the capped model through
`mmtbx.model.manager.process(make_restraints=True)` as today; the total
charge (`charges.py`, `electrons.run` on the hierarchy) counts the caps as
bonded H; riding H and the cctbx restraints path work through the proxies;
the QM engines see the caps as atoms. phenix.refine refines caps as H with
riding or free H per its usual strategy; its geometry statistics exclude
them by origin id.

Stripping caps from a refined model is `phenix.pdbtools remove=capping_h`.
No qr.refine output parameter is added (D3): one way to do it, already
available to every program.

## 11. Tests

cctbx_project, `mmtbx/hydrogens/tst_capping_h.py` (registered in
`mmtbx/run_tests.py`), on the fixtures of `tst_capping.py` plus a
nucleotide missing O3':

1. **Add and recognize.** `cap_model` on a fixture, then
   `find_capping_h` on the capped hierarchy: the selection is exactly the
   atoms `cap_dangling_atoms` added, no more (HB2, HB3 on a truncated CB
   are not selected), no fewer. Every cap has occupancy 0; no other atom
   changed occupancy.
11. **Occupancy survives.** On a processed capped model with an
    alternate conformation, `reset_occupancy_for_hydrogens_simple`,
    `reset_occupancies_for_hydrogens` and `occupancy_selections` leave
    the caps at 0 and out of every occupancy group; the structure factors
    of the capped and uncapped models are identical.
2. **Process.** `mmtbx.model.manager(...).process(make_restraints=True)`
   with default `stop_for_unknowns` succeeds; every cap has one bond proxy
   and its angle proxies with origin id `capping`; the nonbonded registry
   has no unknown types.
3. **Statistics.** Bond and angle rmsd, and proxy counts, of the capped
   model equal those of the uncapped model; the .geo file has the "Capping
   H" header.
4. **Idempotent.** `find_dangling_atoms` on the processed capped model
   returns no record for the capped parents.
5. **Remove.** `model.selection("capping_h")` and
   `hierarchy.atom_selection_cache().selection("capping_h")` equal the
   selection of test 1; `phenix.pdbtools remove=capping_h` on the written
   file gives the pre-capping atom list in the same order;
   `remove="not capping_h"` keeps exactly the caps.
6. **Names.** The nucleotide missing O3' gets a cap not named HO3'; a
   hand-written file with the old HO3' cap is recognized as a cap and
   bonded to C3'.
7. **Still stops.** A wrongly named H on a complete residue, and one too
   many on a truncated residue, keep the unknown-type Sorry.
8. **Riding and minimization.** `mmtbx.geometry_minimization` with riding
   H on the capped model runs; caps end bonded to their parents at the
   X-H length.
9. **Neutron.** With `use_neutron_distances=True` the bond ideal is the
   neutron X-H length; a D cap is recognized.
10. **Speed.** `find_capping_h` on a 20 000-atom model without caps costs
    less than the dictionary classification already paid (timed, asserted
    loosely like `tst_reduce2_capping_timing.py`).

qrefine, `tests/unit/tst_53.py` (own folder, as every qrefine test):
`mmtbx.reduce2 qmready=true capping=true` on a small protein with truncated
side chains, then `qr.refine` on the output with `restraints=cctbx` and with
the QM engine the existing tests use (`mopac`), a few cycles each; the
refined output keeps the caps; `phenix.pdbtools remove=capping_h` on it
returns the atom set of the reduce2 input plus the H reduce2 added. The
total charge qr.refine reports equals reduce2's Capping report.

phenix_regression (location per the plan, with the other phenix.refine
tests): `phenix.refine` one macro cycle on the same capped model, default
strategy; runs, keeps the caps at occupancy 0, reports bond and angle rmsd
equal to the uncapped model's at cycle 0, R factors at cycle 0 equal to the
uncapped model's, and `remove=capping_h` strips them.

phenix_regression, a second test (user amendment 2026-11-04): on the same
capped model, `phenix.geometry_minimization`, `phenix.model_statistics`,
`phenix.real_space_refine`, `phenix.dynamics` and `phenix.maps` each run;
the minimizers and dynamics keep the caps at occupancy 0 and bonded to their
parents; model_statistics reports the bond and angle rmsd of the stripped
model; phenix.maps writes map coefficients identical to the stripped
model's.

## 12. Implementation outline

cctbx_project:

1. `auto_linking_types.py`: the `capping` origin id.
2. `capping.py`: the valence-budget helper split out of
   `find_dangling_atoms`; `find_capping_h`; `cap_names` change; caps at
   occupancy 0; module docstring updated (Limitations: capped models can
   now be processed; old-file HG case; occupancy).
3. `pdb_interpretation.py`: the rescue step of section 5, before the pH and
   enol hooks.
4. `iotbx/pdb/atom_selection.py`: the `capping_h` keyword.
5. `mmtbx/model/model.py` (two occupancy resets) and
   `mmtbx/refinement/occupancies.py`: skip capping H.
6. Tests of section 11, registered.

qrefine: `tests/unit/tst_53.py`, `run_tests.py` registration, CHANGES.md
line. phenix_regression (`phenix_regression/refinement`): the phenix.refine
test.

Order: 1, 2 with tests 1, 4, 6, 10; then 3 with tests 2, 3, 7, 8, 9; then 4
with test 5; then 5 with test 11; then the qrefine and phenix.refine tests.

## 13. Decisions to confirm

- D1. Keyword name `capping_h`.
- D2. Cap angle ideals by steric number with sigma 3 degrees, bond sigma
  0.02 A, no dihedral or planarity proxies.
- D3. No qr.refine output parameter; `phenix.pdbtools remove=capping_h` is
  the one way to strip caps.
- D4. Candidates that fail the valence test stay unknown atoms and stop the
  run as today (no new leniency).
- D5. Residues without a dictionary have no caps by definition; the reduce2
  description file is the record for guess-mode output.
- D6. Old files with a cap under a dictionary H name whose partner is
  present are not rescued (bond outlier instead); new files cannot produce
  them.
- D7. Caps get occupancy 0 (user amendment 2026-11-04); the two H occupancy
  resets and `occupancy_selections` skip them. Files capped before this
  change keep the parent occupancy until re-capped or edited
  (`phenix.pdbtools` can set it: `modify.selection=capping_h` with
  occupancies set to 0).
