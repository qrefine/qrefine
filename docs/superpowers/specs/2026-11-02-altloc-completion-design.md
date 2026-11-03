# Altloc completion for qmready models and qr.refine: design

Date: 2026-11-02. Status: approved 2026-11-02 (D1-D6 settled).

## 1. Goal

`mmtbx.reduce2 qmready=true` and `qr.refine` must accept the alternative
conformations (altlocs) real PDB entries carry. Today both stop on nearly
every high-resolution entry: reduce2 at its atom-complete check, qr.refine at
`conformers._check_complete`. The 2026-10-28 design (decision D3) required
every residue with altlocs to carry every label of the model and left the
general case as a planned extension. This document is that extension.

What this design does: define, per residue group and per label X, which atoms
make up conformer X, complete the model where that needs atoms, and build
qr.refine's conformer models from that definition. What it does not do: water
H orientation, altloc labels for occupancy-encoded alternates (overlapping
unlabelled atoms), and the cost of one full gradient evaluation per label.
Those are separate work.

## 2. Findings (1US0, 1PWM, 1R6J and ten other entries)

Every flag of reduce2's check on the three entries falls into one of four
situations (`rt/run_1793639952_reduce2/classify.log`):

| Situation | 1US0 | 1PWM | 1R6J | Example |
|---|---|---|---|---|
| Residue lacks a label the chain has elsewhere ("inherited") | 42 of 44 | 1 of 5 | 18 of 23 | ARG 40 has A/B; SER 282 has C on CB and OG |
| An atom exists under exactly one label ("single-source") | 0 | 4 | 0 | LYS 85 CB is A-only, CG..NZ are A/B |
| An atom exists under several labels ("multi-source") | 0 | 0 | 5 | SER 205 OG is A/B/C/D on an A/B backbone |
| Atoms absent under every label ("truncated") | 2 | 0 | 0 | HIS 312 ring, GLU 313 (N only) |

The inherited situation is the common case: eight of ten further entries
(1EJG, 1GCI, 1X6Z, 1YK4, 2OV0, 2VB1, 3NIR, 2B97) have residues with fewer
labels than the model, and in every entry with three or more labels. The
in-residue situations are rarer but recur (2VB1 has one). Two further
observations:

- In all three entries every C or D label is a side-chain sub-state of A by
  occupancy: 1R6J SER 205 OG A+C+D = 0.43+0.18+0.20 = 0.81 = backbone A,
  OG B = 0.10 against backbone B 0.19; HIS 202 ring C 0.42 + ring A 0.42 =
  backbone A 0.84; ASP 195 CG..OD2 A 0.41 + C 0.41 = backbone A 0.82.
- Bond length alone does not resolve multi-source atoms: SER 205 OG-D is
  1.42 A from CB-A and 1.31 A from CB-B. Occupancy bookkeeping does.
- 1R6J ASP 195 is genuinely truncated in conformer B: CB exists only in A
  and is 3.0 A from CA-B. Making it shared would be wrong.
- `capping=True` is not a substitute: it caps the CG of a conformer whose CB
  belongs to another label, or stops (1PWM PRO 218, 1R6J SER 205).
- `mmtbx.altloc_remediate` crashes on 1PWM and is not maintained.
- 3NIR and 2B97 carry their extra labels in other chains: reduce2 passes them
  (conformers are per chain) but qr.refine's model-wide check rejects them.

## 3. Definitions

- *Conformer labels* of a model: altloc labels that have at least one
  non-hydrogen atom (`conformers.conformer_labels`; hydrogen-only labels are
  H/D exchange sites, not conformers). Written L(model).
- *Own labels* of a residue group R: the conformer labels present in R.
  Written L(R). A residue group with L(R) empty is a plain residue.
- *Conformer X of R* (X in L(model)): the atoms of R that belong to the
  single-conformer model of label X.
- *Bonded* (atoms a, b): d(a,b) < r(a) + r(b) + 0.4 A with the covalent
  radii (C 0.76, N 0.71, O 0.66, S 1.05 A). This is the criterion
  `mmtbx.hydrogens.capping` already uses to find bonds by distance
  (`GUESS_BOND_TOLERANCE = 0.4`, `reduce_hydrogen._covalent_radius`), reused,
  not redefined. Limits: C-C 1.92, C-O 1.82, C-N 1.87 A. The test is purely
  geometric: no monomer library, so it applies to ligands and waters as it
  does to amino acids. Examples: 1R6J SER 205 OG-C is 1.41 A from CB-A
  (bonded) and 1.95 A from CB-B (not); OG-D is 1.42 and 1.31 A from them
  (both bonded, so the budget test decides); ASP 195 CB-A is 3.0 A from
  CA-B (not bonded).
- *Occupancy of label Y in R*, occ_R(Y): the smallest occupancy among R's
  non-hydrogen atoms labelled Y. It is the occupancy of Y's own atoms, the
  ones no other label borrows (SER 205: A 0.43 from OG-A, not 0.81 from the
  backbone; B 0.10; C 0.18; D 0.20). The budget test and the occupancy
  split use it for the borrowing label (the "anchor" occupancy) and for the
  source label (the "counterpart" occupancy).
- Candidates come only from atom groups with the same resname as
  conformer X's atom group (microheterogeneity: a SER conformer never
  borrows from a PRO conformer of the same residue group).

## 4. The rule: borrowing within a residue group

For every residue group R with L(R) not empty and every X in L(R):

1. Conformer X of R starts as R's blank atoms plus R's X atoms.
2. Every atom name that R has under some label but neither blank nor X is a
   *gap*. Gaps are filled iteratively: a candidate (the atom of that name
   under another label P) is accepted when it passes both tests below; an
   accepted atom joins conformer X and can anchor further candidates. The
   iteration stops when no candidate is accepted.
3. A gap that stays open means R is *truncated in conformer X* at that atom.
   This is not an error of the completion: it is reported, and the existing
   truncated-residue handling applies (reduce2: Sorry, or H caps with
   `capping=True`).

The two tests:

- *Bond test*: the candidate is bonded to at least one atom already in
  conformer X. (1PWM LYS 85 CB-A is 1.53 A from CA, shared: passes for B.
  1R6J ASP 195 CB-A is 3.0 A from CA-B: fails, B is truncated at CA.)
- *Budget test*, only when more than one candidate of the same name passes
  the bond test: the candidate must have occupancy to give. It passes when
  occ(candidate) >= occ_R(X) + occ_R(P) + (occ_R of every other label
  already borrowing this atom from P) - 0.05. (1R6J SER 205, conformer D
  borrowing CB: CB-A 0.81 >= 0.20 + 0.43 passes; CB-B 0.19 >= 0.20 + 0.10
  fails. Conformer C then: CB-A 0.81 >= 0.18 + 0.43 + 0.20 passes, exactly.) Exactly one candidate must survive; otherwise the
  case is *ambiguous* and the run stops with a Sorry naming the residue, the
  atom, the candidates and their distances and occupancies. Nothing is
  guessed.

Hydrogens are atoms like any other: an H under label A bonded to a borrowed
heavy atom is borrowed with it. Borrowing never crosses residue groups.

## 5. Model operation: `complete_altlocs`

Pure hierarchy work in a new module `iotbx/pdb/altloc_completion.py` (the
pattern of `altloc_reconciliation.py`), exposed as
`iotbx.pdb.hierarchy.root.complete_altlocs(log=None)` and wrapped by
`mmtbx.model.manager.complete_altlocs(log=None)` which also resets the X-ray
structure and the missing-atoms cache. It applies section 4 to every residue
group and edits the hierarchy under one invariant:

*Scattering neutrality.* The completed model has the same Fcalc as the input
to rounding: every edit keeps each atom site's total occupancy, and copies
sit at the same coordinates with the same B (and ANISOU). A structure factor
is linear in occupancy and the bulk-solvent mask depends on positions only,
so R-factors are unchanged before refinement starts. The refinement that
follows is the only thing that moves R. (reduce2's own H and OXT change
Fcalc far more than the completion does; that is unchanged by this design.)

- An atom accepted into *every* other own label of R is relabelled blank
  (shared). Its occupancy is left as deposited (1PWM LYS 85 CB stays 0.46):
  the shared atom is what the depositor meant, and raising it would change
  Fcalc. The log notes the occupancy was kept.
- Otherwise the accepted atom is *copied* under label X, at the same
  coordinates and B. The source atom's occupancy is divided among the labels
  that end up holding the atom (the source's own and every borrower) in
  proportion to their occ_R: 1R6J SER 205 CB-A 0.81 over A 0.43, C 0.18,
  D 0.20 (sum 0.81) gives A 0.43, C 0.18, D 0.20, with CB-B 0.19 untouched;
  the per-name total is 1.00 before and after. When the deposited
  bookkeeping is consistent (the occ_R sum to the source within 0.05) these
  are the occ_R themselves; when it is not, the proportional split still
  sums exactly to the source and no occupancy is negative, and the residue
  is reported as having inconsistent occupancies.
- Atom order: a copy is appended to the atom group of label X, which is
  created after the last existing atom group of R if R has no X atom group.
  Serial numbers are reset. i_seq order of existing atoms is unchanged.

The inherited situation needs no edit: R is complete in its own conformers.

Returns `group_args(shared, copied, truncated)`, each a list of
`(residue id_str, atom name, labels)` for the log, and raises `Sorry` for an
ambiguous case. It is idempotent: a second call changes nothing and reports
nothing.

## 6. reduce2 qmready

In `mmtbx/programs/reduce2.py` `run()`, right after `add_terminal_oxt()` and
before the atom-complete check:

```
r = self.model.complete_altlocs(log=self.logger)
```

printed as

```
Altloc completion (atoms borrowed between conformers of a residue):
  LYS A  85 CB: altloc A -> shared (needed by B)
  SER A 205 CB, CA, N, C, O: copied A -> C, A -> D
  ASP A 195 CB: A copy is 3.0 A from CA B, not borrowed; truncated in B
```

Then `missing_atoms_by_residue()` reports only truncated residues, because of
the change in section 8. Riding H are then placed per conformer on complete
conformers. The output model carries the completion, and the description
file lists it. No new PHIL parameter: completion is part of what qmready
means (as OXT and NH3 are). An ambiguous case is a Sorry before anything is
written, as today's atom-complete Sorry is.

## 7. qr.refine

### 7.1 Loading

qr.refine does not edit the model. Completion is reduce2 qmready's job, the
mandatory preparation step. At load, after `reconcile_altlocs`, `qr.py` runs
the rule of section 4 in check-only mode; if it would change anything
(an atom to share or copy) qr.refine stops with a Sorry naming the residues
and saying to run `mmtbx.reduce2 qmready=true`. A truncated conformer is not
a change and passes (its H caps come from reduce2 `capping=True`).

### 7.2 Conformer models: fallback for residues without X

`conformers.split` builds conformer X from, for every residue group R:

- L(R) empty: all atoms of R.
- X in L(R): blank atoms plus X atoms (after completion these are complete,
  or R is truncated in X, which is then true in the model too).
- X not in L(R): blank atoms plus the atoms of the *parent* label, the
  first label of sorted L(R) (A in practice).

So conformer C of 1US0 is conformer A of every residue plus the C side chain
of SER 282, and `_check_complete` goes away. The parent choice could only
matter where R's own link atom (N or C) is labelled and a neighbour's X
conformer has its own labelled link atom: a wrong parent would then break
the peptide bond in conformer X. In the thirteen entries examined
(`rt/run_1793639952_reduce2/parent_choice.py`) this never occurs: every
residue lacking a label has shared link atoms, or the neighbour's extra
label is not on its backbone. So no geometric choice is made. A guard
covers the hypothetical case: if the parent's labelled link atom is not
bonded (section 3) to the neighbour's labelled X link atom, qr.refine stops
with a Sorry naming both residues. Not a guess, not a warning.

### 7.3 Multiplicity and gradients

`multiplicity(ph)` becomes the number of conformer models each atom was
selected into (computed from the selections, not from labels): ARG 40's A
atoms in 1US0 are in A and C (m = 2), its B atoms in B only (m = 1), blank
atoms in all three. `from_conformers.target_and_gradients` is unchanged:
g_i = (1/m_i) sum over the conformer models containing atom i. Occupancies
are still never used in gradients.

### 7.4 What stays

`reconcile_altlocs` runs first and still owns the pairing of labels between
neighbours. Scratch folders, the fragment manager proxy, the forced lbfgs
and gradient_only, `altlocs.handling=native`: unchanged.

## 8. Change to `missing_atoms_by_residue`

`_missing_atoms_records` (mmtbx/model/model.py) drops a monomer mapping whose
conformer label is not among the residue group's own labels when the residue
group has own labels (1US0 ARG 40 in chain conformer C), and reports a plain
residue once, under its blank key, instead of once per chain conformer
(today GLU 313 of 1US0 appears as 'AGLU', 'BGLU' and 'CGLU'). Own conformers
are reported exactly as today, so the existing
`tst_model_atom_completeness.exercise_alt_confs` (TYR OH in A only: conformer
B incomplete) keeps its expectation. The key change for plain residues in
chains with altlocs is an interface change to flag.

## 9. Tests

All qrefine tests follow `tests/unit/tst_00.py`: models as PDB strings in
the test file, a `run(prefix=...)` that creates and enters its own folder
(`qrefine_tst_XX`), `run_fmodel` on the reference model, `run_qrefine_in` on
a perturbed copy (coordinates shaken in process, written to a file), plain
asserts, registered in `run_tests.py`. `restraints=cctbx clustering=False`,
reciprocal space: cctbx on the conformer path is an exact stand-in for QM
and keeps each test under a minute. Each scenario is one file with a
docstring that names the PDB entry it comes from, what the model looks like,
what the change does to it and what is asserted. Every test also asserts
R-work of the refined model below a small value, so each is an actual
refinement.

Made-up data: the seven-residue peptide with H of `tst_00.pdb_str_good`,
turned into each scenario in process by a shared helper module
(`tests/unit/altloc_scenarios.py`): side-chain conformers are chi1
rotations of the original, backbone conformers are rotations of one
residue's backbone about the C(i-1)...N(i+1) axis (peptide bonds to the
neighbours keep their lengths), occupancies follow the originating entry.
Each test writes the model it built, so it can be inspected.

- `tst_49` inherited conformer (1US0 ARG 40 / SER 282 pattern). Residues 1
  and 3 have A/B side chains, residue 2 has A/B/C where C is a third
  rotamer. Asserts: completion changes nothing; qr.refine runs; log says
  residues without label C take their first own label; output keeps all
  labels and occupancies; R-work drops. In process, with the C atoms placed
  exactly on A: the combined cctbx gradient of every labelled atom equals
  that of the A/B-only model (conformer C is a repeat of conformer A), and
  blank atoms satisfy 3 g_ABC = 2 g_AB + g_A with g_A the gradient of the
  single-conformer A model.
- `tst_50` single-source atoms (1PWM LYS 85 and ASN 129 patterns). Residue 2
  has CB under A only with CG.. under A/B; residue 3 has its backbone under A
  only with the side chain under A/B. Asserts: the log lists both as shared;
  in the output CB and the backbone are blank with their deposited
  occupancies; structure factors of the completed and the input model agree
  (scattering neutrality); the completed model is atom for atom (label,
  occupancy, coordinates) the model written with those atoms blank by hand;
  `mmtbx.reduce2 qmready=true` on the input writes that completed model;
  refinement runs and R-work drops.
- `tst_51` multi-source atoms (1R6J SER 205 pattern). Residue 2 is a SER
  with backbone and CB under A (0.81) and B (0.19) and OG under A 0.43,
  B 0.10, C 0.18, D 0.20, OG-D placed 1.42 A from CB-A and 1.31 A from CB-B.
  Asserts: the log says CB, CA, N, C, O copied A -> C and A -> D; output has
  CB under A/B/C/D with occupancies 0.43/0.19/0.18/0.20; structure factors
  of the completed and the input model agree (scattering neutrality);
  `mmtbx.reduce2 qmready=true` on the input writes that completed model;
  `conformers.split` on it gives four complete conformer models; refinement
  runs and R-work drops. A second model in the same file with the D side
  chain at 0.45 (budget of A exceeded, B too small, no candidate survives)
  must give the Sorry naming the residue, CB, and the A and B candidates.
- `tst_52` truncated in one conformer (1R6J ASP 195 pattern). Residue 2 has
  backbone A/B, CB under A only 3.0 A from CA-B, CG..OD2 under A/C. Asserts:
  CB is not borrowed into B and the log says so with the distance, so B
  stays truncated at CA and CG..OD2 are not borrowed into B either (no
  anchor); conformer C borrows CB-A (anchored by CG-C) and then the A
  backbone (CA-B is not bonded to CB), each a single candidate, copied with
  occupancies 0.41/0.41; the output has conformer B truncated at CA (cctbx
  restraints tolerate it, a QM run would need reduce2 capping=True);
  refinement runs and R-work drops.
- the check-only mode of 7.1 and the guard of 7.2 are unit tests in
  `tst_31` (the conformers test): a model needing completion gives the Sorry
  that names reduce2; a two-residue model whose first-label parent is not
  bonded to the neighbour's labelled C link atom gives the guard's Sorry.

cctbx (mmtbx) tests, in the same change:

- `iotbx/pdb/tst_altloc_completion.py`: section 4 on hierarchies built in
  the test, one function per situation of section 2 plus idempotence, the
  occupancy bookkeeping, H borrowed with a heavy atom, a water or ligand with
  partial labels, the Sorry text, and scattering neutrality: structure
  factors of the completed and the input hierarchy agree to 1e-6.
- `mmtbx/regression/model/tst_model_atom_completeness.py`: the inherited
  conformer is not reported; a plain residue in a chain with altlocs is
  reported once under its blank key.
- `mmtbx/regression/tst_reduce2_capping.py` (or a sibling): `qmready=True`
  on the tst_50 and tst_51 models places all H and writes the completed
  model; the tst_52 model needs `capping=True`; the ambiguous model gives the
  Sorry.

Existing tests: `tst_31.exercise_incomplete_conformers` (expects the old
Sorry) is replaced by the fallback and guard checks above; the other altloc tests
(`tst_31`, `tst_40`-`tst_44`) must pass unchanged.

## 10. Implementation outline

1. `iotbx/pdb/altloc_completion.py` and its test (section 4, 5).
2. `mmtbx.model.manager.complete_altlocs`, the `_missing_atoms_records`
   change and its test (section 8).
3. reduce2 qmready integration and its test (section 6).
4. `conformers.split` fallback, membership multiplicity, `_check_complete`
   removal, `qr.py` check-only load; `tst_49`-`tst_52`; `run_tests.py`.
5. The three entries end to end: `mmtbx.reduce2 qmready=true` (with
   `capping=true` for 1US0 and 1R6J) followed by `qr.refine restraints=cctbx`
   on each, as the acceptance run, recorded in the plan but not as a
   regression test (the water clash gate still stops reduce2 until that work
   is done; the acceptance run uses the clash listing patch of
   `rt/run_1793639952_reduce2/reduce2_list_clashes.py` until then).

Each step ends with the full suites green. Changes stay uncommitted.

## 11. Decisions to confirm

- D1 (settled) the bond test is the capping criterion (covalent radii sum + 0.4 A,
  section 3). The occupancy tolerance is the slack in the budget inequality
  of section 4, occ(candidate) >= occ(anchor) + occ(counterpart) - 0.05: it
  absorbs the two-decimal rounding of deposited occupancies (sums off by
  0.01 or 0.02 are common) without accepting a real shortfall (1R6J SER 205
  CB-B is short by 0.11). Both are constants, not parameters.
- D2 (settled) completion is always on under `qmready=True`, with no PHIL switch.
  qr.refine never edits the model; it checks and points to reduce2 (7.1).
- D3 (settled) an ambiguous borrow is a Sorry (manual fix), never a guess.
- D4 (settled 2026-11-02: physical copies, scattering-neutral) completion
  edits the model: shared atoms keep their deposited occupancy, copies split
  the source occupancy in proportion to anchor occupancies, so Fcalc and R
  are unchanged before refinement. The virtual alternative (borrow only
  inside H placement and split, write the input model as is) was rejected:
  its written model is only usable by qr.refine.
- D5 (settled) the parent for a residue without X is its first own label. No
  geometric choice: no entry among the thirteen examined needs one (7.2).
  The guard against a broken peptide bond in conformer X is a Sorry.
- D6 (settled) `tst_53` is dropped. The guard of 7.2 gets a unit test in the
  conformers test (a two-residue model built to trigger it), not a
  refinement test, because no real entry exercises it.
