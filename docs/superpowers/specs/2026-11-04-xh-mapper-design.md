# X-H mapper: nuclear coordinates in qr.refine, X-ray positions for scattering

Date: 2026-11-04. Status: approved in conversation 2026-11-04 (umbrella
decisions D5, D6). Sub-project 1 of `2026-11-04-qrefine-in-phenix-refine-
design.md`. Repositories: cctbx_project (mmtbx), qrefine.

## 1. Goal

QM restraints are natural for nuclear geometry: the X-H bond lengths they
favour are the nuclear ones (about 1.09 A for C-H). X-ray scattering comes
from electrons, and the centroid of the H electron cloud sits closer to the
parent (about 0.97 A for C-H). phenix.refine's convention is: X-ray
scattering tables carry X-ray X-H lengths in the model, neutron and electron
tables carry nuclear ones. qr.refine today refines nuclear positions and
computes Fcalc from them, so against X-ray data every H scatters from a
position about 0.1 A too far out, and the output model breaks the
convention.

After this sub-project, qr.refine refines nuclear coordinates and the
scattering side (fmodel: Fcalc, R factors, maps, data target and gradients)
always sees X-ray H positions when the scattering table is an X-ray one.
Output follows the convention. For neutron and electron tables nothing
changes. The component lives in mmtbx so phenix.refine can expand a model
to nuclear positions before handing it to qrefine (sub-project 3).

## 2. Findings

- `mmtbx.model.manager.set_hydrogen_bond_length(use_neutron_distances)`
  (model.py ~3836) moves every H along its bond to the X-ray or nuclear
  target length. It reprocesses the model when the current
  pdb_interpretation flag differs, reads the lengths from
  `geometry.get_all_bond_proxies` (`distance_ideal`), and carries no
  gradient. It is the forward move, not a mapper.
- `mmtbx.model.xh_connectivity_table(geometry, xray_structure)` (model.py
  287) lists `[i_x, i_h, const_vect, distance_ideal, distance_model]` for
  every simple bond proxy with an H or D end. The ideal length is the one
  the restraints manager was processed with. Both X-ray and nuclear lengths
  are therefore available from two restraints managers of the same
  hierarchy. Cap bonds (origin id `capping`) are simple proxies and appear.
- The riding-H manager (`mmtbx/hydrogens/riding.py`) has gradient reduction
  but re-imposes ideal angles and torsions from the parent frame. It is a
  rigid model, not a length map; it is never used inside qrefine (D5).
- `mmtbx.reduce2 qmready=true` switches on use_neutron_distances: a model
  prepared for qr.refine carries nuclear X-H lengths. A model coming from
  phenix.refine's standard macro-cycles carries X-ray lengths (riding
  idealization). The mapper must handle both on entry.
- qr.py processes the model with `use_neutron_distances=True`,
  `restraints_library.cdl=False`, `sort_atoms=False` and builds
  `calculator.cctbx_geometry` (bond and angle rmsd of the restraints
  manager) from it. Statistics and weight acceptance are on nuclear
  geometry and stay so.
- `calculator.sites` (calculator.py 203): parameters `x` are shifts on
  `sites_cart`, taken from `fmodel.xray_structure` in `_initialize`.
  `update(x)` sets `sites_plus_x` into `fmodel.xray_structure`, updates
  f_calc, calls `target_and_gradients`, which evaluates restraints on
  `sites_plus_x` (through `hdm.shrink` and the keep selection) and the data
  term through `fmodel.target_functor()`, then combines with the weights.
  `hd_mapper` (refine.py) is the precedent: shrink/expand/average with
  explicit index bookkeeping.
- `driver.refine` round-trips state through fmodel copies: the weight loop
  deep-copies the start fmodel per trial, `monitor.update(fmodel)` reads
  `fmodel.xray_structure.sites_cart()`, trial results are kept as
  `monitor.sites_cart` and the best is restored with
  `xrs.set_sites_cart(sites_cart[index])`. `results.manager` holds `model`
  and writes it (`model_as_pdb`) in `finalize` and in the per-cycle PDBs.
  The pre-optimization block copies `monitor.model.get_sites_cart()` into
  fmodel. `mmtbx.hydrogens.fit_rotatable2(model, fmodel)` runs on the input
  objects before refinement when d_min < 1.2.
- `mmtbx.hydrogens.fit_rotatable2(model, fmodel)` (hydrogens/__init__.py
  392), called by refine.run when d_min < 1.2, first calls
  `model.set_hydrogen_bond_length` with the scattering-table convention
  (X-ray lengths for X-ray tables) and `model.idealize_h_riding()` on the
  whole model, then syncs fmodel from `model.get_sites_cart()` before and
  after the OMIT-map rotation scan. So qr.refine today runs a riding-H
  idealization on high-resolution X-ray data, contrary to D5. Its other
  caller is phenix.fit_h (default behaviour must not change).
- `calculator.sites_real_space` and `driver.opt` (map or no data) do not
  use fmodel; they are out of scope (umbrella, deferred).
- Scattering tables in qr.py: `wk1995 it1992 n_gaussian electron neutron`,
  default wk1995; phenix.refine's default is n_gaussian. X-ray tables are
  the first three. qr.py calls `setup_scattering_dictionaries(d_min=1.0)`
  regardless of the data; phenix.refine uses the data resolution, which
  matters for n_gaussian (umbrella D15).

## 3. Why nuclear coordinates are the parameters

Two maps are possible. With nuclear parameters and a data-side map, the QM
energy restrains the X-H length directly and the data gradient is
chain-ruled back; the data term's along-bond component is projected out,
which is harmless because the data say little about H bond lengths. With
X-ray parameters and a restraint-side map, the QM term sees H at a fixed
nuclear length whatever the parameter length is, so nothing restrains the
parameter X-H length and it drifts under the data term alone. The first
is used.

## 4. The map and its gradients

For every mapped H with parent P, nuclear position H, bond vector v = H - P,
length r = |v|, unit vector u = v/r and X-ray target length d:

    H_x = P + d u

The map leaves P and every non-mapped atom unchanged. With s = d/r and
Q = I - u u^T (projector onto the plane normal to the bond):

    dH_x/dH = s Q
    dH_x/dP = I - s Q

Both are symmetric 3x3 matrices. For a target T evaluated at the mapped
positions, with g_H = dT/dH_x and g_P the gradient the data term already
assigns to P as an atom:

    dT/dH  = s Q g_H
    dT/dP  = g_P + (I - s Q) g_H

Implementation note: `Q g = g - u (u . g)`, so per H the cost is a few
vector operations. `expand` is the same map with the nuclear target length
d_n in place of d. Both maps set the length along the current direction, so
they are idempotent and independent of the incoming length: entering with
nuclear or with X-ray lengths gives the same nuclear model.

## 5. Component: `mmtbx.hydrogens.xh_mapper`

File `mmtbx/hydrogens/xh_mapper.py`, class `mapper`. Name open.

Construction: `mapper(model, scattering_table=None)`.
- `scattering_table` defaults to `model.get_scattering_table()`; it must be
  set (Sorry otherwise). For `neutron` and `electron` the mapper is the
  identity: `is_identity` is True, `shrink`, `expand` and `gradients`
  return their input unchanged (a copy for the site arrays), and nothing
  else is built. For every other table the mapper is active.
- Active construction needs the X-ray and nuclear ideal lengths of every
  X-H bond. The model's own restraints manager gives one set, according to
  its `pdb_interpretation.use_neutron_distances`; the other set comes from a
  deep copy of the model processed with the flag flipped (same hierarchy,
  same atom order; `make_restraints=True`, log to null_out). Both are read
  through `xh_connectivity_table`. An H or D that appears in exactly one
  simple bond proxy in both tables, with the same parent, is mapped; an H
  with no bond or more than one (bridging) is not. The two tables must
  agree on the set of mapped H and parents (assert).
- Stored: `i_h` (flex.size_t), `i_x` (flex.size_t), `d_xray` and
  `d_nuclear` (flex.double), `n_atoms`.

Methods (all take and return full-size `flex.vec3_double`, never modify
the input):
- `shrink(sites_cart)`: mapped H set to `P + d_xray u`.
- `expand(sites_cart)`: mapped H set to `P + d_nuclear u`.
- `gradients(sites_cart_nuclear, gradients_xray)`: the chain rule of
  section 4 applied to a gradient computed at `shrink(sites_cart_nuclear)`;
  returns the gradient with respect to the nuclear sites.
- `selection()`: flex.bool of the mapped H, for tests and diagnostics.
- A bond of zero length (r < 1e-6) is left unmapped in all three methods
  for that call; it cannot occur in a sane model and must not raise.

Cost: one extra pdb_interpretation on a model copy at construction, seconds
for a protein, once per qrefine call. The maps are O(number of H).

Not a riding model: no angles, no torsions, no rigid geometry; the H keep
the direction the refinement gives them.

## 6. Use in qrefine (reciprocal space, mode=refine)

Invariant: `fmodel.xray_structure.sites_cart()` holds the scattering view,
`mapper.shrink(sites_nuclear)`; the nuclear sites travel explicitly through
the calculator, the monitor and the result. Nothing reads nuclear
coordinates back from fmodel.

- `refine.run` builds the mapper right after the model is processed
  (identity for neutron and electron, or in mode=opt and real-space mode
  where it is not used) and logs whether it is active and how many H it
  maps. The model is expanded once on entry:
  `model.set_sites_cart(mapper.expand(...))`, then fmodel receives the
  shrunk sites and updates f_calc. `fit_rotatable2` gains two keyword
  arguments whose defaults preserve its behaviour: `reset_h=True` (False
  skips the bond-length reset and the riding idealization) and
  `mapper=None` (when given, both fmodel syncs set
  `mapper.shrink(model.get_sites_cart())`). refine.run calls it with
  `reset_h=False, mapper=mapper` after the expansion, so it only rotates
  H on the nuclear model; afterwards the calculator is re-synced with the
  model's nuclear sites. The "starting r-factors" are therefore computed
  on X-ray H positions (a change for X-ray runs; the log says so).
- `calculator.sites` gains `mapper`. `_initialize(fmodel, sites_cart)` takes
  the nuclear sites explicitly instead of reading fmodel; `reset_fmodel`
  likewise. `update(x)` sets `mapper.shrink(sites_plus_x)` into
  `fmodel.xray_structure`. `target_and_gradients`: restraints on the nuclear
  `sites_plus_x` as now (hdm.shrink, keep selection); data gradient `dg`
  from the target functor, then `dg = mapper.gradients(sites_plus_x, dg)`,
  then the existing exclude zeroing and `hdm.average`. The order matters:
  the mapper acts on the full-size gradient before any selection or
  averaging. (hd_mapper is only active for neutron data, where the mapper
  is the identity, so the two never both act; the code path is still
  general.) The calculator exposes `sites_cart_nuclear` (its current
  `sites_plus_x`).
- `driver.refine`: `monitor.update(sites_cart=calculator.sites_cart_nuclear,
  fmodel=fmodel)` replaces reading sites from fmodel; the trial bookkeeping
  stores nuclear sites; restoring the best trial sets nuclear sites into
  the calculator and shrunk sites into fmodel through one helper
  `calculator.set_sites(sites_cart_nuclear)` that does both. The
  pre-optimization block uses the same helper.
- `results.manager.update(fmodel=..., sites_cart=...)`: sites are passed
  explicitly; `self.model` stays nuclear; bond and angle rmsd stay nuclear.
  The manager takes an `output_mapper` (None or the mapper) and applies
  `shrink` to the model's sites in `write_pdb_file` and `finalize`; qr.py
  passes the mapper unless `write_nuclear_h_positions=True`.
- New qr.py parameter: `write_nuclear_h_positions = False` (bool, help:
  "Write the refined nuclear H positions. Default: H positions follow the
  scattering-table convention, X-ray X-H lengths for X-ray tables, nuclear
  for neutron and electron."). Note for sub-project 2: the API result
  carries both `sites_cart` (convention) and `sites_cart_nuclear`.
- Real-space mode (`calculator.sites_real_space`) and `driver.opt` are
  unchanged.
- qr.py (umbrella D15b, done here so the X-ray tests are re-baselined
  once): default `scattering_table` becomes `*n_gaussian`;
  `setup_scattering_dictionaries` uses `fmodel.f_obs().d_min()` when
  diffraction data are present (`maptbx.d99` of the map in real-space mode,
  1.0 A without data as now).

## 7. Tests

mmtbx, `mmtbx/hydrogens/tst_xh_mapper.py` (registered in
`mmtbx/run_tests.py`):
1. Lengths. A small model with explicit H (protein residues with
   methyl, OH and NH3+, one water with H, one cap from the tst_173 ARG
   fragment, one residue with altloc H on a blank-altloc parent, one D
   atom, one stray unbonded H), X-ray table. `shrink` gives every mapped H
   the X-ray `distance_ideal` of its bond, `expand` the nuclear one, heavy
   atoms and the stray H are untouched, directions are preserved,
   `expand(shrink(x))` and `shrink(expand(x))` are idempotent, and the
   mapped set equals all H but the stray one.
2. Gradients. Finite differences: T(nuclear) = F(shrink(nuclear)) with F a
   smooth test function of all sites (a sum of pairwise terms so parents
   and H couple), central differences with step 1e-5 on every coordinate
   of every atom, compared with `mapper.gradients` to 1e-6 absolute. Then
   the same with F the fmodel least-squares target of a small P1 model, to
   check the integration with `gradients_wrt_atomic_parameters`.
3. Identity. The same model with `neutron` and `electron` tables:
   `is_identity`, all methods return copies equal to their input.
4. Entry from either convention. The model with X-ray lengths and the same
   model with nuclear lengths both expand to the same nuclear model.

qrefine (next free numbers at writing time: tst_54, tst_55):
5. `tst_54`: finite-difference test of `calculator.sites.target_and_gradients`
   with the mapper active, restraints=cctbx, expansion=False, clustering
   off, X-ray fmodel from `run_fmodel` on a small model with H: the total
   gradient equals the central difference of the total target over all
   coordinates of a sample of H and parent atoms to 1e-5 relative. The same
   calculator with scattering_table=neutron reproduces the current code
   path (no mapper) to machine precision.
6. `tst_55`: end-to-end `qr.refine` with restraints=cctbx against X-ray
   data: in the output file every X-H length equals the X-ray target of its
   bond (within 1e-3 A); with `write_nuclear_h_positions=True` every length
   equals the refined nuclear one, longer than the X-ray target; the final
   R factors in the log equal those of an fmodel rebuilt from the output
   file within 6e-4 (hence the convention output is what scattered; 1e-4
   is out of reach because the refine cycles never recompute f_mask, a
   pre-existing behaviour, and the rebuilt model carries PDB rounding; a
   negative control shows the nuclear output fails this check), and the
   convention output sites equal `mapper.shrink` of the nuclear output
   sites within 2e-3 A; with
   `scattering_table=neutron` against the same data the output lengths are
   nuclear and the log reports no mapper.
   tst_55 runs with the new default table and, once, with
   `scattering_table=wk1995` to show the mapper result does not depend on
   the table.
7. Existing suite: every qrefine test with neutron or electron data, or
   without data, gives unchanged numbers. (Found while implementing: a grep
   for neutron, exchangeable and hd_mapper matches tst_04, 10, 11, 25, 26,
   27, 29, 46, but tst_04, 10, 11, 26, 27, 29 and 46 use no data and set
   use_neutron_distances only for restraint geometry, and tst_25 is an
   X-ray run with exchangeable H/D, so hd_mapper and the X-H mapper act
   together there; it asserts no numbers. Only tst_55's own neutron case
   uses neutron data. All passed unchanged. Only tst_00 and tst_01 moved,
   both because their synthetic data came from H at nuclear lengths; their
   data are now generated from H at X-ray lengths and the original bound
   stands. The n_gaussian default moved no asserted number.) X-ray tests with
   explicit H whose asserted numbers move (start R factors, refined
   coordinates, scattering table) are re-baselined deliberately, each with
   a note of which of the two changes (mapper, table) moved it.

## 8. Implementation outline

1. `mmtbx/hydrogens/xh_mapper.py` and `tst_xh_mapper.py`, registered.
2. `calculator.sites`: mapper argument, explicit nuclear sites, `set_sites`
   helper, `sites_cart_nuclear`, gradient mapping.
3. `results.manager`: explicit sites in `update`, `write_pdb_file`,
   `finalize`.
4. `driver.refine`: use the helpers; no site read-back from fmodel.
5. `refine.run`: build the mapper, expand on entry, sync around
   fit_rotatable2, pass to the calculator, choose output sites;
   `write_nuclear_h_positions` in qr.py's PHIL; default table n_gaussian
   and dictionaries at the data resolution in qr.py.
6. qrefine tst_54, tst_55; run the whole unit suite and re-baseline.

## 9. Decisions confirmed

- Nuclear parameters, data-side chain rule (D6).
- Output follows the scattering-table convention; `write_nuclear_h_positions`
  forces nuclear (D5).
- No riding H inside qrefine, ever (D5): the riding idealization that
  `fit_rotatable2` performed for qr.refine at d_min < 1.2 is switched off
  through its new `reset_h` argument (found 2026-11-04 while planning).
- Statistics, weight acceptance, clustering and expansion see nuclear
  geometry only.
- Real-space mode unchanged in this sub-project.
- qr.py default scattering table n_gaussian, dictionaries at the data
  resolution (umbrella D15b, approved 2026-11-04).
