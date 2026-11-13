# Changes

All notable changes to this project are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
   - Alternative conformations (altlocs) are supported by every route
     (plain, expansion, clustering) and every restraints source (cctbx and all
     QM engines). New scope altlocs { handling = *conformers native ;
     reconcile = True }. handling=conformers (default) splits the model into
     one single-conformer model per conformer label (an altloc label with a
     non-hydrogen atom; the model holds the blank atoms plus the atoms of
     that label; hydrogen-only labels such as H/D exchange sites belong to
     every conformer), builds one ordinary restraints manager per conformer
     (scratch folder ase/<label>/) and combines the gradients (per atom, the
     mean over the conformers containing it; occupancies are not used).
     handling=native hands the model with altlocs to cctbx as is.
   - altlocs.reconcile=True relabels altlocs of interacting residues
     (iotbx.pdb.hierarchy.reconcile_altlocs) before restraints are built, so
     that the same label means the same conformer; the output model carries
     the reconciled labels. Unresolved groups are reported and refinement
     continues.
   - On the conformer path the energy is undefined, so minimizer=lbfgs and
     gradient_only=True are forced.
   - altlocs.handling=native with restraints=qm stops with a Sorry.
   - A residue group lacking a conformer's label takes the atoms of its first
     own label for that conformer; a Sorry stops the run only when this
     breaks a peptide bond (conformers.check_links). Exception: a residue
     group with no blank atoms, bonded to no other residue group, exists
     only in its own conformers (see Changed, clustering with altlocs).
   - tst_45: finite-difference check of the gradients refinement uses with
     altlocs (cctbx, both handlings; MOPAC, conformers). The sum of the
     conformer energies is the test function: its finite differences must
     equal multiplicity times the combined gradient.
   - tst_46: the same for the expansion and clustering routes, whose
     gradients are partial derivatives with the symmetry copies, buffer
     atoms and caps frozen: the frozen-environment energies are
     differentiated (cctbx, both handlings; MOPAC, conformers), and the
     cctbx routes are compared with plain cctbx restraints on an explicit
     2x2x2 supercell in a box where the images clash with the model.
   - MOPAC calculator keyword EPS (COSMO dielectric, default 78.4 as before;
     None runs in vacuum). With COSMO, MOPAC's energy is not smooth at the
     0.01 kcal/mol level, so its gradients agree with finite differences only
     to a few percent; in vacuum they agree to 0.2 %.
   - tst_47: the auxiliary command lines qr.cluster, qr.charges, qr.restraint
     and qr.gtest run (cctbx restraints).
   - tst_53: qr.refine refines a model capped by mmtbx.reduce2 capping=True
     (cctbx restraints and MOPAC); the caps keep occupancy 0 and
     phenix.pdbtools remove=capping_h strips exactly the caps.
### Fixed
   - Fragments were keyed by atom serial numbers converted with int():
     serial numbers above 99999 are hybrid-36 strings, so a super sphere of
     more than 99999 atoms stopped qr.refine. Clusters and fragments are now
     atom positions in the super sphere (fragment.py).
   - Clusters larger than cluster.maxnum_residues_in_cluster that betweenness
     centrality clustering cannot split (dense blocks of residues) were cut
     in half by residue index: twice the maximum was still too large, and a
     half could be residues with no interaction between them. They are now
     split by clustering.split_cluster: the pieces of a cluster that is not
     connected, otherwise a bisection along the Fiedler vector of the
     interaction graph cutting the fewest covalent bonds (residue pairs
     bonded by distance, fragment.bonded_residue_pairs, any residue type)
     and then the fewest interactions per pair of residues separated,
     repeated until every cluster has at most the maximum (tst_48).
   - refine.hd_mapper.expand no longer asserts that every atom has a
     non-zero gradient (a legitimately zero gradient stopped the run).
   - Re-clustering during refinement (cluster.re_calculate_rmsd_tolerance)
     never happened in mode=refine: the check looked for an attribute the
     clustered restraints manager does not have, and the call it guarded
     would have failed (an unknown keyword argument). mode=refine now
     re-clusters like mode=opt once the model moved by more than
     cluster.re_calculate_rmsd_tolerance.
   - qr.cluster failed with a NameError (undefined log) and its dispatcher
     imported the wrong module; qr.gtest failed with a TypeError (missing
     argument of create_calculator); qr.restraint called functions that no
     longer exist and is now a thin wrapper of qr.refine (mode=opt,
     number_of_micro_cycles=1, max_iterations_refine=1, clustering=False).
   - Real-space refinement with data_weight=None failed with a NameError when
     printing the initial weight estimate.
   - Charge embedding: the point charges of all cut-bond atoms but the first
     were not scaled (an exhausted filter iterator), and the Turbomole charge
     file writer (write_pdb_hierarchy_xyzq_file) did not exist.
   - MOPAC plugin: the subprocess pipes are drained with communicate() (a
     chatty MOPAC could block on a full pipe); nproc is a proper keyword.
   - The backbone-connection list added to the interaction graph before
     clustering treated every residue of a chain that starts with an amino
     acid as an amino acid (waters numbered consecutively were "bonded"); the
     list is gone, bonded residues are found by the interaction analysis.
### Changed
   - updating documentation.
   - Fragmentation setup time (mmtbx.pair_interaction, outside this
     repository; results are identical): the grid search for interacting
     residue pairs keeps atoms in a cell list instead of scanning all atoms
     per grid point, evaluates the density only for pairs with a cluster
     residue (the others were never used) and only over the grid points
     around the cluster, and the wave functions are extracted from Python
     once per call; the 15 Bohr core neighbourhood is a flex selection
     instead of a Python double loop. Fragment setup on 1i07 (1942 atoms,
     15 clusters) went from 49 s to 17 s, the whole-model interaction search
     scales with the number of atoms instead of its square.
   - One monomer library server for all fragment models
     (utils.hierarchy_utils.get_mon_lib_srv; mmtbx.model.manager.set_mon_lib_srv
     and process() reusing it, outside this repository): the library was
     read again for every fragment capper.
   - The cluster/fragment integrity check (fragment.check_residue_integrity)
     asserts that every residue group is taken whole, for any residue type,
     instead of looking for CA/CB atom names.
   - The cctbx restraints manager of the super sphere is built on first use
     only (super_cell.expand.super_sphere_geometry_restraints_manager): QM
     restraints (expansion and clustering) no longer interpret the super
     sphere with pdb_interpretation for nothing.
   - charges.charges_class computes the total charge with
     mmtbx.ligands.electrons directly (no second pdb_interpretation of the
     model, no PDB serialization); the hard-wired per-residue charge
     bookkeeping, utils/electrons.py and the unused helpers of
     utils/hierarchy_utils.py are gone. qr.charges lost the debug options
     electrons and assert_correct_chain_terminii.
   - Fragments: pair_interaction.run gets hierarchy.deep_copy() instead of
     copy.deepcopy(); selection index lookups use sets.
   - dump_gradients stops the program cleanly after writing the files.
   - Dead code removed: debug_rss, the Girvan-Newman clustering (needed a
     plugin that does not exist; clustering_method keeps bcc only), the
     unused driver helpers, commented-out blocks and unused imports.
   - cctbx (outside this repository): mmtbx.validation.rama_z loads its
     database once per process (it was reloaded and converted on every call,
     4 s each, most of the run time of a small refinement), and
     libtbx.introspection queries the macOS memory size with sysctl instead
     of system_profiler (10 s on recent macOS). Both shorten every qr.refine
     run and the test suite.
   - qr.finalise runs mmtbx.reduce2 qmready=True (output modelH_charge<q>.pdb).
   - QM fragment capping uses mmtbx.hydrogens.capping (qrefine/fragment_capping.py,
     shared with mmtbx.reduce2 capping=True): one H per cut bond (NH2 on N, one H
     on C, HG on a cut SG), truncated side chains capped, charge from the capped
     fragment; caps sit inside their residues, the cut carbonyl cap is named HOXT.
   - Models with H/D exchange sites go through hd_mapper whenever exchangeable
     H/D pairs exist (also together with heavy-atom altlocs, and with A/B
     H/D labels shared with disordered side chains); the exchanged H is placed
     in the blank atom group, so it belongs to every conformer.
   - The AIMNet2 refusal of models with altlocs is removed (the conformer path
     handles them).
   - A model whose altlocs are relabeled by altlocs.reconcile is rebuilt from
     its hierarchy, so input-file LINK/SSBOND records are not used.
   - Unit tests: qr.refine runs go through one helper, tests.unit.run_qrefine_in
     (folder=None runs in the test's folder, folder=<name> in a subfolder of it
     so that several runs of one test do not share pdb/ and ase/); it returns
     absolute output paths (folder, std_out, pdb).
   - Clustering with altlocs (cluster.clustering=True, altlocs.handling=
     conformers) clusters the altloc model once and evaluates every fragment
     in its conformer variants (fragment.fragments conformer_selections;
     cluster_restraints.from_cluster): a variant is the fragment restricted
     to one conformer, relabeled to a single-conformer fragment and capped
     against that conformer's super cell; variants with the same atoms
     collapse into one job of weight w; the gradient of an atom is
     sum(w*g)/sum(w) over the jobs having it as a cluster atom, buffer atoms
     get zero. All jobs of all conformers run in one parallel pool per
     evaluation. Before, every conformer was clustered and evaluated on its
     own (from_conformers, one after the other). Job folders and engine
     labels are ase/<index> for a fragment with one variant of full weight
     and ase/<index>_<label> otherwise; the ase/<label> folders of
     conformers are used by the plain and expansion routes only, which keep
     from_conformers. A single-conformer model is unchanged: same
     gradients (bit for bit), same ase/<index> layout. The energy is the
     mean over the conformers of the sum of the fragment energies, not a
     model energy (gradient_only is forced with clustering).
     1PWM (8452 atoms, A/B): 190 fragments, 373 jobs, cctbx setup 198 s.
   - Interactions between residues that both carry altlocs were lost by
     mmtbx.pair_interaction on the altloc hierarchy (a grid point near an
     altloc residue is assigned to the A and B copies of that residue, which
     is the same residue, so no pair). The interaction list, the bonded
     residue pairs and the buffer search of the clustering route therefore
     run on single-conformer views of the model, one per conformer label
     (conformers.single_conformer_hierarchy), and are united; a fragment is
     whole residues of any view. 1PWM: 191 clusters (at that time) instead
     of the 155 of the altloc-hierarchy search.
   - Residue groups that exist in some conformers only (every route): a
     residue group with no blank atoms (atoms whose only label is on
     hydrogens count as blank), whose own labels are a proper subset of the
     model's, and which is bonded to no other residue group (cctbx
     distance_based_connectivity, hydrogens excluded, peptide links by
     iotbx.pdb.altloc_completion.bonded, bonds across symmetry not seen),
     is absent from the other conformers instead of taking the atoms of its
     first own label there; its atoms have multiplicity 1 (a water labelled
     B only). Such residue groups bonded to each other by atoms of the
     same label (a B-only two-residue ligand or glycan chain) are one
     entity, absent or falling back as a whole.
     1PWM: 10 residue groups absent from A, 19 from B (CL 1321
     among them) and no capping failure (23 of 373 variants failed before).
     A fragment whose cluster is absent from a conformer has no variant for
     it, so its weights sum to less than the number of conformers.
     Known limitation: with cluster.charge_embedding=True a fragment whose
     variants merge into one job (label "") takes its point charges from
     the smallest label's conformer sphere, so its field is conformer A's
     even where the environment within charge_cutoff differs.
   - cctbx fragment restraints are evaluated in the main process whatever
     parallel.nproc (QM jobs keep parallel_map): a process per job
     re-pickled the super-sphere restraints at every evaluation (1PWM, 373
     jobs: 21.6 s per evaluation, 0.17 s in process after a 20 s first
     call). Exceptions of cctbx fragments propagate as they are (no Sorry
     wrapper, no error folder). qr.refine on 1PWM (cctbx, clustering, one
     weight-search and one refinement cycle, nproc=4): 2344 s -> 536 s.
     R-free 0.1350 -> 0.1349 because the cached fragment restraints keep
     their updated pair proxies where the old path reselected them from the
     super sphere at every call.
   - QM capper (fragment_capping.near_selection): cut bonds are searched
     among the parent atoms in the fragment's box widened by the bond
     cutoff, not in the whole super cell: 2.35 s -> 1.4 s per capper on 1PWM
     (373 cappers, identical cut bonds on 75 compared variants). The
     default pdb_interpretation scope is parsed once per process.
   - Memory of the clustered route (1PWM): QM setup held 1394 MB (macOS
     footprint), now 391 MB (the capper drops its processed fragment model
     once the caps are placed; the relabeled conformer super cells exist only
     while the cappers are built; fragment_capped_initial is gone); cctbx
     first evaluation 5.8 s -> 1.7 s and heap after two evaluations 942 ->
     534 MB (the super-sphere restraints are built without
     plain_pairs_radius, which only ADP restraints use); the model is copied
     for the expansion route only. Results are bit-identical. The peak
     (1.6-1.8 GB footprint) is cctbx's operator search in
     super_cell.manager, unchanged.
   - Tests: tst_64-69 (variants and capper parents, weighted combination,
     duplicate-conformer identity, MOPAC super toy, qr.refine on the toy
     data, 6uf7 with real data, cctbx only), tst_70-73 registered in
     run_tests.py; tst_05, 11, 40, 42, 43, 46, 56, 59 follow the new route
     and folder layout. tst_53, 56, 57 and 59 use r_free=0.05 (qrefine's
     run_fmodel takes an r_free override) so that the generated flags pass
     cctbx's suitability test (2.5 A data with the default 0.01 give
     exactly 1 % of flags). phenix_regression tst_177 has a case with an
     altloc model and qrefine.clustering=True (phenix.refine equal to
     standalone qr.refine); tst_174, 177 and 178 were split into tests of
     1-2 minutes (tst_180-188).
### Removed
   - finalise.py, completion.py and the use_reduce parameter.
   - debug_rss, cluster.clustering_method=gnc, utils/electrons.py.
   - cluster.altloc_method, cluster.bond_with_altloc,
     restraints.from_altlocs2/from_cctbx_altlocs/from_altlocs and the
     per-altloc fragment code. Old .eff files that set the removed parameters
     now fail with an unknown-parameter error.
### Known limitations
   - The buffer of a cluster is the residues it has a density interaction
     with (one layer). A cctbx nonbonded term between residues without
     such an interaction (an N...C contact at the van der Waals distance
     two residues apart on a helix) is then missing from the fragment,
     and the clustering-route gradient of that atom is off by ~1e-4 for
     very small clusters (maxnum_residues_in_cluster of 2 on helix.pdb).
     tst_10 uses 3.
   - Altloc labels that differ only by case (e.g. a and A) share one
     ase/<label>/ scratch folder on case-insensitive file systems.
   - altlocs.handling=native with clustering (cctbx only) still runs the
     interaction search on the altloc hierarchy and can miss interactions
     between altloc residues; use handling=conformers.
   - With symmetry images in the super sphere, the one-cluster target of the
     conformers route (mean over conformers of the fragment energies) equals
     the mean of the conformer energies only without images: the altloc
     super sphere is a superset of each conformer's sphere. Gradients agree
     either way; the target is not a model energy.
   - The test that a residue group is bonded (see Changed) does not see
     bonds across symmetry: a B-only ligand bonded only to a symmetry mate
     counts as absent from A.
   - The peak memory of every clustered or expansion run is the operator
     search of cctbx's super_cell (1.6-1.8 GB on 1PWM, growing with
     atoms * (select_within_radius + 1)^3). The whole-model cctbx
     restraints that the clustering route builds and does not use (64 MB on
     1PWM) are kept for bit identity, and every job holds flex.bool
     selections over the sphere (quadratic in model size).
   - qr.refine adopts the unit cell of the data through another conversion
     path than phenix.refine: with a cell that is not exact in float32 (MTZ
     storage) the phenix hook and standalone qr.refine start 1-4 ULP apart
     and the gradient-only minimizer can amplify that to ~0.01 A with equal
     R factors (tst_177 uses a float32-exact cell).

## [1.0-0] - 2019.1.19
### 
   - Optimization  
   - Refinement
   - QM Interfaces
   - Semi-empirical Interfaces
   - Neural Network interfaces
   - Finalize models by adding missing atoms and hydrogrens
   - Fragments for divide and conquer
   - Capping of unsaturated groups
   - Point charges for modelling environment
   - Buffer region for modelling environment 
   - Symmetry is handled
   - Alternate locations are handled 
   - Works on a workstation or a cluster
   - Gradient only minimizer   
