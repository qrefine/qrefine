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
   - A residue whose altloc labels are a non-empty proper subset of the
     model's labels stops with a Sorry (its conformers would be incomplete).
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
### Changed
   - updating documentation.
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
### Removed
   - finalise.py, completion.py and the use_reduce parameter.
   - cluster.altloc_method, cluster.bond_with_altloc,
     restraints.from_altlocs2/from_cctbx_altlocs/from_altlocs and the
     per-altloc fragment code. Old .eff files that set the removed parameters
     now fail with an unknown-parameter error.
### Known limitations
   - Altloc labels that differ only by case (e.g. a and A) share one
     ase/<label>/ scratch folder on case-insensitive file systems.

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
