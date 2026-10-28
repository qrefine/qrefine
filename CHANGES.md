# Changes

All notable changes to this project are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Changed
   - updating documentation.
   - qr.finalise runs mmtbx.reduce2 qmready=True (output modelH_charge<q>.pdb).
   - QM fragment capping uses mmtbx.hydrogens.capping (qrefine/fragment_capping.py,
     shared with mmtbx.reduce2 capping=True): one H per cut bond (NH2 on N, one H
     on C, HG on a cut SG), truncated side chains capped, charge from the capped
     fragment; caps sit inside their residues, the cut carbonyl cap is named HOXT.
### Removed
   - finalise.py, completion.py and the use_reduce parameter.

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
