# Quantum Refinement Example

## This tutorial example is still under development. 


## Installation

Before we start please install PHENIX, see https://www.phenix-online.org/

After installing, please don't forget to run:
```source phenix_env.sh```

For this tutorial we used phenix-dev-3386.

Next we need to install qrefine. 

```phenix.python modules/cctbx_project/libtbx/auto_build/bootstrap.py --builder=qrefine  ```
 
 Note: adding --nproc=N can speedup the compilation step.
 
 For this tutorial we are using qrefine version: v1.0-37-gb8d9a6  
 
 More information can be found 
 
 Please make sure you run the unit test before you continue with this example. 
 Then you know that everything is ready to go, and less likely to get a surprise. 
 
 
## Getting model and data

To obtain the PDB file and data file please visit: https://www.rcsb.org/

Alternatively, you can fetch them using a command line tool:

  
 ```phenix.fetch_pdb 1uso --mtz ```
 
 
## Pre-process the model 


0) a87_99_h.pdb is derived from 1us0.pdb by:
  - extracting 87-99 residues from chain A;
  - removing alternative conformations and setting occupancies to 1;
  - adding H atoms;
  - placing model into a P1 box: iotbx.pdb.box_around_molecule;
  - setting B-factors to 80.
  - converting all residues to GLY.
  
  This model is the reference ("answer"). It contains information about helix
  H-bonding that is defined from actual experimental diffraction data.

## 1. Compute simulated diffraction data 
  
  ```phenix.fmodel a87_99_h.pdb add_random_error_to_amplitudes_percent=5 type=real high_res=4 low_res=6 r_free=0.1```
  
  ```mv a87_99_h.pdb.mtz data.mtz```
  
  Given that we introduce 5% error into Fobs, we expect the best possible
  R-factor after refinement to be ~5%.
  
  
## 2. Perturb the model
  
  
 Generates many perturbed models - starting points for refinements.
   
   ```run_perturb.py ```
  
  
## 3. Refine the structure
  
 Run Q|R refinement, sample command:
  
  ```python ../qr-core/qrefine.py a87_99_h.pdb.mtz perturbed/1.5/4.pdb restraints=cctbx update_all_scales=False```
  
  
## 4. Analysis
  
```run_analyze.py``` is to analyze the results.

example analyze.dat 

| Pert. Dose        |Min dist.          |  Max dist  | Av. dist.         |  % recovered |
| ------------- |:-------------:| -----:|:-------------:| -----:|
| 0.3       | 1.880 | 2.734 |2.239   |   45.56 
 


