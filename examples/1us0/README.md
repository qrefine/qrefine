0) a87_99_h.pdb is derived from 1us0.pdb by:
  - extracting 87-99 residues from chain A;
  - removing alternative conformations and setting occupancies to 1;
  - adding H atoms;
  - placing model into a P1 box: iotbx.pdb.box_around_molecule;
  - setting B-factors to 80.
  - converting all residues to GLY.
  
  This model is the reference ("answer"). It contains information about helix
  H-bonding that is defined from actual expetimenatl diffraction data.

1) Compute simulated diffraction data:
  
  phenix.fmodel a87_99_h.pdb add_random_error_to_amplitudes_percent=5 type=real high_res=4 low_res=6 r_free=0.1
  
  mv a87_99_h.pdb.mtz data.mtz
  
  Given that we introduce 5% error into Fobs, we expect the best possible
  R-factor after refinement to be ~5%.
  
2) run_perturb.py generates many perturbed models - starting points for 
   refinements.
  
3) Run Q|R refinement, sample command:
  
  python ../qr-core/qrefine.py a87_99_h.pdb.mtz perturbed/1.5/4.pdb restraints=cctbx update_all_scales=False
  
4) run_analyze.py is to analyze the results.

example analyze.dat 

| Pert. Dose        |Min dist.          |  Max dist  | Av. dist.         |  % recovered |
| ------------- |:-------------:| -----:|:-------------:| -----:|
| 0.3       | 1.880 | 2.734 |2.239   |   45.56 
 


