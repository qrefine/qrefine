# Quantum Refinement Example


The structure we will use here is HIV-1 capsid C-terminal domain mutant (3dtj).

This tutorial example is still under development. 


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

  
 ```phenix.fetch_pdb 3dtj --mtz ```

Normally this would result in 3dtj.pdb and 3dtj.mtz, 
but the conversion to mtz file fails because the structure factors are
saved as 3dtj-sf.cif.

```phenix.cif_as_mtz 3dtj-sf.cif```  
 
 This results in a new file 
 
 ```3dtj-sf.mtz```
 
 Now we have the two files we need 3dtj.pdb and 3dtj-sf.mtz we can run a cctbx refinement. 
 

## CCTBX Refinement 

 To run the first refinement, please run the following command:

 ```qr.refine 3dtj.pdb 3dtj-sf.mtz restraints=cctbx```
 
 This will start a refinement on your machine starting statistics are:
 
 ``` Rw: 0.3230 Rf: 0.3113 Rf-Rw: -0.0117 rmsd(b):  0.0151 rws:  1.000 n_fev: 0 ```
 
 The refinement on our machine took around 140 seconds.
 
 At end of the refinement the best statistics are printed:
 
 ```Best r_work: 0.2526 r_free: 0.2966 ```
 
 The refined structure is saved so you can view it in VMD or PyMOl:
 
 ```3dtj_refined.pdb```
 
 
 We cannot use quantum engines on a large system, 3dtj has 2291 atoms and therefore intractable. 
 Instead we need to break the protein down in small manageable pieces. What we will do is perform 
 a separate qm calculation on each of the small manageable pieces. 
 
   

## Clustering

This command line tool will cluster the protein up into pieces 
    
 ```qr.cluster 3dtj.pdb ```


The output from the script indicates the default value of 15 being the maximum number of residues allowed
in any of the clusters. This parameter controls the size of large proteins.

```max number of residues in each cluster: 15```

The residue indices for each cluster:

```[[120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 146], [189, 190, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202], [102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113], [220, 221, 222, 223, 224, 225, 226, 260, 262, 263, 264, 265], [232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 261], [6, 7, 8, 9, 10, 11, 12, 15, 16, 19, 47], [61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71], [74, 75, 77, 94, 95, 96, 97, 98, 99, 100, 101], [86, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143], [276, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290], [20, 21, 22, 23, 24, 25, 26, 27, 28, 29], [48, 50, 51, 52, 53, 54, 55, 56, 57, 58], [73, 203, 208, 209, 210, 211, 212, 213, 214, 215], [148, 166, 167, 168, 169, 170, 171, 172, 173, 174], [153, 158, 159, 160, 161, 162, 163, 164, 165, 188], [175, 176, 177, 178, 179, 180, 181, 182, 183, 185], [249, 250, 251, 252, 253, 254, 255, 256, 257, 258], [267, 268, 269, 270, 271, 272, 273, 274, 275, 292], [78, 79, 80, 88, 89, 90, 91, 92, 93], [44, 76, 114, 115, 116, 117, 118, 119], [147, 149, 150, 151, 152, 187, 191], [1, 2, 3, 4, 5, 45], [31, 34, 35, 37, 38, 39], [227, 228, 229, 230, 231, 266], [40, 41, 42, 43, 46], [81, 82, 83, 84, 85], [243, 244, 246, 247, 248], [13, 14, 17, 18], [30, 32, 33, 36], [154, 155, 156, 157], [204, 205, 206, 207], [216, 217, 218, 219], [277, 278, 279, 280], [87, 144, 145], [59, 60], [132, 133], [49], [72], [184], [186], [245], [259], [291]]```


And finally the time taken is printed.

```Time: 53.6436```
 
 
 
## Finalize 
 
 However, another challenge with refining a protein with quantum chemistry is that we need a complete atomic model. 
 
 What that means is we need to add missing hydrogens. We also want to add missing atoms from residues.  
 
 ```qr.finalise 3dtj.pdb```
 
 This results in two new files:
 
 ```3dtj_complete.pdb ```
 
 and
 
 ```3dtj_readyset_input.pdb```
 
  The 3dtj_complete.pdb contains the complete model that we will use for the rest of this example.
  
  Just to make sure that the structure has no serious issues, we will run another refinement:
 
  ```qr.refine 3dtj_complete.pdb 3dtj-sf.mtz restraints=cctbx```
  
  
  The statistics for the refinement are  slightly different to the original refinement because we 
  added some missing atoms.
  
  ```Rw: 0.3271 Rf: 0.3041 Rf-Rw: -0.0229``` 
  
  The final refinement 
  
  ```Best r_work: 0.2832 r_free: 0.3028```
  
  The refined structure is written to your working directory:
  
  ```3dtj_complete_refined.pdb```
  
  ##Charges
  
  As another check we can use the qr.charges command line tool to get the overall charge of the model.
  
  ```qr.charges 3dtj.pdb verbose=true```
  
  resulting in:
  
  ```Charge: -1292```   
  
  Clearly, this is a nonsensical result. We will run the same script on the 3dtj_complete.pdb model.
  
  ```qr.charges 3dtj_complete.pdb verbose=true ```
  
  resulting in:
  ```Charge: 0```
  
  This shows the importance of running the qr.finalise command line tool.
  

  ## Fragmentation
  
  Another command line tool that becomes useful for debugging is qr.fragmentation.
  
  ```qr.fragmentation 3dtj_complete.pdb```
  
  This tool performs a clustering of the system, but then adds an additional layer of atoms (called the buffer)
  to the outside of each cluster. Then you can inspect each of the clusters in your viewer (e.g. VMD or PyMol).
  
  The residue indices are printed out again, but now the output also contains the 
  
  ```
  capping frag: 0
  point charge file: 0
  charge_cutoff:  8.0
  write mm pdb file: 0
  capping frag: 1
  point charge file: 1
  charge_cutoff:  8.0
  write mm pdb file: 1
  capping frag: 2
  point charge file: 2
  charge_cutoff:  8.0
  .
  .
  .
  capping frag: 36
  point charge file: 36
  charge_cutoff:  8.0 
  write mm pdb file: 36
```
  
  A bunch of new files are written to disk, they can be inspected visually to see, for example the files for the zeroth cluster: 
  
  ```0.pdb 0_frag.pdb  0_cluster.pdb ```

  The time taken is significant:

```Time: 809.7237```
  
  ## Restraint
  
  The qr.restraint command line tool is used to perform a single energy and gradient evaluation. 
  
  This tool is useful for trying different qm engines:
  
  ```qr.restraint 3dtj_complete.pdb``` 
  
  ** an issue was raised on qr.restraint that is being resolved.
  
  You can also use this tool for trying different methods and basis set combinations. 
  
  The benefit is you can get an estimate on how long a refinement will take. 
  
  This can help when making a decision on whether a refinement is feasible on your hardware.
  
  
  ## QM Refinement
  
  We use a cluster of computing nodes to perform refinement when a QM engine is needed. 
  
  The is because each of the QM calculations are expensive even when using clustering.
  
  ```qr.refine 3dtj_complete.pdb 3dtj-sf.mtz restraints=qm engine_name=terachem parallel.method=slurm``` 
  
  This will submit and run a series of TeraChem jobs on your cluster.

  

  
     
  
  

