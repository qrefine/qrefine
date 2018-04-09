# Quantum Refinement Module


Quantum Chemistry can improve bio-macromolecular structures,
especially when only low-resolution data derived from crystallographic
or cryo-electron microscopy experiments are available. Quantum-based
refinement utilizes chemical restraints derived from quantum chemical
methods instead of the standard parameterized library-based restraints
used in experimental refinement packages. The motivation for a quantum
refinement is twofold: firstly, the restraints have the potential to
be more accurate, and secondly, the restraints can be more easily
applied to new molecules such as drugs or novel cofactors.

However, accurately refining bio-macromolecules using a quantum
chemical method is challenging due to issues related to
scaling. Quantum chemistry has proven to be very useful for studying
bio-macromolecules by employing a divide and conquer type approach. We
have developed a new fragmentation approach for achieving a
quantum-refinement of bio-macromolecules.

### Citations:

Zheng, M., Moriarty, N.W., Xu, Y., Reimers, J.R., Afonine, P.V. & Waller, M.P. 
Solving the scalability issue in quantum-based refinement: Q|R#1
(2017). Acta Cryst. D73, 
DOI: 10.1107/S2059798317016746

Min Zheng, Jeffrey Reimers, Mark P. Waller, and Pavel Afonine,
Q|R: Quantum-based Refinement , 
(2017) Acta Cryst. D73, 45-52.
DOI: 10.1107/S2059798316019847


#### Clustering

Min Zheng, Mark P. Waller, Yoink: An interaction‐based partitioning API,
Journal of Computational Chemistry 2018, 39, 799–806 

Min Zheng, Mark P. Waller, Toward more efficient density-based
adaptive QM/MM methods, Int J. Quant. Chem (2017) e25336 DOI:
qua.25336

Min Zheng, Mark P. Waller, Adaptive QM/MM Methods, WIREs
Comput. Mol. Sci. (2016) DOI: 10.1002/wcms.1255

Mark P. Waller, Sadhana Kumbhar, Jack Yang,
A Density‐Based Adaptive Quantum Mechanical/Molecular Mechanical Method
ChemPhysChem 2014, 15, 3218 – 3225. 

### Quickstart

Please first install PHENIX, see https://www.phenix-online.org/
 
Once you have PHENIX installed, go to the directory where you installed PHENIX.

```
 source phenix_env.sh
 phenix.python modules/cctbx_project/libtbx/auto_build/bootstrap.py --builder=qrefine
 ```
 Note: adding --nproc=N can speedup the compilation step.

 Note: you may need to use sudo depending on the permissions of your PHENIX installation.

 ###### Using the Git repository of *cctbx*.

To remain up-to-date with the changes in the *cctbx* project that contains many
of the functions used in Q|R, remove the **cctbx_project** directory in the
modules directory. The above command will clone it from *GitHub*.

 ###### In case the quickstart command fails

 Clone the qrefine repo in the modules directory of the Phenix installation.
```
  phenix.python -m pip install ase
  phenix.python -m pip install JPype1
  phenix.python -m pip install pymongo
  libtbx.configure qrefine
```

 ### Run Tests 

``` 
 qr.test
```
If any of the tests fail, please raise and issue here: [issue tracker](https://github.com/qrefine/qr-core/issues)

### Run Example 

If tests run successfully, then try and run an example: 


```
 qr.finalise --example
``` 

```
 qr.chunk --example
``` 
 
```
 qr.restraint --example
```

```
 qr.refine --example
```

### Help 

If you run into any trouble please ask for help:
```
 qr.refine --help
```

### Contact us 

The best way to get a hold of us is via the  [issue tracker](https://github.com/qrefine/qr-core/issues)


### Developers

* [Min Zheng](https://github.com/zhengmin317)
* [Pavel Afonine](https://github.com/pafonine)
* [Mark Waller](https://github.com/mpwaller)
* [Nigel Moriarty](https://github.com/nwmoriarty)


