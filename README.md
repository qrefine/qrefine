# Quantum Refinement Module

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Code Health](https://landscape.io/github/qrefine/qrefine/master/landscape.svg?style=flat)](https://landscape.io/github/qrefine/qrefine/master)


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


### Quickstart

Please first install PHENIX, see https://www.phenix-online.org/
 
Once you have PHENIX installed, go to the directory where you installed PHENIX.

```
 source phenix_env.sh
 phenix.python modules/cctbx_project/libtbx/auto_build/bootstrap.py --builder=qrefine --use-conda --nproc=8
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
  phenix.python -m pip install ase==3.17.0
  phenix.python -m pip install pymongo
  git clone https://github.com/qrefine/qrefine.git
  libtbx.configure qrefine
```

 ### Run Tests 

``` 
 mkdir tests
 cd tests
 qr.test
```
If any of the tests fail, please raise an issue here: [issue tracker](https://github.com/qrefine/qrefine/issues)

 ### Documentation
 
 Can be found at: https://qrefine.com/qr.html
 

### Help 

If you run into any trouble please ask for help:
```
 qr.refine --help
```

### Commandline options

If you want to see the available options and default values please type:
```
 qr.refine --defaults or qr.refine --show
``` 


 

### Contact us 

The best way to get a hold of us is by sending us an email: qrefine@googlegroups.com


### Developers

* [Min Zheng](https://github.com/zhengmin317)
* [Pavel Afonine](https://github.com/pafonine)
* [Mark Waller](https://github.com/mpwaller)
* [Nigel Moriarty](https://github.com/nwmoriarty)
* [Malgorzata Biczysko](https://github.com/biczysko)
* [Xu Yanting](https://github.com/yanting0928)
* [Lum Wang](https://github.com/Mooooony)
* [Holger Kruse](https://github.com/hokru)


### Citations:

Min Zheng, Jeffrey Reimers, Mark P. Waller, and Pavel Afonine,
Q|R: Quantum-based Refinement, 
(2017) Acta Cryst. D73, 45-52.
DOI: [10.1107/S2059798316019847](http://scripts.iucr.org/cgi-bin/paper?S2059798316019847)

Min Zheng, Nigel W. Moriarty, Yanting Xu, Jeffrey Reimers,  Pavel Afonine, and Mark P. Waller,
Solving the scalability issue in quantum-based refinement: Q|R#1
(2017). Acta Cryst. D73, 1020-1028.
DOI: [10.1107/S2059798317016746](http://scripts.iucr.org/cgi-bin/paper?S2059798317016746)

Min Zheng, Malgorzata Biczysko, Yanting Xu, Nigel W. Moriarty, Holger Kruse, Alexandre Urzhumtsev, Mark P. Waller, and Pavel V. Afonine
Including Crystallographic Symmetry in Quantum-based Refinement: Q|R#2
(2020). Acta Cryst. D76, 41-50.
DOI: [10.1107/S2059798319015122](http://scripts.iucr.org/cgi-bin/paper?S2059798319015122)

Lum Wang, Holger Kruse, Oleg V. Sobolev, Nigel W. Moriarty, Mark P. Waller, Pavel V. Afonine, and Malgorzata Biczysko
Real-space quantum-based refinement for cryo-EM: Q|R#3
bioRxiv 2020.05.25.115386; 
DOI:[0.1101/2020.05.25.115386](https://www.biorxiv.org/content/10.1101/2020.05.25.115386v1)


#### Clustering

Min Zheng, Mark P. Waller, 
Yoink: An interaction‐based partitioning API,
(2018) Journal of Computational Chemistry, 39, 799–806.
DOI: [10.1002/jcc.25146](https://doi.org/10.1002/jcc.25146)

Min Zheng, Mark P. Waller, 
Toward more efficient density-based adaptive QM/MM methods, 
(2017)Int J. Quant. Chem  e25336 
DOI: [10.1002/qua.25336](https://doi.org/10.1002/qua.25336)

Min Zheng, Mark P. Waller, Adaptive QM/MM Methods,
(2016) WIREs Comput. Mol. Sci., 6, 369–385.
DOI: [10.1002/wcms.1255](https://doi.org/10.1002/wcms.1255)

Mark P. Waller, Sadhana Kumbhar, Jack Yang,
A Density‐Based Adaptive Quantum Mechanical/Molecular Mechanical Method
(2014) ChemPhysChem  15, 3218–3225. 
DOI: [10.1002/cphc.201402105](https://doi.org/10.1002/cphc.201402105 )



