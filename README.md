# Quantum Refinement Module

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

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

Min Zheng, Nigel W. Moriarty, Yanting Xu, Jeffrey Reimers,  Pavel Afonine, and Mark P. Waller,
Solving the scalability issue in quantum-based refinement: Q|R#1
(2017). Acta Cryst. D73, 1020-1028.
DOI: [10.1107/S2059798317016746](http://scripts.iucr.org/cgi-bin/paper?S2059798317016746)

Min Zheng, Jeffrey Reimers, Mark P. Waller, and Pavel Afonine,
Q|R: Quantum-based Refinement, 
(2017) Acta Cryst. D73, 45-52.
DOI: [10.1107/S2059798316019847](http://scripts.iucr.org/cgi-bin/paper?S2059798316019847)


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
  git clone https://github.com/qrefine/qrefine.git
  libtbx.configure qrefine
```

 ### Run Tests 

``` 
 mkdir tests
 cd tests
 qr.test
```
If any of the tests fail, please raise and issue here: [issue tracker](https://github.com/qrefine/qrefine/issues)


### adding dispersion or BSSE (gCP) correction to arbitrary methods
Sometimes quantum chemistry programs do not, or insufficiently, provide the D3 dispersion or the gCP BSSE correction. This can be remedied with the following interface:

Examples:

adds "-D3(BJ)" to B3LYP calculations. Functional needs to parametrized and "-bj or -zero" needs to be specified:
```
qr.refine [usual options] quantum.qm_addon=dftd3 quantum.qm_addon_method="b3-lyp -bj"
```

adds "gCP correction for HF/6-31G" to any QM results. Basis sets needs to parametrized. Will output an error if a basis is not available (overview of parameter sets available if calling "gcp -h"):
```
qr.refine [usual options] quantum.qm_addon=gcp quantum.qm_addon_method="hf/631g"
```

The *dftd3* and *gcp* programs needed to be installed and can be obtained from Prof. Grimme [homepage](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software)
The gCP and D3 corrections are described herein:

H. Kruse,S. Grimme J. Chem. Phys. 136 , 154101 (2012) 
DOI: [10.1063/1.3700154](https://doi.org/10.1063/1.3700154)

S. Grimme, J. Antony, S. Ehrlich, H. Krieg J. Chem. Phys. 132, 154104 (2010);
DOI: [10.1063/1.3382344](https://doi.org/10.1063/1.3382344)


### Help 

If you run into any trouble please ask for help:
```
 qr.refine
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
* [Zhang Hongli](https://github.com/zhangholly)
* [Lum Wang](https://github.com/Mooooony)
* [Holger Kruse](https://github.com/hokru)




