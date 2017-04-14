# Quantum Refinement Module

[![Build Status](https://travis-ci.org/qrefine/qrefine.svg?branch=master)](https://travis-ci.org/qrefine/qrefine)

Quantum Chemistry can improve bio-macromolecular structures, especially when only low-resolution data derived from crystallographic or cryo-electron microscopy experiments are available. Quantum-based refinement utilizes chemical restraints derived from quantum chemical methods instead of the standard parameterized library-based restraints used in experimental refinement packages. The motivation for a quantum refinement is twofold: firstly, the restraints have the potential to be more accurate, and secondly, the restraints can be more easily applied to new molecules such as drugs or novel cofactors.

However, accurately refining bio-macromolecules using a quantum chemical method is challenging due to issues related to scaling. Quantum chemistry has proven to be very useful for studying bio-macromolecules by employing a divide and conquer type approach. The fragmentation approaches we developed for achieving a quantum-refinement of bio-macromolecule will be presented.

### Citations:
Min Zheng, Jeffrey Reimers, Mark P. Waller, and Pavel Afonine, Q|R: Quantum-based Refinement , Acta Crystallographica Section D (2017) D73, 45-52. DOI: 10.1107/S2059798316019847

Min Zheng, Mark P. Waller, Toward more efficient density-based adaptive QM/MM methods, Int J. Quant. Chem (2017) e25336
DOI: qua.25336

Min Zheng, Mark P. Waller, Adaptive QM/MM Methods, WIREs Comput. Mol. Sci. (2016) DOI: 10.1002/wcms.1255


### Quickstart

Please first install PHENIX, see https://www.phenix-online.org/
 
Once you have PHENIX installed, go to the directory where you installed PHENIX.

```
 source phenix_env.sh
 cd modules 
 git clone https://github.com/qrefine/qrefine.git
 cd qrefine
 ./patch.sh
 ```
 Note: you may need to use sudo depending on the permissions of your PHENIX installation.
 
 ### Run Tests 

``` 
 qr.test
 
```
If any of the tests fail, please raise and issue here:

### Run Example 

If tests run successfully, then try and run an example: 


```
 qr.finalise --example
``` 

```
 qr.cluster --example
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


