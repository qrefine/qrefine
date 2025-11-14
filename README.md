# Quantum Refinement Module

[![CI pipeline on Mamba](https://github.com/qrefine/qrefine/actions/workflows/ci-mamba.yaml/badge.svg)](https://github.com/qrefine/qrefine/actions/workflows/ci-mamba.yaml)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Dockerhub](https://img.shields.io/badge/dockerhub-images-important.svg?logo=Docker")](https://hub.docker.com/r/qrefine/cctbx-qr)

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

## Installation

Depending on your use case, installation of qrefine follows 3 paths:

 - [cctbx-only (open-source)](cctbx.md#cctbx-installation)
 - [phenix user](phenix.md#fully-automated-recommended)


**Requirements**:
 - python >= 3.9
 - conda binary, e.g., miniconda. (A conda environment is not needed for Phenix)
 - For Apple Silicon architecture please see the [additional notes](#apple-silicon)!


**AQuaRef notes**:
AQuaRef is fully integrated into Phenix  starting dev-5395 version.

To use AQuaRef it is recommended to install [recent Phenix version](https://phenix-online.org/download/nightly_builds.cgi).

Usage instructions are [here](https://phenix-online.org/version_docs/2.0-5867/reference/AQuaRef.html)

A few extra notes for open-souce only, also for [performance](aimnet2.md#performance), are provided here: [AQuaRef notes](aimnet2.md)


### Apple Silicon

We cannot recommend to run qrefine with *clustering* on Apple Silicon machines as the pair interaction is unreliable (unknown cause).
When following the cctbx installation route use the following to create the env:
   ```
    conda env create -n qrefine -f config/arm64-osx.yaml
   ```

For Phenix installations we recommend to switch the blas implementation to apple's accelerate in 
   ```
    conda install -p <phenix_conda> libblas=*=*accelerate 
   ```


### Run Tests

Tests need to be run in an empty directory.

```
 mkdir tests
 cd tests
 qr.test
```

If any of the tests fail, please raise an issue here: [issue tracker](https://github.com/qrefine/qrefine/issues)

### Documentation

Unfortunately the HTML documentation has not been updated yet.
It can be found at: https://qrefine.com/qr.html

### Commandline options

If you want to see the available options and default values please type:

```
 qr.refine --show-defaults
```

### Example

command line options are added like this:

```
qr.refine tests/unit/data_files/helix.pdb engine=mopac clustering=0 gradient_only=1
```
for AQuaRef run
```
qr.aquaref your_pdb.pdb your_map.map
```

### Contact us

The best way to get a hold of us is by sending us an email: qrefine@googlegroups.com

### Developers

- [Pavel Afonine](https://github.com/pafonine)
- [Malgorzata Biczysko](https://github.com/biczysko)
- [Mark Waller](https://github.com/mpwaller)
- [Nigel Moriarty](https://github.com/nwmoriarty)
- [Holger Kruse](https://github.com/hokru)

### Citations:

#### AQuaRef

Roman Zubatyuk, Malgorzata Biczysko, Kavindri Ranasinghe, Nogel W. Moriarty, Hatice Gokcan, Holger Kruse, Billy K. Poon, Paul D. Adams, Mark P. Waller, Adrian E. Roitberg, Olexandr Isayev, and Pavel V. Afonine
AQuaRef: machine learning accelerated quantum refinement of protein structures. 
(2025) Nat. Commun. 16, 9224. 
DOI: [10.1038/s41467-025-64313-1](https://doi.org/10.1038/s41467-025-64313-1)

#### Q|R

Min Zheng, Jeffrey Reimers, Mark P. Waller, and Pavel V. Afonine,
Q|R: Quantum-based Refinement,
(2017) Acta Cryst. D73, 45-52.
DOI: [10.1107/S2059798316019847](http://scripts.iucr.org/cgi-bin/paper?S2059798316019847)

Min Zheng, Nigel W. Moriarty, Yanting Xu, Jeffrey Reimers, Pavel V. Afonine, and Mark P. Waller,
Solving the scalability issue in quantum-based refinement: Q|R#1
(2017) Acta Cryst. D73, 1020-1028.
DOI: [10.1107/S2059798317016746](http://scripts.iucr.org/cgi-bin/paper?S2059798317016746)

Min Zheng, Malgorzata Biczysko, Yanting Xu, Nigel W. Moriarty, Holger Kruse, Alexandre Urzhumtsev, Mark P. Waller, and Pavel V. Afonine,
Including Crystallographic Symmetry in Quantum-based Refinement: Q|R#2
(2020) Acta Cryst. D76, 41-50.
DOI: [10.1107/S2059798319015122](http://scripts.iucr.org/cgi-bin/paper?S2059798319015122)

Lum Wang, Holger Kruse, Oleg V. Sobolev, Nigel W. Moriarty, Mark P. Waller, Pavel V. Afonine, and Malgorzata Biczysko,
Real-space quantum-based refinement for cryo-EM: Q|R#3
(2020) Acta Cryst. D76, 1184-1191.
DOI：[10.1107/S2059798320013194](https://doi.org/10.1107/S2059798320013194)
bioRxiv 2020.05.25.115386.
DOI:[0.1101/2020.05.25.115386](https://www.biorxiv.org/content/10.1101/2020.05.25.115386v1)

Yaru Wang, Holger Kruse, Nigel W. Moriarty, Mark P. Waller, Pavel V. Afonine, and Malgorzata Biczysko,
Optimal clustering for quantum refinement of biomolecular structures: Q|R#4
(2023) Theor. Chem. Acc. 142, 100.
DOI: [10.1007/s00214-023-03046-0](https://doi.org/10.1007/s00214-023-03046-0)
bioRxiv 2022.11.24.517825
DOI:[10.1101/2022.11.24.517825](https://doi.org/10.1101/2022.11.24.517825)

#### Clustering

Min Zheng, Mark P. Waller,
Yoink: An interaction‐based partitioning API,
(2018) Journal of Computational Chemistry, 39, 799–806.
DOI: [10.1002/jcc.25146](https://doi.org/10.1002/jcc.25146)

Min Zheng, Mark P. Waller,
Toward more efficient density-based adaptive QM/MM methods,
(2017)Int J. Quant. Chem e25336
DOI: [10.1002/qua.25336](https://doi.org/10.1002/qua.25336)

Min Zheng, Mark P. Waller, Adaptive QM/MM Methods,
(2016) WIREs Comput. Mol. Sci., 6, 369–385.
DOI: [10.1002/wcms.1255](https://doi.org/10.1002/wcms.1255)

Mark P. Waller, Sadhana Kumbhar, Jack Yang,
A Density‐Based Adaptive Quantum Mechanical/Molecular Mechanical Method
(2014) ChemPhysChem 15, 3218–3225.
DOI: [10.1002/cphc.201402105](https://doi.org/10.1002/cphc.201402105)
