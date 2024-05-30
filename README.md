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
 - cctbx user
 - phenix installer user
 - phenix source user

Requirements:
 - python >= 3.9
 - conda (miniconda)

We require a minimum python version of 3.9! 

### Phenix route
Note: you may need to use sudo depending on the permissions of your Phenix installation.

#### Phenix source

```
 source phenix_env.sh
 git clone https://github.com/qrefine/qrefine modules/qrefine
 # request to install aimnet2 is optional
 sh modules/qrefine/config/update_phenix.sh [aimnet2]
 qrefine.python -m pip install git+https://github.com/zubatyuk/aimnet2calc.git 
```

#### Phenix installer 

```
  source path_to_phenix/phenix_env.sh
  git clone https://github.com/qrefine/qrefine 
  cd qrefine
  # request to install aimnet2 is optional
  sh build_into_phenix_install.sh [aimnet2]
```

#### cctbx (open-source)


2.  Clone this repo and enter it's directory.

```
git clone https://github.com/qrefine/qrefine.git qrefine && cd qrefine
```

3.  Use the provided `environment.yaml` to generate a new conda environment called `qrefine`. After installation activate the enviroment. The activation (unless set to automatic) has to be done for every new shell:

    - full installation with pytorch/cuda
      ```
      conda env create -n qrefine -f environment.yaml
      conda activate qrefine
      ```

4.  A couple of configuration steps are needed to setup `qrefine` within `cctbx`.

```
bash ./build_into_conda.sh
```

5.  run the given `source <path>/setpaths.sh` command at the end of the script. Also this needs to sourced for every new shell.


#### AIMNET2 (and torchani) plugins

For cctbx:

```
  conda env update -f config/aimnet2.yaml
  qrefine.python -m pip install git+https://github.com/zubatyuk/aimnet2calc.git
```

For Phenix either use the provided script `config/update_phenix.sh` as described above or install into the phenix conda env

```
  # identify your phenix directory
  conda env list
  conda env update -f -p /path/to/phenix/conda_base config/aimnet2.yaml
  qrefine.python -m pip install git+https://github.com/zubatyuk/aimnet2calc.git
```

Set the following in your Terminal for optimal performance. Save it to your .bashrc (or similar).

```
  export NUMBA_CUDA_USE_NVIDIA_BINDING=1
```

To check if the cuda components are working run:

```
  qrefine.python -c "import numba.cuda; print(numba.cuda.is_available())"
  qrefine.python -c "import torch; print(torch.cuda.is_available())"
```

(Optional) Check if the cuda AEV version of torchani was installed:

    ```
    mamba list | grep torchani
    ls $(qrefine.python -c 'import site; print(site.getsitepackages()[0])')/torchani/cuaev
    ```

    It should say `torchani=*=cuda...` and the `cuaev` directory is present. If not you can try the the pip/wheel installation:

    ```
    mamba remove torchani
    pip install torchani
    ```

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

### conda packages (work in progres)

[![Anaconda](https://anaconda.org/qrefine/qrefine/badges/latest_release_date.svg)](https://anaconda.org/qrefine/qrefine)
[![Anaconda](https://anaconda.org/qrefine/qrefine/badges/version.svg)](https://anaconda.org/qrefine/qrefine)
[![Anaconda](https://anaconda.org/qrefine/qrefine/badges/platforms.svg)](https://anaconda.org/qrefine/qrefine)

A conda package is provided for qrefine. We currently make use of the nightly build of cctbx. Use `conda >=23.10` or `mamba`
**The conda-forge setup https://github.com/conda-forge/miniforge#miniforge3 is recommended.**

During quick development cycles the conda packages will lag behind as they are build manually.

```
conda create -n QR qrefine -c qrefine -c cctbx-nightly
conda activate QR
qr.test
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

Can be found at: https://qrefine.com/qr.html

### Commandline options

If you want to see the available options and default values please type:

```
 qr.refine --show-defaults
```

### Example

command line options are added ike this:

```
qr.refine tests/unit/data_files/helix.pdb engine=mopac clustering=0 gradient_only=1
```

### Contact us

The best way to get a hold of us is by sending us an email: qrefine@googlegroups.com

### Developers

- [Pavel Afonine](https://github.com/pafonine)
- [Malgorzata Biczysko](https://github.com/biczysko)
- [Mark Waller](https://github.com/mpwaller)
- [Nigel Moriarty](https://github.com/nwmoriarty)
- [Holger Kruse](https://github.com/hokru)
- [Min Zheng](https://github.com/zhengmin317)
- [Lum Wang](https://github.com/Mooooony)

### Citations:

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
