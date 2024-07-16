### CCTBX installation

1.  Clone this repo and enter it's directory.

```
git clone https://github.com/qrefine/qrefine.git qrefine && cd qrefine
```

2.  Use the provided `environment.yaml` to generate a new conda environment called `qrefine`. After installation activate the enviroment. The activation (unless set to automatic) has to be done for every new shell:

    - full installation with pytorch/cuda
      ```
      conda env create -n qrefine -f environment.yaml
      conda activate qrefine
      ```

3.  A couple of configuration steps are needed to setup `qrefine` within `cctbx`. These are automated in a bash script:

```
bash ./build_into_conda.sh
```

4.  run the given `source <path>/setpaths.sh` command at the end of the script. This needs to sourced for every new shell.



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

