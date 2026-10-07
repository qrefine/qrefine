### AIMNET2-based AQua plugins

Two AIMNet2 python packages are supported:

* `aimnet` from [isayevlab/aimnetcentral](https://github.com/isayevlab/aimnetcentral) (current, preferred)
* `aimnet2calc` from [zubatyuk/aimnet2calc](https://github.com/zubatyuk/aimnet2calc) (legacy)

If both are installed, `aimnet` is used.

#### Using the `aimnet` package

```
  qrefine.python -m pip install git+https://github.com/isayevlab/aimnetcentral.git
```

Notes:

* The AQuaRef model `aimnet2-qr` is a legacy TorchScript model that is not in
  the `aimnet` model registry. Q|R looks for the file
  `aimnet2-qr_b97md4_qzvp_2.jpt` in the `aimnet` cache directory
  (`$AIMNET_CACHE_DIR`, default `~/.cache/aimnet`) and downloads it from the
  original model zoo on first use. To work offline, copy the file there yourself.
* On CUDA, the `aimnet` neighbour-list kernels are compiled with Triton, which
  needs the CUDA `ptxas` binary. Q|R sets `TRITON_PTXAS_PATH` to the `ptxas`
  shipped in the Phenix conda environment (from the `cuda-nvcc-tools` package)
  if it is not already set.

#### Using the legacy `aimnet2calc` package

in general please follow the general installation notes. To update your qrefine installation use:

For cctbx

```
  conda env update -f config/aimnet2.yaml
  qrefine.python -m pip install git+https://github.com/zubatyuk/aimnet2calc.git
```

For Phenix either use the provided script `config/update_phenix.sh` as described above or install into the phenix conda env like  this:

```
  conda env update -f -p /path/to/phenix/conda_base config/aimnet2.yaml
  qrefine.python -m pip install git+https://github.com/zubatyuk/aimnet2calc.git
```

#### Performance
Set the following in your Terminal for optimal performance. Save it to your .bashrc (or similar).

```
  export NUMBA_CUDA_USE_NVIDIA_BINDING=1
```

To check if the cuda components are working run:

```
  qrefine.python -c "import numba.cuda; print(numba.cuda.is_available())"
  qrefine.python -c "import torch; print(torch.cuda.is_available())"
```


### torchani

torchani is not installed by default. Try to install it via conda

```
conda install -p <phenix_conda> torchani -c conda-forge
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
