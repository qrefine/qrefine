### Phenix installation
Note: you may need to use sudo depending on the permissions of your Phenix installation.

#### fully automated, recommended

The following should work for both the typical user install and a developer source build.
The installation type is auto-detected.
```
  # activate phenix, e.g. with
  source path_to_phenix/phenix_env.sh
  cd qrefine
  # choose one of the below:
  # default, minimal installation
  sh build_into_phenix.sh 
  # installation for AQua (aimnet2, GPU) CUDA12
  sh build_into_phenix.sh aimnet2
  # installation for AQua (aimnet2, GPU) and CUDA11
  sh build_into_phenix.sh aimnet2
```


#### semi-manual install (developer build)

Below is ideal to update the phenix conda env, e.g. to install aimnet2 components at a later stage or during qrefine development.

```
 # activate phenix
 cd <phenix_installation>
 git clone https://github.com/qrefine/qrefine modules/qrefine
 # optional arguments to install AQua (aimnet2, GPU), flag to request cuda 11 instead of 12 and "-y" skips the confirmation question.
 sh modules/qrefine/config/update_phenix.sh [aimnet2] [cuda11] [-y]
```
