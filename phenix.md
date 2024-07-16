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
  # installation for AQua (aimnet2, GPU)
  sh build_into_phenix.sh aimnet2
```


#### semi-manual install (developer build)

Below is ideal to update the phenix conda env, e.g. to install aimnet2 components at a later stage or during qrefine development.

```
 # activate phenix
 cd <phenix_installation>
 git clone https://github.com/qrefine/qrefine modules/qrefine
 # request to install aimnet2 for AQua is optional
 sh modules/qrefine/config/update_phenix.sh [aimnet2]
```
