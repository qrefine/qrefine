### Phenix installation
Note: you may need to use sudo depending on the permissions of your Phenix installation.

#### Phenix installer (typical user installation)

```
  source path_to_phenix/phenix_env.sh
  git clone https://github.com/qrefine/qrefine 
  cd qrefine
  # choose one of the below:
  # default, minimal installation
  sh build_into_phenix_install.sh 
  # installation including aimnet2
  sh build_into_phenix_install.sh aimnet2
```


#### Phenix source (developer installation)

```
 # activate phenix
 cd <phenix_installation>
 git clone https://github.com/qrefine/qrefine modules/qrefine
 # request to install aimnet2 is optional
 sh modules/qrefine/config/update_phenix.sh [aimnet2]
```
