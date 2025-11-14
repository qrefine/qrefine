### Phenix installation
Note: you may need to use sudo depending on the permissions of your Phenix installation.

#### fully automated, recommended

AQuaRef is fully integrated into Phenix  starting dev-5395 version.

It is recommended to install recent Phenix version available at https://phenix-online.org/download/nightly_builds.cgi

For using AQuaRef with GPU it is needed to download version for CUDA 11 or CUDA 12, depending on the user system set-up.

Usage instructions are available at https://phenix-online.org/version_docs/2.0-5867/reference/AQuaRef.html



#### semi-manual install (developer build) - NOT TESTED ! with recent Phenix

Below is ideal to update the phenix conda env, e.g. to install aimnet2 components at a later stage or during qrefine development.

```
 # activate phenix
 cd <phenix_installation>
 git clone https://github.com/qrefine/qrefine modules/qrefine
 # optional arguments to install AQua (aimnet2, GPU), flag to request cuda 11 instead of 12 and "-y" skips the confirmation question.
 sh modules/qrefine/config/update_phenix.sh [aimnet2] [cuda11] [-y]
```
