#!/bin/bash

export PHENIX=`pwd`

#ToDo: parse cmd line arguments to switch on/off aimnet and set phenix dir
# Call getopt to validate the provided input. 
# if [ -z $1]; then
if [[ "$1" == *"pytorch"* ]]; then
    echo "installing aimnet2 dependencies!"
    TORCH=true
fi
# fi
exit

echo "####################################"
echo "#### QREFINE INSTALLER FOR PHENIX ##"
echo "####################################"
echo ""
echo "Phenix location:  $PHENIX"

# update phenix conda_base
echo "Updating Phenix's conda-base"
conda env update -p $PHENIX/conda_base -f config/phenix.yaml > conda_update1.out
if [ TORCH ]; then
    conda env update -p $PHENIX/conda_base -f config/aimnet2.yaml > conda_update1.out
fi
conda clean --all

### set up build dir and exes
echo "updating phenix/cctbx"
cd $PACKAGES/build
libtbx.configure qrefine
