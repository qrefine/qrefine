#!/bin/bash

export PHENIX=`pwd`

#ToDo: parse cmd line arguments to switch on/off aimnet and set phenix dir
# Call getopt to validate the provided input. 
options=$(getopt --long pytorch: -- "$@")
[ $? -eq 0 ] || { 
    echo "Incorrect options provided"
    exit 1
}
eval set -- "$options"
while true; do
    case "$1" in
    --pytorch)
        TORCH=true
    esac
    shift
done



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
