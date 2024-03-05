#!/bin/bash

export PHENIX=`pwd`

if [[ "$1" == *"pytorch"* ]]; then
    TORCH=true
fi

echo "####################################"
echo "#### QREFINE INSTALLER FOR PHENIX ##"
echo "####################################"
echo ""
echo "Phenix location :  $PHENIX"
echo "Pytorch/Aimnet2 :  $TORCH"
echo ""

QR=$PHENIX/modules/qrefine
# update phenix conda_base
echo "Updating Phenix's conda-base .."
conda env update -p $PHENIX/conda_base -f $QR/config/phenix.yaml 
if [ $TORCH ]; then
echo "Installing pytorch and aimnet2 depedencies .."
    conda env update -p $PHENIX/conda_base -f $QR/config/aimnet2.yaml 
fi

### set up build dir and exes
echo "updating phenix/cctbx"
cd $PHENIX/build
libtbx.configure qrefine > qrefine_configure.log
