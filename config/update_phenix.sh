#!/bin/bash

export PHENIX=`pwd`
if [ -d "${PHENIX}/conda_base" ]; then
 export PHENIX="$PHENIX"
elif  [ -f "${PHENIX}/qr.py" ]; then
 export PHENIX="${PHENIX}/../../"
fi
export PACKAGES=`phenix.python -c 'import site; print(site.getsitepackages()[0])'`

if [[ "$1" == *"aimnet2"* ]]; then
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

# libtbx.configure expects iotbx in the python site-packages dir
cp -r $PHENIX/modules/cctbx_project/iotbx $PACKAGES/.

### set up build dir and exes
echo "updating phenix/cctbx"
cd $PHENIX/build
libtbx.configure qrefine > qrefine_configure.log
