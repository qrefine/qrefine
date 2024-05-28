#!/bin/bash

if [[ -z ${PHENIX_PREFIX} ]]; then
    echo "source phenix_env.sh first!"
    exit
fi

PHENIX=${PHENIX_PREFIX}

# execute from phenix base dir
if [[ -d "${PHENIX_PREFIX}/conda_base" ]]; then
 CONDA=${PHENIX_PREFIX}"/conda_base"
 echo "Found developer version of phenix."
elif  [[ -d "${PHENIX_PREFIX}/conda-meta" ]]; then
 CONDA=${PHENIX_PREFIX}
 INSTALLER=true
 echo "Found installer version of phenix."
else
    echo "something wrong"
    echo "conda: $CONDA"
fi
export PACKAGES=`phenix.python -c 'import site; print(site.getsitepackages()[0])'`

if [[ "$1" == *"aimnet2"* ]]; then
    TORCH=true
else
    TORCH=false
fi

if [[ $(uname -m) == "arm64" ]]; then
    ARM64=true
    echo ""
    echo "Found osx-arm64 (Apple Silicon) machine."
    echo "installation of aimnet2 will be skipped if requested"
    echo "Contact a developer to use aimnet2 on Apple Silicon"
    echo ""
    TORCH=false
fi

echo ""
echo "####################################"
echo "#### QREFINE INSTALLER FOR PHENIX ##"
echo "####################################"
echo ""
echo "Phenix location :  $PHENIX"
echo "Phenix conda    :  $CONDA"
echo "Pytorch/Aimnet2 :  $TORCH"
echo "osx-arm64       :  $ARM64"
echo ""
echo ""
read -p "Check above locations/settings. Continue (y/n)?" choice
case "$choice" in 
  y|Y ) echo "yes";;
  n|N ) echo "no"; exit;;
  * ) echo "invalid"; exit;;
esac

QR=$PHENIX/modules/qrefine
# update phenix conda_base
echo "Updating Phenix's conda-base .."

conda env update -p $CONDA -f $QR/config/phenix.yaml 
if [ $TORCH == "true" ]; then
echo "Installing pytorch and aimnet2 depedencies .."
    conda env update -p $CONDA  -f $QR/config/aimnet2.yaml
    phenix.python -m pip install git+https://github.com/zubatyuk/aimnet2calc.git 
fi

# libtbx.configure expects iotbx in the python site-packages dir
if [[ $INSTALLER == "false" ]]; then
    cp -r $PHENIX/modules/cctbx_project/iotbx $PACKAGES/.
fi

# run configure. Skip if qrefine_configure.log exists, this means this script was called by build_into_phenix_install.sh
if [[ ! -e "$PHENIX/build/qrefine_configure.log" ]]; then
    echo "updating phenix/cctbx"
    cd $PHENIX/build
    libtbx.configure qrefine > qrefine_configure.log
fi
