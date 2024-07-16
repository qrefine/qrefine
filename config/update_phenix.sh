#!/bin/bash


# phenix_prefix env can be found in the installer version, while
# the source install just has phenix..
# check both, first installer

PHENIX_PREFIX=`libtbx.printenv | grep 'PHENIX_PREFIX='`
PHENIX_PREFIX=${PHENIX_PREFIX#"PHENIX_PREFIX="}

if [[ -z ${PHENIX_PREFIX} ]]; then
    # source install
    # echo "debug: hit"
    PHENIX=`libtbx.printenv | grep 'PHENIX='`
    PHENIX=${PHENIX#"PHENIX="}
    libtbx.printenv | grep 'PHENIX'
else
    PHENIX=$PHENIX_PREFIX
fi

# if PHENIX is still empty then no phenix was activated.
if [[ -z ${PHENIX} ]]; then
    echo "activate the phenix installation first!"
    exit
fi

# echo "DEBUG"
# echo "${PHENIX}/conda_base"

if [[ -d "${PHENIX}/conda_base" ]]; then
 CONDA=${PHENIX}"/conda_base"
 echo "Found developer version of phenix."
elif  [[ -d "${PHENIX}/conda-meta" ]]; then
 CONDA=${PHENIX}
 INSTALLER=true
 echo "Found installer version of phenix."
else
    echo "something wrong"
    echo "conda: $CONDA"
fi
export PACKAGES=`libtbx.python -c 'import site; print(site.getsitepackages()[0])'`

if [[ "$1" == *"aimnet2"* ]]; then
    TORCH=true
else
    TORCH=false
fi

ARM64=false
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
# echo "####################################"
# echo "#### QREFINE UPDATER FOR PHENIX   ##"
# echo "####################################"
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

# run configure. Skip if qrefine_configure.log exists already
if [[ ! -e "$PHENIX/build/qrefine_configure.log" ]]; then
    echo "updating phenix/cctbx"
    cd $PHENIX/build
    libtbx.configure qrefine > qrefine_configure.log
fi
