#!/bin/bash


# phenix_prefix env can be found in the installer version, while
# the source install just has phenix..
# check both, first installer

PHENIX_PREFIX=`libtbx.printenv | grep 'PHENIX_PREFIX='`
PHENIX_PREFIX=${PHENIX_PREFIX#"PHENIX_PREFIX="}

if [[ -z ${PHENIX_PREFIX} ]]; then
    # source install
    PHENIX=`libtbx.printenv | grep 'PHENIX='`
    PHENIX=${PHENIX#"PHENIX="}
else
    PHENIX=$PHENIX_PREFIX
fi

# if PHENIX is still empty then no phenix was activated.
if [[ -z ${PHENIX} ]]; then
    echo "activate the phenix installation first!"
    exit
fi

if [[ -d "${PHENIX}/conda_base" ]]; then
 CONDA=${PHENIX}"/conda_base"
 echo "Found developer version of phenix."
 INSTALLER=false
elif  [[ -d "${PHENIX}/conda-meta" ]]; then
 CONDA=${PHENIX}
 INSTALLER=true
 echo "Found installer version of phenix."
else
    echo "something wrong"
    echo "conda: $CONDA"
fi
export PACKAGES=`libtbx.python -c 'import site; print(site.getsitepackages()[0])'`

SKIP=0
if [[ "$1" == "-y" ]] || [[ "$2" == "-y" ]] || [[ "$3" == "-y" ]] ; then
SKIP=1
fi

if [[ "$1" == *"aimnet2"* ]] || [[ "$2" == *"aimnet2"* ]] || [[ "$3" == *"aimnet2"* ]] ; then
    TORCH=true
else
    TORCH=false
fi

CUDA11=false
if [[ "$1" == "cuda11" ]] || [[ "$2" == "cuda11" ]] || [[ "$3" == "cuda11" ]] ; then
CUDA11=true
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

if [[ $SKIP == 0 ]]; then
    read -p "Check above locations/settings. Continue (y/n)?" choice
    case "$choice" in 
      y|Y ) echo "yes";;
      n|N ) echo "no"; exit;;
      * ) echo "invalid"; exit;;
    esac
fi

if [[ ! -e "$PHENIX/build/qrefine_configure.log" ]]; then
    echo "updating phenix/cctbx"
    cd $PHENIX/build
    libtbx.configure qrefine > qrefine_configure.log
fi


QR=$PHENIX/modules/qrefine

# update phenix conda_base
echo "Updating Phenix's conda-base .."
# env update sometimes installs cctbx-base by itself!? changing to install 
# conda env update -p $CONDA --file $QR/config/phenix.yaml > --json > $PHENIX/build/conda_update.json
# conda install -p $CONDA -y --file $QR/config/phenix_req.txt -c conda-forge > $PHENIX/build/conda_qr.out
qrefine.python -m pip install -r $QR/config/phenix.txt > $PHENIX/build/qrefine_req.out

if [ $TORCH == "true" ]; then
echo "Installing pytorch and aimnet2 depedencies .."

    if [ $CUDA11 == "true" ]; then
        # conda install -y -p $CONDA --file $QR/config/aimnet2_cuda11.txt -c nvidia  -c pytorch  -c pyg -c conda-forge > $PHENIX/build/conda_aimnet2.out
        conda env update -p $CONDA --file $QR/config/cuda11.yaml
    else
        # conda install -y -p $CONDA --file $QR/config/aimnet2_cuda12.txt -c conda-forge -c pytorch -c nvidia -c pyg > $PHENIX/build/conda_aimnet2.out
        # conda install -y -p $CONDA --file $QR/config/aimnet2.txt
        conda env update -p $CONDA --file $QR/config/cuda12.yaml
    fi
    
    phenix.python -m pip install git+https://github.com/zubatyuk/aimnet2calc.git 
fi