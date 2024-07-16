#!/bin/bash
QREFINE=`pwd`

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

TORCH=false
if [[ "$1" == *"aimnet2"* ]]; then
    TORCH=true
fi

echo "######################################"
echo "#### QREFINE INSTALLER FOR PHENIX  ###"
echo "######################################"
echo ""
echo "QREFINE location: $QREFINE"
echo "Phenix location: $PHENIX"
echo "Aimnet2 install?: $TORCH"
echo ""

if [[ ! -d "${PHENIX}/modules" ]]; then
    mkdir $PHENIX/modules
fi
if [[ ! -d "${PHENIX}/build" ]]; then
    mkdir $PHENIX/build
fi
# mkdir $PHENIX/modules $PHENIX/build

#### QREFINE
echo "Copying qrefine module ..."
mkdir -p $PHENIX/modules/qrefine
cp -r $QREFINE/../qrefine/* $PHENIX/modules/qrefine/.

# dont know if that should be here.
### run configure (still in ./build)
# cd $PHENIX/build
# echo "Running libtbx.configure for qrefine in ./build ..."
# libtbx.configure qrefine &> qrefine_configure.log
# echo "Running: source $PHENIX/build/setpaths.sh"
# source $PHENIX/build/setpaths.sh

# install packages
echo "Updating phenix conda with QR packages ..."
echo "  Running $QREFINE/config/update_phenix.sh"
cd $PHENIX
if [[ $TORCH == "true" ]]; then
    sh $QREFINE/config/update_phenix.sh aimnet2
    
else
    sh $QREFINE/config/update_phenix.sh
fi
    

echo "Setup QR+Phenix in the future with:"
echo "   source $PHENIX/build/setpaths.sh"
