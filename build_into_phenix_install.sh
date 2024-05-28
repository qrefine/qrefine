#!/bin/bash
QREFINE=`pwd`
PHENIX=$PHENIX_PREFIX

# if [[ "$1" == "" ]]; then
#     echo "usage: sh ./build_into_phenix_install.sh <phenix_location>"
#     echo ""
#     echo "E.g. sh ./build_into_phenix_install.sh /opt/phenix-1.21.1"
# fi

if [[ -z ${PHENIX_PREFIX} ]]; then
    echo "source phenix_env.sh first!"
    exit
fi

TORCH=false
if [[ "$1" == *"aimnet2"* ]]; then
    TORCH=true
fi

echo "##################################################"
echo "#### QREFINE INSTALLER FOR PHENIX USER INSTALL ###"
echo "##################################################"
echo ""
echo "QREFINE location: $QREFINE"
echo "Phenix location: $PHENIX"
echo "Aimnet2 install?: $TORCH"
echo ""

mkdir $PHENIX/modules $PHENIX/build

#### QREFINE
echo "Copying qrefine module ..."
cp -r $QREFINE/../qrefine $PHENIX/modules/qrefine

### run configure (still in ./build)
cd $PHENIX/build
echo "Running libtbx.configure for qrefine in ./build ..."
libtbx.configure qrefine &> qrefine_configure.log
echo "Running: source $PHENIX/build/setpaths.sh"
source $PHENIX/build/setpaths.sh

# install packages
echo "Updating phenix conda with QR packages ..."
echo "  Running $QREFINE/config/update_phenix.sh"
cd $PHENIX
if [[ $TORCH == "true" ]]; then
    echo "that"
    sh $QREFINE/config/update_phenix2.sh aimnet2
    
else
    echo "this"
    sh $QREFINE/config/update_phenix2.sh
fi
    

echo "Setup QR+Phenix again with:"
echo "   source $PHENIX/build/setpaths.sh"
