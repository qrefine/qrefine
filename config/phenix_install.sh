#!/bin/bash

export PHENIX=`pwd`

echo "####################################"
echo "#### QREFINE INSTALLER FOR PHENIX ##"
echo "####################################"
echo ""
echo "Phenix location:  $PHENIX"

# add missing modules
echo "Downloading qrefine"
mkdir $PACKAGES/modules
cd $PACKAGES/modules
git clone https://github.com/qrefine/qrefine
cd qrefine

# update phenix conda_base
echo "Updating Phenix's conda-base"
conda env update -p $PHENIX/conda_base -f config/phenix.yaml

### set up build dir and exes
echo "updating phenix/cctbx"
cd $PACKAGES/build
libtbx.configure qrefine
