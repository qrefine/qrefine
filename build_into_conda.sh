#!/bin/bash

# HowTo
# 0. download/install latest condaforge
# 1. git clone qrefine repository && cd qrefine
# 2. conda env create -n qrefine -f environment.yaml
# 3. conda activate qrefine
# 3. sh build_into_conda.sh
# 4. set paths with the command given at the end of the script output

export PACKAGES=`python -c 'import site; print(site.getsitepackages()[0])'`
export QREFINE=`pwd`

echo "####################################"
echo "#### QREFINE INSTALLER FOR CONDA ###"
echo "####################################"
echo ""
echo "QR location:  $QREFINE"
echo "python modules location: $PACKAGES"


# add missing modules
mkdir $PACKAGES/modules
cd $PACKAGES/modules

#### QREFINE
mkdir -p $PACKAGES/modules/qrefine
cp -r $QREFINE/* $PACKAGES/modules/qrefine/.

#### PROBE
# dummy directory
mkdir probe

### REDUCE
echo "Downloading reduce"
git clone https://github.com/rlabduke/reduce
cp reduce/reduce_src/hybrid_36_c.c $PACKAGES/iotbx/pdb/hybrid_36_c.c

### set up build dir and exes
echo "CCTBX-install packages"
mkdir $PACKAGES/build
cd $PACKAGES/build
mkdir -p probe/exe/
cp $CONDA_PREFIX/bin/probe probe/exe/
mkdir -p reduce/exe
cp $CONDA_PREFIX/bin/reduce reduce/exe/

### run configure (still in ./build)
libtbx.configure probe qrefine reduce
mmtbx.rebuild_rotarama_cache
mmtbx.rebuild_cablam_cache

echo "run:"
echo "   source $PACKAGES/build/setpaths.sh"
