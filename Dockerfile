FROM condaforge/mambaforge:23.3.1-1
SHELL ["/bin/bash", "--login", "-c"]

# Base environment setup
COPY environment.yaml .
RUN mamba env create --name cctbx-cuda -f environment.yaml && mamba clean --all


# Activate conda and clean up
RUN echo "conda activate cctbx-cuda" >> ~/.bashrc && echo "export NUMBA_CUDA_USE_NVIDIA_BINDING=1" >> ~/.bashrc
ENV PATH=/opt/conda/envs/cctbx-cuda/bin:${PATH}

# Add restraints libraries
WORKDIR /opt/conda/envs/cctbx-cuda/lib/python3.10/site-packages
RUN mkdir chem_data &&\
    cd chem_data &&\
    svn --quiet --non-interactive --trust-server-cert co svn://svn.code.sf.net/p/geostd/code/trunk geostd &&\
    svn --quiet --non-interactive --trust-server-cert co https://github.com/rlabduke/mon_lib.git/trunk mon_lib &&\
    svn --quiet --non-interactive --trust-server-cert co https://github.com/rlabduke/reference_data.git/trunk/Top8000/Top8000_rotamer_pct_contour_grids rotarama_data &&\
    rm -rf rotarama_data/.svn &&\
    svn --quiet --non-interactive --trust-server-cert --force co https://github.com/rlabduke/reference_data.git/trunk/Top8000/Top8000_ramachandran_pct_contour_grids rotarama_data &&\
    svn --quiet --non-interactive --trust-server-cert co https://github.com/rlabduke/reference_data.git/trunk/Top8000/Top8000_cablam_pct_contour_grids cablam_data &&\
    cd ..

# Add reduce and probe programs
RUN mkdir modules
WORKDIR /opt/conda/envs/cctbx-cuda/lib/python3.10/site-packages/modules
RUN svn --quiet --non-interactive --trust-server-cert co https://github.com/rlabduke/probe.git/trunk probe &&\
    cd probe &&\
    make &&\
    cp hybrid_36_c.c /opt/conda/envs/cctbx-cuda/lib/python3.10/site-packages/iotbx/pdb/hybrid_36_c.c

RUN svn --quiet --non-interactive --trust-server-cert co https://github.com/rlabduke/reduce.git/trunk reduce &&\
    cd reduce &&\
    make

# QREFINE

# install qrefine itself
RUN git clone https://github.com/qrefine/qrefine.git qrefine

# clean up qrefine from java
RUN rm -rf qrefine/plugin/yoink

WORKDIR /opt/conda/envs/cctbx-cuda/lib/python3.10/site-packages
# ## add cctbx config directory
RUN mkdir build &&\
    cd build &&\
    mkdir -p probe/exe/ && cp ../modules/probe/probe probe/exe/. &&\
    libtbx.configure probe &&\
    libtbx.configure qrefine &&\
    libtbx.configure reduce &&\
    mmtbx.rebuild_rotarama_cache &&\
    mmtbx.rebuild_cablam_cache

RUN apt-get update && apt-get install -y vim curl

ENV PATH=/opt/conda/envs/cctbx-cuda/bin:/opt/conda/envs/cctbx-cuda/lib/python3.10/site-packages/build/bin:${PATH}
WORKDIR /mnt
