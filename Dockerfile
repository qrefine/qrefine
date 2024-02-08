FROM condaforge/mambaforge:23.3.1-1
SHELL ["/bin/bash", "--login", "-c"]

# Base environment setup
COPY environment.yaml .
RUN mamba env create --name cctbx-cuda -f environment.yaml && mamba clean --all

# Activate conda and clean up
RUN echo "conda activate cctbx-cuda" >> ~/.bashrc && echo "export NUMBA_CUDA_USE_NVIDIA_BINDING=1" >> ~/.bashrc
ENV PATH=/opt/conda/envs/cctbx-cuda/bin:${PATH}

# currently conda does not want to install the cuaev version. We use the wheel
RUN mamba remove -p /opt/conda/envs/cctbx-cuda/ torchani && pip install torchani

# Add reduce and probe programs
COPY build_into_conda.sh .
RUN bash build_into_conda.sh

# QREFINE

# install qrefine itself
WORKDIR /opt/conda/envs/cctbx-cuda/lib/python3.10/site-packages
RUN git clone https://github.com/qrefine/qrefine.git qrefine

# clean up qrefine from java
RUN rm -rf qrefine/plugin/yoink

# for debugging
RUN apt-get update && apt-get install -y vim curl

ENV PATH=/opt/conda/envs/cctbx-cuda/bin:/opt/conda/envs/cctbx-cuda/lib/python3.10/site-packages/build/bin:${PATH}
WORKDIR /mnt
