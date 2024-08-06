FROM condaforge/mambaforge:24.3.0-0
SHELL ["/bin/bash", "--login", "-c"]

# Base environment setup
RUN mkdir -p /opt/qrefine
COPY . /opt/qrefine
WORKDIR /opt/qrefine

# clean up qrefine from java
RUN rm -rf plugin/yoink

RUN conda env create --name cctbx-cuda -f environment.yaml && mamba clean --all
RUN conda env update --name cctbx-cuda -f config/cuda12.yaml && mamba clean --all

# debug
RUN conda install vim

# Activate conda and clean up
RUN echo "conda activate cctbx-cuda" >> ~/.bashrc && echo "export NUMBA_CUDA_USE_NVIDIA_BINDING=1" >> ~/.bashrc
ENV PATH=/opt/conda/envs/cctbx-cuda/bin:${PATH}

RUN qrefine.python -m pip install git+https://github.com/zubatyuk/aimnet2calc.git

# run installer
RUN bash build_into_conda.sh

ENV PATH=/opt/conda/envs/cctbx-cuda/bin:/opt/conda/envs/cctbx-cuda/lib/python3.10/site-packages/build/bin:${PATH}
ENV OMP_MAX_ACTIVE_LEVELS=1
ENV OMP_STACKSIZE="4G"
WORKDIR /mnt
