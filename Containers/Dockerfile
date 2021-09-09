#INFORMATION
#MAINTAINER Julie Orjuela & Francois Sabot (Docker)
#version="1.0"
#software="all Culebront dependencies"
#description="All dependencies used to launch CulebrONT in LOCAL mode"
#website="https://culebront-pipeline.readthedocs.io/en/latest/"


FROM ubuntu:bionic
USER root

#ENVIRONMENT
RUN export CULEBRONT="/usr/local/CulebrONT_pipeline/"
RUN export PATH="/usr/local/miniconda/miniconda3/bin:$PATH"
RUN export PATH="/usr/local/miniconda/miniconda3/envs/snakemake/bin:$PATH"

RUN apt update
RUN apt install -y git vim curl wget less locate build-essential openssh-server python3-all-dev python3-pip python3-venv graphviz build-essential libssl-dev uuid-dev libgpgme11-dev squashfs-tools libseccomp-dev pkg-config
