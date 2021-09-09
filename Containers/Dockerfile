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

#INSTALL, Global
RUN apt update
RUN apt install -y git vim curl wget less locate build-essential openssh-server python3-all-dev python3-pip python3-venv graphviz build-essential libssl-dev uuid-dev libgpgme11-dev squashfs-tools libseccomp-dev pkg-config

RUN export DEBIAN_FRONTEND=noninteractive
RUN ln -fs /usr/share/zoneinfo/Europe/Paris /etc/localtime
RUN apt install -y tzdata
RUN dpkg-reconfigure --frontend noninteractive tzdata

#INSTALL, R for reporting
RUN apt install -y r-recommended r-doc-html util-linux zip

#Install for Python
RUN echo 'export LC_ALL=C.UTF-8' >> /environment
RUN echo 'export LANG=C.UTF-8' >> /environment
RUN export LC_ALL=C.UTF-8
RUN export LANG=C.UTF-8
RUN apt install -y rsync gzip libcairo2-dev libxt-dev
RUN apt install -y build-essential zlib1g-dev cmake curl
RUN python3 -m pip install PyYAML
RUN python3 -m pip install pandas
RUN python3 -m pip install seaborn
RUN python3 -m pip install matplotlib
RUN python3 -m pip install tabulate
RUN python3 -m pip install rpy2
RUN python3 -m pip install ipython
RUN python3 -m pip install biopython
RUN python3 -m pip install tqdm pyyaml pysam docopt==0.6.2
RUN python3 -m pip install numpy
RUN python3 -m pip install argparse

#INSTALL miniconda & conda dependencies
RUN mkdir /usr/local/miniconda && cd /usr/local/miniconda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN chmod +x Miniconda3-latest-Linux-x86_64.sh
RUN ./Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/miniconda/miniconda3
RUN cd /usr/local/miniconda/miniconda3/bin
RUN export PATH="/usr/local/miniconda/miniconda3/bin:$PATH" >> $SINGULARITY_ENVIRONMENT

