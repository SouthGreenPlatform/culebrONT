#INFORMATION
#MAINTAINER Julie Orjuela & Francois Sabot (Docker)
#version="1.0"
#software="all Culebront dependencies"
#description="All dependencies used to launch CulebrONT in LOCAL mode"
#website="https://culebront-pipeline.readthedocs.io/en/latest/"

FROM ubuntu:bionic
USER root

# ENVIRONMENT
ENV CULEBRONT="/usr/local/CulebrONT_pipeline"
ENV PATH="/usr/local/miniconda/miniconda3/bin:$PATH"
ENV PATH="/usr/local/miniconda/miniconda3/envs/snakemake/bin:$PATH"

# INSTALL, Global
RUN apt update
RUN DEBIAN_FRONTEND=noninteractive apt -y install keyboard-configuration git vim curl wget less locate openssh-server python3-all-dev python3-pip python3-venv graphviz build-essential libssl-dev uuid-dev libgpgme11-dev squashfs-tools libseccomp-dev pkg-config rsync gzip libcairo2-dev libxt-dev zlib1g-dev cmake cryptsetup
RUN ln -fs /usr/share/zoneinfo/Europe/Paris /etc/localtime
RUN apt install -y tzdata
RUN dpkg-reconfigure --frontend noninteractive tzdata

# INSTALL, R for reporting
RUN apt install -y r-recommended r-doc-html util-linux zip

# Install for Python
RUN echo 'export LC_ALL=C.UTF-8' >> /environment
RUN echo 'export LANG=C.UTF-8' >> /environment
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
RUN python3 -m pip install PyYAML
RUN python3 -m pip install pandas
RUN python3 -m pip install seaborn
RUN python3 -m pip install matplotlib
RUN python3 -m pip install tabulate
RUN python3 -m pip install rpy2
RUN python3 -m pip install ipython
RUN python3 -m pip install biopython
RUN python3 -m pip install tqdm 
RUN python3 -m pip install pyyaml 
RUN python3 -m pip install pysam 
RUN python3 -m pip install docopt==0.6.2
RUN python3 -m pip install numpy
RUN python3 -m pip install argparse

# INSTALL miniconda & conda dependencies
RUN mkdir /usr/local/miniconda && cd /usr/local/miniconda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN chmod +x Miniconda3-latest-Linux-x86_64.sh
RUN ./Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/miniconda/miniconda3
RUN cd /usr/local/miniconda/miniconda3/bin
ENV PATH="/usr/local/miniconda/miniconda3/bin:$PATH"
RUN echo $PATH

RUN /usr/local/miniconda/miniconda3/bin/conda config --add channels defaults
RUN /usr/local/miniconda/miniconda3/bin/conda config --add channels bioconda
RUN /usr/local/miniconda/miniconda3/bin/conda config --add channels conda-forge

RUN /usr/local/miniconda/miniconda3/bin/conda create -n snakemake
RUN /usr/local/miniconda/miniconda3/bin/conda install mamba
RUN /usr/local/miniconda/miniconda3/bin/mamba install snakemake

# installing Singularity
RUN wget https://dl.google.com/go/go1.16.4.linux-amd64.tar.gz
RUN tar -C /usr/local -xzvf go1.16.4.linux-amd64.tar.gz 
RUN rm go1.16.4.linux-amd64.tar.gz 
ENV PATH=/usr/local/go/bin:$PATH

## Install Singularity itself
RUN wget https://github.com/sylabs/singularity/releases/download/v3.8.1/singularity-ce-3.8.1.tar.gz 
RUN tar -xzf singularity-ce-3.8.1.tar.gz
ENV PKG_CONFIG_PATH=/usr/lib/pkgconfig:/usr/lib/x86_64-linux-gnu/pkgconfig/
RUN cd /singularity-ce-3.8.1 && ./mconfig -p /usr/local && make -C builddir && make -C builddir install

# Clone CulebrONT ... yeah!
RUN mkdir -p /usr/local
RUN git clone  --recursive https://github.com/SouthGreenPlatform/CulebrONT_pipeline.git /usr/local/CulebrONT_pipeline
RUN cd /usr/local/CulebrONT_pipeline/

# Download build singularity
RUN cd /usr/local/CulebrONT_pipeline/Containers && wget --no-check-certificate -rm -nH --cut-dirs=2 --reject="index.html*" --no-parent https://itrop.ird.fr/culebront_utilities/singularity_build/

# sed toolpath.yaml
RUN cd /usr/local/CulebrONT_pipeline/ && sed -i "s|/path/to/|/usr/local/CulebrONT_pipeline/|g" /usr/local/CulebrONT_pipeline/tools_path.yaml
ENV PATH="/usr/local/CulebrONT_pipeline":$PATH

# Adding singularity binding
ENV SINGULARITY_BIND=/usr/local/CulebrONT_pipeline/

# horrible sed to rule_graph
RUN cd /usr/local/CulebrONT_pipeline/ && sed -i "s|{params.cmd}|snakemake -s /usr/local/CulebrONT_pipeline/Snakefile --use-singularity --rulegraph |g" /usr/local/CulebrONT_pipeline/Snakefile


