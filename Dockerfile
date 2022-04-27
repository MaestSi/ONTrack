FROM continuumio/miniconda3

########### set variables
ENV DEBIAN_FRONTEND noninteractive

########## generate working directories
RUN mkdir /home/tools

######### dependencies
RUN apt-get update -qq \
    && apt-get install -y \
    build-essential \
    wget \
    unzip \
    bzip2 \
    git \
    libidn11* \
    nano \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

############################################################ install ONTrack
WORKDIR /home/tools/

RUN git clone https://github.com/MaestSi/ONTrack.git
WORKDIR /home/tools/ONTrack
RUN chmod 755 *

RUN sed -i 's/PIPELINE_DIR <- .*/PIPELINE_DIR <- \"\/home\/tools\/ONTrack\/\"/' config_MinION_mobile_lab.R
RUN sed -i 's/MINICONDA_DIR <- .*/MINICONDA_DIR <- \"\/opt\/conda\/\"/' config_MinION_mobile_lab.R

RUN conda config --add channels bioconda && \
conda config --add channels anaconda && \
conda config --add channels r && \
conda config --add channels conda-forge
RUN conda create -n ONTrack_env -c bioconda bioconductor-biostrings
RUN conda install -n ONTrack_env python blast emboss vsearch seqtk mafft minimap2 samtools=1.15 nanopolish bedtools ncurses ont_vbz_hdf_plugin
ENV HDF5_PLUGIN_PATH=/opt/conda/envs/ONTrack_env/hdf5/lib/plugin
RUN /opt/conda/envs/ONTrack_env/bin/pip install pycoQC

WORKDIR /home/
