FROM continuumio/miniconda3

########### set variables
ENV DEBIAN_FRONTEND noninteractive
ENV GUPPY_VERSION=6.0.1

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

################################################################ install Guppy
RUN mkdir /home/tools/guppy_inst
WORKDIR /home/tools/guppy_inst

# installation
RUN wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_${GUPPY_VERSION}_linux64.tar.gz \
  && tar -xf ont-guppy-cpu_${GUPPY_VERSION}_linux64.tar.gz
ENV PATH $PATH:/home/tools/guppy_inst/ont-guppy-cpu/bin/

############################################################ install ONTrack
WORKDIR /home/tools/

RUN git clone https://github.com/MaestSi/ONTrack.git
WORKDIR /home/tools/ONTrack
RUN chmod 755 * \
  && ./install.sh

WORKDIR /home/
