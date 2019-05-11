#
# Copyright 2019 Simone Maestri. All rights reserved.
# Simone Maestri <simone.maestri@univr.it>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

####################################################################################################
#load BioStrings package
suppressMessages(library(Biostrings))
#if do_subsampling_flag <- 1, subsampling of num_fast5_files fast5 files is performed; otherwise set do_subsampling_flag <- 0
do_subsampling_flag <- 0
#num_fast5_files is the number of fast5 files to be subsampled/analysed (if do_subsampling_flag <- 1)                                                      
num_fast5_files <- 100
#BC_int are the barcodes used in the experiment
#BC_int <- c("BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08", "BC09", "BC10", "BC11", "BC12")
BC_int <- c("BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07")
#barcode kits
#barcode_kits <- c("EXP-NBD103", "EXP-NBD114", "EXP-PBC001", "EXP-PBC096", "SQK-16S024", "SQK-LWB001", "SQK-RAB201", "SQK-RBK001", "SQK-RBK004", "SQK-RLB001")
barcode_kits <- c("EXP-PBC001")
#kit (1D/1D^2 reads/rapid 16S)
kit <- "SQK-LSK109"
#flowcell chemistry (R9.4.1 --> FLO-MIN106, R9.5 --> FLO-MIN107 chemistry)
flowcell <- "FLO-MIN106"
#fast_basecalling_flag <- 1 if you want to use the fast basecalling algorithm; otherwise set fast_basecalling_flag <- 0 if you want to use the accurate but slow one (FLO-MIN106 only)
fast_basecalling_flag <- 1
#pair_strands_flag <- 1 if, in case a 1d2 kit and FLO-MIN107 flow-cell have been used, you want to perform 1d2 basecalling; otherwise set pair_strands_flag <- 0
pair_strands_flag <- 0
#save_space_flag <- 1 if you want temporary files to be automatically deleted; otherwise set save_space_flag <- 0
save_space_flag <- 0
#set the maximum number of threads to be used
num_threads <- 30
#set a mean amplicon length [bp]
amplicon_length <- 710
#set primers length [bp]
primers_length <- 25
#doBlast <- 1 if you want to perform blast analysis of consensus sequences; otherwise set doBlast <- 0
doBlast <- 1
########################################################################################################
PIPELINE_DIR <- "/path/to/ONTrack"
#MINICONDA DIR
MINICONDA_DIR <- "/path/to/miniconda3"
#basecaller_dir (v3.0.3+7e7b7d0)
BASECALLER_DIR <- "/path/to/ont-guppy-cpu/bin/"
#NCBI nt database
NTDB <- "/path/to/NCBI_nt_db/nt"
########################################################################################################
#path to ONTrack.R
ONTrack <- paste0(PIPELINE_DIR, "/ONTrack.R")
#num iterations of ONTrack pipeline (must be odd)
num_iterations <- 3
#path to DecONT.sh
DECONT <- paste0(PIPELINE_DIR, "/decONT.sh")
#path to remove_long_short.pl
remove_long_short <- paste0(PIPELINE_DIR, "/remove_long_short.pl")
#path to subsample fast5
subsample_fast5 <- paste0(PIPELINE_DIR, "/subsample_fast5.sh")
#########################################################################################################
#MAFFT (v7.407)
MAFFT <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/mafft")
#VSEARCH (v2.4.4_linux_x86_64)
VSEARCH <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/vsearch")
#NANOPOLISH (v0.11.0)
NANOPOLISH <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/nanopolish")
#BLASTn (v2.2.28+)
BLASTN <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/blastn")
#EMBOSS cons (v6.6.0.0)
CONS <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/cons")
#SEQTK (v1.3-r106)
SEQTK <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/seqtk")
#MINIMAP2 (v2.1.1-r341)
MINIMAP2 <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/minimap2")
#SAMTOOLS (v1.7)
SAMTOOLS <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/samtools")
#PORECHOP (v0.2.3_seqan2.1.1)
PORECHOP <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/porechop")
#PYCOQC (v2.2.3.3)
PYCOQC <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/pycoQC")
