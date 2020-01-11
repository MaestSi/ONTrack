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
##Note: rows starting with '#' are notes for the user, and are ignored by the software
#if do_subsampling_flag <- 1, subsampling of num_fast5_files fast5 files is performed; otherwise set do_subsampling_flag <- 0
do_subsampling_flag <- 0
#num_fast5_files is the number of fast5 files to be subsampled/analysed (if do_subsampling_flag <- 1)                          
num_fast5_files <- 50000
#BC_int are the barcodes used in the experiment
#BC_int <- c("BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08", "BC09", "BC10", "BC11", "BC12")
BC_int <- c("BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08", "BC09", "BC10", "BC11", "BC12")
#barcode kits
#barcode_kits <- c("EXP-NBD103", "EXP-NBD114", "EXP-PBC001", "EXP-PBC096", "SQK-16S024", "SQK-LWB001", "SQK-RAB201", "SQK-RBK001", "SQK-RBK004", "SQK-RLB001")
barcode_kits <- c("EXP-PBC001")
#kit (1D/1D^2 reads/rapid 16S)
kit <- "SQK-LSK109"
#flowcell chemistry (R9.4/R9.5 chemistry)
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
#fixed_lenfil_flag <- 1 if you want to keep reads in the range [amplicon_length - lenfil_tol/2; amplicon_length + lenfil_tol/2]; otherwise set fixed_lenfil_flag <- 1 if you want to keep reads in the range [mean_length -2*sd; mean_length + 2*sd] where mean_length and sd are evaluated on a sample basis
fixed_lenfil_flag <- 0
#if fixed_lenfil_flag <- 1, lenfil_tol [bp] is the size of the window centered in amplicon_length for reads to be kept
lenfil_tol <- 300
#set primers length [bp]
primers_length <- 25
#if disable_porechop_demu_flag <- 1 porechop is only used for adapters trimming and not for doing a second round of demultiplexing; otherwise set disable_porechop_demu_flag <- 0
disable_porechop_demu_flag <- 0
#do_blast_flag <- 1 if you want to perform blast analysis of consensus sequences; otherwise set do_blast_flag <- 0
do_blast_flag <- 1
#do_clustering_flag <- 1 if you want to perform preliminary clustering for getting rid of contaminants; otherwise set do_clustering_flag <- 0
do_clustering_flag <- 1
#num iterations of ONTrack pipeline (must be odd)
num_iterations <- 3
#majority_rule_full_seq_flag <- 1 if you want to pick the full consensus sequence supported by the highest number of iterations; otherwise set majority_rule_full_seq_flag <- 0
majority_rule_full_seq_flag <- 1
########################################################################################################
PIPELINE_DIR <- "/path/to/ONTrack"
#MINICONDA DIR
MINICONDA_DIR <- "/path/to/miniconda3"
#basecaller_dir
BASECALLER_DIR <- "/path/to/ont-guppy-cpu/bin/"
#NCBI nt database
NTDB <- "/path/to/NCBI_nt_db/nt"
########################################################################################################
#load BioStrings package
suppressMessages(library(Biostrings))
#path to ONTrack.R
ONTrack <- paste0(PIPELINE_DIR, "/ONTrack.R")
#path to DecONT.sh
DECONT <- paste0(PIPELINE_DIR, "/decONT.sh")
#path to remove_long_short.pl
remove_long_short <- paste0(PIPELINE_DIR, "/remove_long_short.pl")
#path to subsample fast5
subsample_fast5 <- paste0(PIPELINE_DIR, "/subsample_fast5.sh")
#########################################################################################################
#MAFFT
MAFFT <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/mafft")
#VSEARCH
VSEARCH <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/vsearch")
#NANOPOLISH
NANOPOLISH <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/nanopolish")
#BLASTN
BLASTN <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/blastn")
#EMBOSS cons
CONS <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/cons")
#SEQTK
SEQTK <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/seqtk")
#MINIMAP2
MINIMAP2 <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/minimap2")
#SAMTOOLS
SAMTOOLS <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/samtools")
#PORECHOP
PORECHOP <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/porechop")
#PYCOQC
PYCOQC <- paste0(MINICONDA_DIR, "/envs/ONTrack_env/bin/pycoQC")
