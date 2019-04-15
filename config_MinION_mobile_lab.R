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
#if do_subsampling_flag is equal to 1, subsampling of num_kilo_reads*1000 reads is performed
do_subsampling_flag <- 0
#num_kilo reads is the number of thousands of reads to be subsampled/analysed (if do_subsampling_flag == 1)                                                      
num_kilo_reads <- 5
#BC_int are the barcodes used in the experiment
#BC_int <- c("BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08", "BC09", "BC10", "BC11", "BC12")
BC_int <- c("BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07")
#barcode kits
barcode_kits <- c("EXP-NBD103", "EXP-NBD104", "EXP-NBD114", "EXP-PBC001", "EXP-PBC096", "SQK-16S024", "SQK-LWB001", "SQK-PBK004", "SQK-RAB201", "SQK-RAB204", "SQK-RBK001", "SQK-RBK004", "SQK-RLB001", "SQK-RPB004", "VMK-VMK001")
#kit (1D/1D^2 reads/rapid 16S)
kit <- "SQK-LSK108"
#pair_strands_flag is a variable that controls, in case a 1d2 kit has been used, if 1d2 basecalling has to be performed (when set equal to 1 ) or not
pair_strands_flag <- 0
#flowcell chemistry (R9.4/R9.5 chemistry)
flowcell <- "FLO-MIN106"
#save_space_flag is a variable that controls deletion of temporary files; if it is equal to 1, temporary files are deleted
save_space_flag <- 0
#set the maximum number of threads to be used
num_threads <- 30
#set a mean amplicon length [bp]
amplicon_length <- 710
#set primers length [bp]
primers_length <- 25
#flip_flop flag is a variable that determines whether basecalling is performed with the standard
#algorithm or with the flip-flop 'large' algorithm (when set equal to 1), slower but more accurate 
flip_flop_flag <- 0
#doBlast flag (if equal to 1, blast consensus sequence against NCBI nt database)
doBlast <- 1
########################################################################################################
PIPELINE_DIR <- "/home/simone/MinION/MinION_scripts/mobile_laboratory_autosetup/ONTrack"
#MINICONDA DIR
MINICONDA_DIR <- "/home/simone/miniconda3"
#basecaller_dir (v2.3.7+e041753)
BASECALLER_DIR <- "/home/simone/MinION/software/ont-guppy-cpu/bin/"
#NCBI nt database
NTDB <- "/home/db/NR_2018_06_01/nt"
########################################################################################################
#path to ONTrack.R
ONTrack <- paste0(PIPELINE_DIR, "/ONTrack.R")
#num iterations of ONTrack pipeline (must be odd)
num_iterations <- 3
#path to DecONT.sh
DECONT <- paste0(PIPELINE_DIR, "/DecONT.sh")
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
