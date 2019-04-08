# ONTrack
A MinION-based pipeline for tracking species biodiversity

**ONTrack.R**

Usage: Rscript ONTrack.R \<home_dir\> \<fast5_dir\> \<sequencing_summary.txt\>

Note: script run by MinION_mobile_lab.R, but can be also run as a main script if you have already basecalled and demultiplexed your sequences; config_MinION_mobile_lab.R must be in the same directory of ONTrack.R

Inputs:
* \<home_dir\>: directory containing fastq and fasta files for each sample
* \<fast5_dir\>: directory containing raw fast5 files for nanopolish polishing, optional
* \<sequencing_summary.txt\>: sequencing summary file generated during base-calling, used to speed-up polishing, optional

Outputs (saved in <home_dir>):
* \<"sample_name".contigs.fasta\>: polished consensus sequence in fasta format
* \<"sample_name".blastn.txt\>: blast analysis of consensus sequence against NCBI nt database (if doBlast flag is set to 1 in config_MinION_mobile_lab.R)
* \<"sample_name"\>: directory including intermediate files

**launch_MinION_mobile_lab.sh**

Usage:
launch_MinION_mobile_lab.sh \<fast5_dir\>

Note: config_MinION_mobile_lab.R and MinION_mobile_lab.R must be in the same directory of launch_MinION_mobile_lab.sh; modify config_MinION_mobile_lab.R before running; the script runs the full pipeline from fast5 files to consensus sequences

Input
* \<fast5_dir\>: directory containing raw fast5 files

Outputs (saved in \<fast5_dir\>_analysis/analysis):
* \<"sample_name".contigs.fasta\>: polished consensus sequence in fasta format
* \<"sample_name".blastn.txt\>: blast analysis of consensus sequence against NCBI nt database (if do_Blast is set to 1 in config_MinION_mobile_lab.R)

**MinION_mobile_lab.R**

Note: script run by launch_MinION_mobile_lab.sh

**config_MinION_mobile_lab.R**

Note: configuration script, must be modified before running launch_MinION_mobile_lab.sh or ONTrack.R

**subsample_fast5.sh**

Note: script run by MinION_mobile_lab.R if do_subsampling_flag is set to 1 in config_MinION_mobile_lab.R

**remove_long_short.pl**

Note: script run by MinION_mobile_lab.R for removing reads shorter than mean - 2\*sd and longer than mean + 2\*sd

**DecONT.sh**

Usage: DecONT.sh \<reads.fasta\> \<VSEARCH\> \<SEQTK\>

Note: script run by ONTrack.R; also fastq reads must be present in the same directory

Inputs:
* \<reads.fasta\>: MinION reads in fasta format
* \<VSEARCH\>: path to VSEARCH executable
* \<SEQTK\>: path to SEQTK executable

Outputs:
* \<reads_decont.fasta\>: reads in most abundant cluster in fasta format
* \<reads_decont.fastq\>: reads in most abundant cluster in fastq format
