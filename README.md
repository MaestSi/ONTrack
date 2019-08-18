
# ONTrack

**ONTrack** is a rapid and accurate MinION-based barcoding pipeline for tracking species biodiversity on site; starting from MinION sequence reads, the ONTrack pipeline is able to provide accurate consensus sequences in ~15 minutes per sample on a standard laptop. Moreover, a preprocessing pipeline is provided, so to make the whole bioinformatic analysis from raw fast5 files to consensus sequences straightforward and simple.

<p align="center">
  <img src="Figures/ONTrack_logo.png" alt="drawing" width="550" title="ONTrack_logo">
</p>

## Getting started

**Prerequisites**

* Miniconda3.
Tested with conda 4.6.11.
```which conda``` should return the path to the executable.
If you don't have Miniconda3 installed, you could download and install it with:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod 755 Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```

Then, after completing _ONTrack_ installation, set the _MINICONDA_DIR_ variable in **config_MinION_mobile_lab.R** to the full path to miniconda3 directory.

* Guppy, the software for basecalling and demultiplexing provided by ONT. Tested with Guppy v2.3 (ONTrack-v1.0), Guppy v3.0 (ONTrack-v1.1) and Guppy v3.1.
If you don't have [Guppy](https://community.nanoporetech.com/downloads) installed, choose an appropriate version and install it.
For example, you could download and unpack the archive with:
```
wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_version_of_interest.tar.gz
tar -xf ont-guppy-cpu_version_of_interest.tar.gz
```
A directory _ont-guppy-cpu_ should have been created in your current directory.
Then, after completing _ONTrack_ installation, set the _BASECALLER_DIR_ variable in **config_MinION_mobile_lab.R** to the full path to _ont-guppy-cpu/bin_ directory.

* NCBI nt database (optional, in case you want to perform a local Blast analysis of your consensus sequences).

For downloading the database (~65 GB):

```
mkdir NCBI_nt_db
cd NCBI_nt_db
echo `date +%Y-%m-%d` > download_date.txt
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*
targz_files=$(find . | grep \.tar\.gz$ | sed 's/\.\///g')
for f in $targz_files; do tar -xzvf $f; done
rm $targz_files
```

Then, after completing the _ONTrack_ installation, set the _NTDB_ variable in **config_MinION_mobile_lab.R** to the full path to NCBI_nt_db/nt

**Installation**

```
git clone https://github.com/MaestSi/ONTrack.git
cd ONTrack
chmod 755 *
./install.sh
```

A conda environment named _ONTrack_env_ is created, where blast, emboss, vsearch, seqtk, mafft, porechop, minimap2, samtools, nanopolish, bedtools, pycoQC and R with package Biostrings are installed.
Then, you can open the **config_MinION_mobile_lab.R** file with a text editor and set the variables _PIPELINE_DIR_ and _MINICONDA_DIR_ to the value suggested by the installation step.

## Overview

<p align="center">
  <img src="Figures/ONTrack_pipeline_flowchart.png" alt="drawing" width="700" title="ONTrack_pipeline_flowchart">
</p>

## Usage

The ONTrack pipeline can be applied either starting from raw fast5 files, or from already basecalled and demultiplexed sequences.
In both cases, the first step of the pipeline requires you to open the **config_MinION_mobile_lab.R** file with a text editor and to modify it according to the features of your sequencing experiment and your preferences.
If you have already basecalled and demultiplexed your sequences, you can run the pipeline using the **ONTrack.R** script.
Otherwise, you can run the pipeline using the **Launch_MinION_mobile_lab.sh** script.

**ONTrack.R**

Usage: Rscript ONTrack.R \<home_dir\> \<fast5_dir\> \<sequencing_summary.txt\>

Note: Activate the virtual environment with ```source activate ONTrack_env``` before running. The script is run by **MinION_mobile_lab.R**, but can be also run as a main script if you have already basecalled and demultiplexed your sequences. If less than 200 reads are available after contaminants removal, a warning message is printed out, but still a consensus sequence is produced.

Inputs:
* \<home_dir\>: directory containing fastq and fasta files for each sample
* \<fast5_dir\>: directory containing raw fast5 files for nanopolish polishing, optional
* \<sequencing_summary.txt\>: sequencing summary file generated during base-calling, used to speed-up polishing, optional

Outputs (saved in <home_dir>):
* \<"sample_name".contigs.fasta\>: polished consensus sequence in fasta format
* \<"sample_name".blastn.txt\>: blast analysis of consensus sequence against NCBI nt database (if _do_blast_flag_ variable is set to 1 in **config_MinION_mobile_lab.R**)
* \<"sample_name"\>: directory including intermediate files

**Launch_MinION_mobile_lab.sh**

Usage:
Launch_MinION_mobile_lab.sh \<fast5_dir\>

Note: modify **config_MinION_mobile_lab.R** before running; the script runs the full pipeline from raw fast5 files to consensus sequences.

Input
* \<fast5_dir\>: directory containing raw fast5 files

Outputs (saved in \<fast5_dir\>_analysis/analysis):
* \<"sample_name".contigs.fasta\>: polished consensus sequence in fasta format
* \<"sample_name".blastn.txt\>: blast analysis of consensus sequence against NCBI nt database (if _do_blast_flag_ variable is set to 1 in **config_MinION_mobile_lab.R**)
* \<"sample_name"\>: directory including intermediate files

Outputs (saved in \<fast5_dir\>_analysis/qc):
* Read length distributions and pycoQC report

Outputs (saved in \<fast5_dir\>_analysis/basecalling):
* Temporary files for basecalling

Outputs (saved in \<fast5_dir\>_analysis/preprocessing):
* Temporary files for demultiplexing, filtering based on read length and adapters trimming

## Auxiliary scripts

In the following, auxiliary scripts run either by **ONTrack.R** or by **Launch_MinION_mobile_lab.sh** are listed. These scripts should not be called directly.

**MinION_mobile_lab.R**

Note: script run by _Launch_MinION_mobile_lab.sh_.

**config_MinION_mobile_lab.R**

Note: configuration script, must be modified before running _Launch_MinION_mobile_lab.sh_ or _ONTrack.R_.

**subsample_fast5.sh**

Note: script run by _MinION_mobile_lab.R_ if _do_subsampling_flag_ variable is set to 1 in _config_MinION_mobile_lab.R_.

**remove_long_short.pl**

Note: script run by _MinION_mobile_lab.R_ for removing reads shorter than mean - 2\*sd and longer than mean + 2\*sd.

**decONT.sh**

Note: script run by _ONTrack.R_ for clustering reads at 70% identity and keeping only reads in the most abundant cluster, if _do_clustering_flag_ variable is set to 1 in _config_MinION_mobile_lab.R_.

## Checking scripts

**Sanger_check.sh**

Usage: Sanger_check.sh \<consensus dir\> \<sanger dir\>

Note: set _BLASTN_ variable to blastn executable inside the script; sample name should contain the sample id (e.g. BC01)

Inputs:
* \<consensus dir\>: directory containing files "sample_name".contigs.fasta obtained with the _ONTrack_ pipeline
* \<sanger dir\>: directory containing fasta files reference_"sample_name".fasta obtained with Sanger sequencing

Output (saved in \<contigs dir\>):
* <results_"sample_name".txt>: file including alignment of MinION consensus sequence to corresponding Sanger sequence
* \<Sanger_check_report.txt\>: file including overall alignment statistics and number of uncertain nucleotides in Sanger sequences

**Calculate_mapping_rate.sh**

Usage: Calculate_mapping_rate.sh \<reads\> \<draft reads\> \<consensus sequence\>

Note: set _MINIMAP2_ and _SAMTOOLS_ variables to minimap2 and samtools executables inside the script

Inputs:
* \<reads\>: MinION reads in fastq or fasta format
* \<draft reads\>: MinION reads in fastq or fasta format used for creating draft consensus sequence, after contaminants removal
* \<consensus sequence\>: polished consensus sequence in fasta format

Output (saved in current directory):
* \<"sample_name"_report_mapping_rate.txt\>: mapping rate statistics

**Calculate_error_rate.sh**

Usage: Calculate_error_rate.sh \<reads\> \<reference\>

Note: set _MINIMAP2_ and _SAMTOOLS_ variables to minimap2 and samtools executables inside the script

Inputs:
* \<reads\>:  MinION reads in fastq or fasta format
* \<reference\>: Sanger sequence corresponding to MinION reads

Outputs:
* \<"sample_name"_error_rate_stats.txt\>: error rate statistics

## Contaminants inspection analysis

When the mapping rate of all reads from a sample is not in the range 95%-100%, you might be interested either in spotting if there is a predominant contaminant, or in trying to rescue the consensus sequence of your sample, if based on Blast analysis you realize that the consensus sequence from the most abundant cluster is not from the sample that you were supposed to sequence.
In these cases, you could try to retrieve the reads that don't map to your consensus sequence, and run the _ONTrack_ pipeline again just on those reads. Remember in these cases to set the _do_clustering_flag_ variable to 1 in the _config_MinION_mobile_lab.R_ file.
As an example, you could use the following code to retrieve unmapped reads for sample BC01 and save them to _contaminants_analysis_ folder.

```
SAMPLE_NAME=BC01
ANALYSIS_DIR=/path/to/fast5_reads_analysis/analysis
PIPELINE_DIR=/path/to/ONTrack

source activate ONTrack_env

cd $ANALYSIS_DIR
$PIPELINE_DIR"/Calculate_mapping_rate.sh" $SAMPLE_NAME".fastq" $SAMPLE_NAME"/"$SAMPLE_NAME"_decont.fastq" $SAMPLE_NAME".contigs.fasta"

if [ ! -d "contaminants_analysis" ]; then
  mkdir contaminants_analysis
fi

samtools view -f4 -b $SAMPLE_NAME"_reads_on_contig.bam" > "contaminants_analysis/"$SAMPLE_NAME"_unmapped.bam"
bedtools bamtofastq -i "contaminants_analysis/"$SAMPLE_NAME"_unmapped.bam" -fq "contaminants_analysis/"$SAMPLE_NAME".fastq"
seqtk seq -A "contaminants_analysis/"$SAMPLE_NAME".fastq" > "contaminants_analysis/"$SAMPLE_NAME".fasta"
```

## Meta-barcoding analysis (experimental)

Although the ONTrack pipeline is not intended for analysing meta-barcoding samples, you might be interested in sorting out sequences coming from different species and running the ONTrack pipeline on the most abundant species separately.
The **MetatONTrack.sh** script reproduces what the EPI2ME 16S workflow does, blasting each read against an NCBI-downloaded database (e.g. 16S Bacterial), and afterwards saving sets of reads matching the different species to separate files. You can then run the **ONTrack.R** script on them, for obtaining a more accurate consensus sequence (set _do_clustering_flag_ to 0 in _config_MinION_mobile_lab.R_). This feature is experimental, and has only been tested on a pool of 7 samples with 80% maximum sequence identity based on pairwise alignment of Sanger sequences.

**MetatONTrack.sh**

Usage: MetatONTrack.sh \<fastq reads\> \<min num reads\>

Note: set _BLASTN_, _SEQTK_ and _DB_ variables to blastn, seqtk executables and to an NCBI Blast-indexed database respectively inside the script

Inputs:
* \<fastq reads\>:  MinION fastq reads from a meta-barcoding experiment
* \<min num reads\>: minumum number of reads supporting the identification of a species 

Outputs:
* \<MetatONTrack_output\>: directory containing fastq and fasta files for running the **ONTrack.R** script
* \<MetatONTrack_output_logs\>: directory containing txt files storing read IDs corresponding to each species, a "sample_name"\_Blast_species_counts.txt file storing the number of reads supporting each species and some other temporary files

## Citation

If this tool is useful for your work, please consider citing our [manuscript](https://www.mdpi.com/2073-4425/10/6/468).

Maestri S, Cosentino E, Paterno M, Freitag H, Garces JM, Marcolungo L, Alfano M, Njunjić I, Schilthuizen M, Slik F, Menegon M, Rossato M, Delledonne M. A Rapid and Accurate MinION-Based Workflow for Tracking Species Biodiversity in the Field. Genes. 2019; 10(6):468.

## Side notes

As a real-life _Pokédex_, the workflow described in our [manuscript](https://www.mdpi.com/2073-4425/10/6/468) will facilitate tracking biodiversity in remote and biodiversity-rich areas. For instance, during a [Taxon Expedition](https://taxonexpeditions.com/) to Borneo, our analysis confirmed the novelty of a [beetle](https://www.theguardian.com/science/2018/apr/30/new-beetle-species-named-after-leonardo-dicaprio) species named after Leonardo DiCaprio.

Last but not least, special thanks to Davide Canevazzi ([davidecanevazzi](https://github.com/davidecanevazzi)) and Luca Marcolungo ([Liukvr](https://github.com/Liukvr)) for helping me out in setting up and debugging the pipeline.
