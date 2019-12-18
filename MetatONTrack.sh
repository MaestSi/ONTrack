#!/bin/bash
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

if [ $# -eq 0 ]; then
  echo -e "\nUsage: $(basename $0) <Fastq_reads> <Min_num_reads_supporting_a_species> \n";
  exit;
fi

FASTQ_READS=$1
MIN_READS=$2
FASTQ_READS_FULL=$(realpath $FASTQ_READS)

######################################################################
THREADS=8
DB=/path/to/NCBI_Blast-indexed_database #e.g. PRJNA33175 BioProject for Bacterial 16S
BLASTN=blastn
SEQTK=seqtk
######################################################################

SAMPLE_NAME=$(echo $(basename $FASTQ_READS_FULL) | sed 's/\.fastq.*//')
WORKING_DIR=$(dirname $FASTQ_READS_FULL)
OUTPUT_DIR=$WORKING_DIR"/MetatONTrack_output"
LOGS_DIR=$WORKING_DIR"/MetatONTrack_output_logs"

$SEQTK seq -A $FASTQ_READS_FULL > $WORKING_DIR"/"$SAMPLE_NAME".fasta"
$BLASTN -db $DB -query $WORKING_DIR"/"$SAMPLE_NAME".fasta" -num_threads $THREADS -outfmt "6 qseqid evalue salltitles" -max_target_seqs 1 -perc_identity 0.77 -qcov_hsp_perc 0.3 \
> $WORKING_DIR"/"$SAMPLE_NAME"_top_Blast_hits_tmp.txt"
cat $WORKING_DIR"/"$SAMPLE_NAME"_top_Blast_hits_tmp.txt" | sort -u -k1,1 -s > $WORKING_DIR"/"$SAMPLE_NAME"_top_Blast_hits.txt"
rm $WORKING_DIR"/"$SAMPLE_NAME"_top_Blast_hits_tmp.txt"
cat $WORKING_DIR"/"$SAMPLE_NAME"_top_Blast_hits.txt" | cut -f3 | sort | uniq -c | sort -nr > $WORKING_DIR"/"$SAMPLE_NAME"_Blast_taxa_counts.txt"
cat $WORKING_DIR"/"$SAMPLE_NAME"_top_Blast_hits.txt" | cut -f3 | cut -d" " -f2,3 | sort | uniq -c | sort -nr > $WORKING_DIR"/"$SAMPLE_NAME"_Blast_species_counts.txt"
cat $WORKING_DIR"/"$SAMPLE_NAME"_top_Blast_hits.txt" | cut -f3 | cut -d" " -f2 | sort | uniq -c | sort -nr > $WORKING_DIR"/"$SAMPLE_NAME"_Blast_genera_counts.txt"
cat $WORKING_DIR"/"$SAMPLE_NAME"_Blast_species_counts.txt" | awk -v var="$MIN_READS" '{if ($1 > var) {print $2" "$3}}' > $WORKING_DIR"/"$SAMPLE_NAME"_detected_species.txt"

mkdir $OUTPUT_DIR
mkdir $LOGS_DIR

COUNTER=0

while IFS= read SPECIES_CURR
do
  let COUNTER=COUNTER+1
  if [ "$COUNTER" -lt 10 ]; then
    COUNTER_TWODIG="0"$COUNTER;
  else
    COUNTER_TWODIG=$COUNTER;
  fi
  SPECIES_CURR_MOD=$(echo $SPECIES_CURR | sed 's/ /_/g')
  cat $WORKING_DIR"/"$SAMPLE_NAME"_top_Blast_hits.txt" | grep "$SPECIES_CURR" | cut -f1  > $LOGS_DIR"/BC"$COUNTER_TWODIG"_"$SPECIES_CURR_MOD"_reads_IDs.txt"
  $SEQTK subseq $FASTQ_READS_FULL $LOGS_DIR"/BC"$COUNTER_TWODIG"_"$SPECIES_CURR_MOD"_reads_IDs.txt" > $OUTPUT_DIR"/BC"$COUNTER_TWODIG".fastq"
  $SEQTK seq -A $OUTPUT_DIR"/BC"$COUNTER_TWODIG".fastq" > $OUTPUT_DIR"/BC"$COUNTER_TWODIG".fasta"
done < $WORKING_DIR"/"$SAMPLE_NAME"_detected_species.txt"

mv $WORKING_DIR"/"$SAMPLE_NAME"_top_Blast_hits.txt" \
$WORKING_DIR"/"$SAMPLE_NAME"_Blast_taxa_counts.txt" \
$WORKING_DIR"/"$SAMPLE_NAME"_Blast_genera_counts.txt" \
$WORKING_DIR"/"$SAMPLE_NAME"_Blast_species_counts.txt" \
$WORKING_DIR"/"$SAMPLE_NAME"_detected_species.txt" \
$LOGS_DIR
