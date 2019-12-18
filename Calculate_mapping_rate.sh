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

READS=$1
DRAFT_READS=$2
CONTIG=$3

MINIMAP2=minimap2 #specify full path if you want to use a version of the program that is not in your PATH
SAMTOOLS=samtools #specify full path if you want to use a version of the program that is not in your PATH

wdir=$(realpath $(pwd))
reads_full=$wdir"/"$READS
draft_reads_full=$wdir"/"$DRAFT_READS
contig_full=$wdir"/"$CONTIG
sample_id=$(echo $READS | sed 's/\.fast.*//')

echo $sample_id > $wdir"/"$sample_id"_report_mapping_rate.txt"
echo "Mapping command: " >> $wdir"/"$sample_id"_report_mapping_rate.txt"
echo "minimap2 -ax map-ont "$contig_full $draft_reads_full" | samtools view -Sbh -F2048 -F256 | samtools sort - -o "$wdir"/"$sample_id"_draft_reads_on_contig.bam" >> $wdir"/"$sample_id"_report_mapping_rate.txt"
$MINIMAP2 -ax map-ont $contig_full $draft_reads_full | $SAMTOOLS view -Sbh -F2048 -F256 | $SAMTOOLS sort - -o $wdir"/"$sample_id"_draft_reads_on_contig.bam"
echo " " >> $wdir"/"$sample_id"_report_mapping_rate.txt"
echo "Flagstat command: " >> $wdir"/"$sample_id"_report_mapping_rate.txt"
echo "samtools flagstat $wdir"/""$sample_id"_draft_reads_on_contig.bam"  >> $wdir"/"$sample_id"_report_mapping_rate.txt"
$SAMTOOLS flagstat $wdir"/"$sample_id"_draft_reads_on_contig.bam"  >> $wdir"/"$sample_id"_report_mapping_rate.txt"
echo " " >> $wdir"/"$sample_id"_report_mapping_rate.txt"
echo "Mapping command: " >> $wdir"/"$sample_id"_report_mapping_rate.txt"
echo "minimap2 -ax map-ont "$contig_full $reads_full" | samtools view -Sbh -F2048 -F256 | samtools sort - -o "$wdir"/"$sample_id"_reads_on_contig.bam" >> $wdir"/"$sample_id"_report_mapping_rate.txt"
$MINIMAP2 -ax map-ont $contig_full $reads_full | $SAMTOOLS view -Sbh -F2048 -F256 | $SAMTOOLS sort - -o $wdir"/"$sample_id"_reads_on_contig.bam"
echo " " >> $wdir"/"$sample_id"_report_mapping_rate.txt"
echo "Flagstat command: " >> $wdir"/"$sample_id"_report_mapping_rate.txt"
echo "samtools flagstat "$wdir"/"$sample_id"_reads_on_contig.bam" >> $wdir"/"$sample_id"_report_mapping_rate.txt"
$SAMTOOLS flagstat $wdir"/"$sample_id"_reads_on_contig.bam" >> $wdir"/"$sample_id"_report_mapping_rate.txt"
