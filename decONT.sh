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

FASTA=$1
VSEARCH=$2
SEQTK=$3

sample_id=$(echo $(basename $FASTA) | sed 's/\.fasta//')
wdir=$(realpath $(dirname $FASTA))
prefix=$wdir"/"$sample_id"_cluster_"
prefix_bn=$(basename $prefix)
$VSEARCH --cluster_smallmem $FASTA --usersort --id 0.7 --iddef 2 --clusterout_sort --fasta_width 0 --strand both --sizeout --consout $wdir"/consensus_"$sample_id".fasta" --clusters $prefix
centroid_mac=$(head -n1 $wdir"/consensus_"$sample_id".fasta")
n_cl=$(ls $wdir | grep $prefix_bn | wc -l)
id_centroid_tmp=$(echo $centroid_mac | sed 's/centroid=//g')
id_centroid=$(echo $id_centroid_tmp | sed 's/;seqs.*$//g')

for (( i=0; i <= $n_cl -1 ; i++ ))
do
  cat $prefix$i | grep $id_centroid && mac_file=$prefix$i && break
done

cat $mac_file | grep "^>" | sed 's/^>//' | sed 's/\;.*$//' > $wdir"/"$sample_id"_ids_mac.txt"

$SEQTK subseq $wdir"/"$sample_id".fasta" $wdir"/"$sample_id"_ids_mac.txt" > $wdir"/"$sample_id"_decont.fasta"
$SEQTK subseq $(realpath $wdir"/"$sample_id".fastq") $wdir"/"$sample_id"_ids_mac.txt" > $wdir"/"$sample_id"_decont.fastq"

tmp_dir=$wdir"/"$sample_id"/decontam_tmp_"$sample_id
mkdir $tmp_dir

mv $wdir"/"$sample_id"_ids_mac.txt" $tmp_dir
mv $wdir"/consensus_"$sample_id".fasta" $tmp_dir

for cluster in $(find $wdir -maxdepth 1 -mindepth 1 | grep $prefix_bn); do
  mv $cluster $tmp_dir;
done
