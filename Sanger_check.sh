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

CONTIGS_DIR=$1
SANGER_DIR=$2

BLASTN=path/to/blastn

contigs_files=$(find $CONTIGS_DIR -maxdepth 1 | grep "\\.contigs\\.fasta")
sanger_files=$(find $SANGER_DIR -maxdepth 1 | grep "reference.*\\.fasta")

for cf in $contigs_files
do
  bc=$(basename $cf | awk -F . '{ print $1 }')
  sf=$(find $SANGER_DIR -maxdepth 1 | grep "reference_"$bc".*\\.fasta")
  $BLASTN -query $cf -subject $sf > $CONTIGS_DIR"/results_"$bc".txt"
done

results_files=$(realpath $(find $CONTIGS_DIR -maxdepth 1 | grep "results_BC"))
report_file=$CONTIGS_DIR"/Sanger_check_report.txt"
echo "***************************" > $report_file

for rf in $results_files
do
  echo $(basename $rf) >> $report_file
  cat $rf | grep "Identities" >> $report_file
  num_unc=$(cat $rf | grep "Sbjct" | sed 's/Sbjct//g' | grep -P "N|B|D|H|K|M|R|S|V|W|Y" -o | grep -P "N|B|D|H|K|M|R|S|V|W|Y" -c)
  echo "Number of uncertain nucleotides in Sanger sequence: $num_unc" >> $report_file
  echo "***************************"  >> $report_file
done


