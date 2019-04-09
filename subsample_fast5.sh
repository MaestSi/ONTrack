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

INPUT_DIR=$1
NUM_K_READS=$2

output_dir=$(dirname $INPUT_DIR)/$(basename $INPUT_DIR)"_"$NUM_K_READS"k_subsampled_reads"
mkdir $output_dir

for f in $(find $INPUT_DIR | grep "\\.fast5" | shuf -n $((1000*$NUM_K_READS)));
  do  cp $f $output_dir;
done
