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
REFERENCE=$2

MINIMAP2=/home/simone/miniconda3/bin/minimap2
SAMTOOLS=/home/simone/miniconda3/bin/samtools

SAMPLE_NAME=$(echo $(basename $READS) | sed 's/\.fast.//g')
WORKING_DIR=$(dirname $(realpath $READS))
BAM=$WORKING_DIR"/"$SAMPLE_NAME"_vs_reference.bam"

$MINIMAP2 -ax map-ont $REFERENCE $READS | $SAMTOOLS view -hSb -F2304 | $SAMTOOLS sort - -o $BAM && $SAMTOOLS index $BAM
$SAMTOOLS stats $BAM | grep "^SN" | cut -f 2- > $WORKING_DIR"/"$SAMPLE_NAME"_error_rate_stats.txt"
