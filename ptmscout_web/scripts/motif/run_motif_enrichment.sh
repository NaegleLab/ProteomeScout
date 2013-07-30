#!/bin/bash
set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd "$DIR"

foreground_file=$1
background_file=$2
output_file=$3

shift 3

perl c_stat_sig_peptide_refpass_depth_memory.phospho_acetyl.pl --foreground $foreground_file --background $background_file $@ | grep "NUM_TESTS_FOR_KRISTEN" > $output_file.log 2> /dev/null
perl process_motifs.pl log.$foreground_file.$background_file > $output_file 2> /dev/null
rm log.$foreground_file.$background_file
