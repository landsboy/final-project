#!/bin/bash
genome="/home/alu/netlandes/MHCpan/bed6_of_genom_38.bed"
chr="$1"
mut_position="$2"
start_seq=$((mut_position-15))
end_seq=$((mut_position+15))
# path to create temp file of BED format for query
out_temp="$PWD/TEMP.bed"
# create the temp file with the BED entry  $genome
printf "%s\t%s\t%s\n" "$chr" "$start_seq" "$end_seq" > "$out_temp"
bedtools intersect -wa -a "$genome" -b "$out_temp" > intersects_output2.bed