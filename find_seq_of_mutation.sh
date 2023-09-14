#!/bin/bash
# path to the genome reference:
genome="/private/dropbox/Genomes/Human/hg38/hg38.fa"
chr="$1"
mut_position="$2"
start_seq=$((mut_position-$3))
end_seq=$((mut_position+$3))
# path to create temp file of BED format for query
out_temp="$PWD/TEMP.bed"
# create the temp file with the BED entry
printf "%s\t%s\t%s\n" "$chr" "$start_seq" "$end_seq" > "$out_temp"
# extract the sequence and print it
total_seq=$(bedtools getfasta -fi "$genome" -bed "$out_temp" | grep -A 1 "^>" | tail -n 1 | tr -d "\n" | awk '{print toupper($0)}')
echo "$total_seq"