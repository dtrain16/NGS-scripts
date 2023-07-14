#!/bin/bash

set -u

## get average peak length from macs2 peak calling output (bed file)

if [ "$#" -lt 1 ]; then
        echo "Missing required arguments!"
        echo "USAGE: get_peak_length.sh <sample>"
        echo "EXAMPLE: get_peak_length.sh col0-r1.merged.bed"
        exit 1
fi

sample=$1

awk '{ $10= $3 - $2} { sum += $10} END {print sum / NR}' $sample


 
