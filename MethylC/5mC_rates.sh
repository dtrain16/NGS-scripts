#!/bin/bash

set -u

# Get methylation rates for all contexts across Chr1-5 as well as % CHH in chloroplast genome as an indication of sodium bisulfite conversion efficiency (unconverted CHH in Cp and Mt genome)

if [ "$#" -lt 2 ]; then
	echo "Missing required arguments!"
	echo "USAGE: methylation_rates.sh <sample> <file>"
	echo "EXAMPLE: methylation_rates.sh col0-r1 bed/cov"
	exit 1
fi

sample=$1
file=$2

cg="${sample}_CG*.${file}"
chg="${sample}_CHG*.${file}"
chh="${sample}_CHH*.${file}"

echo "5mC % in $1"

echo "mCG Chr1-5: "$cg" "
grep -e "Chr1" -e "Chr2" -e "Chr3" -e "Chr4" -e "Chr5" $cg |  awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print ((met / total))}'

echo "mCHG Chr1-5: "$chg" "
grep -e "Chr1" -e "Chr2" -e "Chr3" -e "Chr4" -e "Chr5" $chg |  awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print ((met / total))}'

echo "mCHH Chr1-5: "$chh" "
grep -e "Chr1" -e "Chr2" -e "Chr3" -e "Chr4" -e "Chr5" $chh |  awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print ((met / total))}'

echo "mCHH ChrC: "$chh" "
grep -e "ChrC"  $chh | awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print ((met / total))}'

echo "mCHH ChrM: "$chh" "
grep -e "ChrM"  $chh | awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print ((met / total))}'

echo "DONE"
