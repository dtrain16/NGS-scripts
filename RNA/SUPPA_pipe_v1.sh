#!/bin/bash
set -u

# Performs event-based splicing analysis using SUPPA2 
# https://github.com/comprna/SUPPA#command-and-subcommand-structure
# tutorial: https://github.com/comprna/SUPPA/wiki/SUPPA2-tutorial

if [ "$#" -lt 5 ]; then
echo "Missing arguments!"
echo "USAGE: SUPPA_pipe_v1.sh <annotation> <file dir> <group1> <group2> <name>"
echo "EXAMPLE: SUPPA_pipe_v1.sh $HOME/ref_seqs/AtRTD2_QUASI_19April2016.gtf $HOME/ws/sal1_AS/raw_files/ col0_rep1,col0_rep2,col0_rep3 grp7_rep1,grp7_rep2,grp7_rep3 RTD2-quasi"
exit 1
fi

#### Parameters
## annotation file
I=$1
## events output name
N=$5
## kallisto quant files
S=$2
# group 1 IDs
grp1=$3
# group 2 IDs
grp2=$4

## quantification
mkdir kallisto_output

fls=$(dir $S)
for i in $fls; do 
mkdir kallisto_output/${i%%_kallisto*};
cp $S/${i}/*/abundance.tsv kallisto_output/${i%%_kallisto*}/abundance.tsv; done

python3.5 ~/bin/SUPPA-2.3/multipleFieldSelection.py -i kallisto_output/*/abundance.tsv -k 1 -f 5 -o iso_tpm.txt

### generateEvents

mkdir generateEvents
cd generateEvents

## generate transcript events
python3.5 ~/bin/SUPPA-2.3/suppa.py generateEvents -i $I -o $N -f ioi
M="${N}.ioi"

## generate local AS events
python3.5 ~/bin/SUPPA-2.3/suppa.py generateEvents -i $I -o $N -f ioe -e SE SS MX RI FL

#Put all the ioe events in the same file:
awk '
    FNR==1 && NR!=1 { while (/^<header>/) getline; }
    1 {print}
' *.ioe > ${N}.allevents.ioe
N="${N}.allevents.ioe"

mv $M ../
mv $N ../

awk '
    FNR==1 && NR!=1 { while (/^<header>/) getline; }
    1 {print}
' *.gtf > ${N%%.allevents*}.allevents.gtf
mv *.allevents.gtf ../

cd ../

### PSI per transcript Isoform
python3.5 ~/bin/SUPPA-2.3/suppa.py psiPerEvent -i $N -e iso_tpm.txt -o ${N%%.allevents*}_events

### Differential splicing with local events
## PSI and TPM per condition
Rscript $HOME/scripts/RNA/split_file.R ./iso_tpm.txt $grp1 $grp2 ${grp1%%_rep*}_iso.tpm ${grp2%%_rep*}_iso.tpm

Rscript $HOME/scripts/RNA/split_file.R ./${N%%.allevents*}_events.psi $grp1 $grp2 ${grp1%%_rep*}_events.psi ${grp2%%_rep*}_events.psi

## differential splicing analysis
python3.5 ~/bin/SUPPA-2.3/suppa.py diffSplice -m empirical -gc -i $N -p ${grp2%%_rep*}_events.psi ${grp1%%_rep*}_events.psi -e ${grp2%%_rep*}_iso.tpm ${grp1%%_rep*}_iso.tpm -o ${N%%.ioe}_${grp2%%_rep*}-${grp1%%_rep*}_diffSplice

## .dpsi gives the DeltaPSI as the difference of the mean PSI between conditions, and the p_value of this difference.
## .psivec gives psi per sample

## differential trascript usage

### compute psi of transcript isoforms
python3.5 ~/bin/SUPPA-2.3/suppa.py psiPerIsoform -g $I -e ./iso_tpm_formatted.txt -o ./iso

### Split PSI between 2 conditions:
Rscript ~$HOME/scripts/RNA/split_file.R ./iso_isoform.psi $grp1 $grp2 ${grp1%%_rep*}_iso.psi ${grp2%%_rep*}_iso.psi

### diffsplice
python3.5 ~/bin/SUPPA-2.3/suppa.py diffSplice -m empirical -gc -i $M -p  ${grp2%%_rep*}_iso.psi ${grp1%%_rep*}_iso.psi -e ${grp2%%_rep*}_iso.tpm ${grp1%%_rep*}_iso.tpm -o ${M%%.ioi}_${grp2%%_rep*}-${grp1%%_rep*}_diffSplice_iso


