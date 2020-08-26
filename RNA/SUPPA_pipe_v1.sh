#!/bin/bash
set -eu

# Performs event-based splicing analysis using SUPPA2 
# https://github.com/comprna/SUPPA#command-and-subcommand-structure
# tutorial: https://github.com/comprna/SUPPA/wiki/SUPPA2-tutorial

#### Parameters
## annotation file
I="$HOME/ref_seqs/AtRTD2/AtRTD2_QUASI_19April2016.gtf"
# events output name
N="RTD2-quasi"
# locate processed data abundance data, organised into individual folder per sample (see kallisto_v1.sh)
S="$HOME/ws/sal1_AS/raw_files/"

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

## generate local AS events
python3.5 ~/bin/SUPPA-2.3/suppa.py generateEvents -i $I -o $N -f ioe -e SE SS MX RI FL

#Put all the ioe events in the same file:
awk '
    FNR==1 && NR!=1 { while (/^<header>/) getline; }
    1 {print}
' *.ioe > ${N}.allevents.ioe
N="${N}.allevents.ioe"
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
# split samples
Rscript $HOME/bin/SUPPA-2.3/scripts/split_file.R iso_tpm.txt col0_rep1,col0_rep2,col0_rep3 grp7_rep1,grp7_rep2,grp7_rep3 col0_iso.tpm grp7_iso.tpm -i

~/bin/SUPPA-2.3/scripts/split_file.R ${N%%.allevents*}_events.psi col0_rep1,col0_rep2,col0_rep3 grp7_rep1,grp7_rep2,grp7_rep3 col0_events.psi grp7_events.psi -e


