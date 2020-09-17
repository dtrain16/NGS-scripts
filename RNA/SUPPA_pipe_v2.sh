#!/bin/bash
set -u

# Performs event-based splicing analysis using SUPPA2 using stringtie TPM output 

if [ "$#" -lt 5 ]; then
echo "Missing arguments!"
echo "USAGE: SUPPA_pipe_v2.sh <annotation> <file dir> <group1> <group2> <name>"
echo "EXAMPLE: SUPPA_pipe_v2.sh AtRTD2_QUASI_19April2016.gtf raw_files/ col0_rep1,col0_rep2,col0_rep3 grp7_rep1,grp7_rep2,grp7_rep3 RTD2-quasi"
exit 1
fi

#### Parameters
## annotation file
I=$1
## events output name
N=$5
## extracted TPM (see stringtie_extract_tpm.r)
S=$2
# group 1 IDs
grp1=$3
# group 2 IDs
grp2=$4

## quantification
mkdir tpm_output

fls=$(dir $S)

for i in $fls; do 
mkdir tpm_output/${i%%_stringtie.tpm*};
cp $S/$i tpm_output/${i%%_stringtie.tpm*}/abundance.tpm; done

python3.5 ~/bin/SUPPA-2.3/multipleFieldSelection.py -i tpm_output/*/abundance.tpm -k 1 -f 3 -o iso_tpm.txt

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

### PSI per event
python3.5 ~/bin/SUPPA-2.3/suppa.py psiPerEvent -i $N -e iso_tpm.txt -o ${N%%.allevents*}_events

### Differential splicing with local events
## PSI and TPM per condition
Rscript $HOME/scripts/RNA/split_file.R ./iso_tpm.txt $grp1 $grp2 ${grp1%%_rep*}_iso.tpm ${grp2%%_rep*}_iso.tpm

Rscript $HOME/scripts/RNA/split_file.R ./${N%%.allevents*}_events.psi $grp1 $grp2 ${grp1%%_rep*}_events.psi ${grp2%%_rep*}_events.psi

## differential splicing analysis
python3.5 ~/bin/SUPPA-2.3/suppa.py diffSplice -m empirical -gc -i $N -p ${grp2%%_rep*}_events.psi ${grp1%%_rep*}_events.psi -e ${grp2%%_rep*}_iso.tpm ${grp1%%_rep*}_iso.tpm -o ${N%%.ioe}_${grp2%%_rep*}-${grp1%%_rep*}_diffSplice_events

## differential trascript usage
### PSI per isoform
python3.5 ~/bin/SUPPA-2.3/suppa.py psiPerIsoform -g $I -e iso_tpm.txt -o ${M%%.ioi}

### Split PSI between 2 conditions:
Rscript $HOME/scripts/RNA/split_file.R ./${M%%.ioi}_isoform.psi $grp1 $grp2 ${grp1%%_rep*}_iso.psi ${grp2%%_rep*}_iso.psi

### diffsplice
python3.5 ~/bin/SUPPA-2.3/suppa.py diffSplice -m empirical -gc -i $M -p  ${grp2%%_rep*}_iso.psi ${grp1%%_rep*}_iso.psi -e ${grp2%%_rep*}_iso.tpm ${grp1%%_rep*}_iso.tpm -o ${M%%.ioi}_${grp2%%_rep*}-${grp1%%_rep*}_diffSplice_iso

## collect output
mkdir suppa2_output
mv ${N%%.ioe}_${grp2%%_rep*}-${grp1%%_rep*}_diffSplice_events* -t suppa2_output
mv ${M%%.ioi}_${grp2%%_rep*}-${grp1%%_rep*}_diffSplice_iso* -t suppa2_output


