#!/bin/bash
#SBATCH --job-name=bbmap
#SBATCH --account=project_2001175
#SBATCH --time=0-04:00:00
#SBATCH --partition=small
#SBATCH --mem=4GB
#SBATCH --ntasks=4

REF="ref/HAMBI-24.16S-amplicon-mapdb.ffn"

cat $1 | while read MYLIB;
do

mkdir mapqc/${MYLIB}

# map to HAMBI 16S sequences
bbmap.sh ow=t int=f maxindel=20 minid=0.85 ambiguous=best \
ref=${REF} \
ehist=mapqc/${MYLIB}/${MYLIB}.ehist \
qahist=mapqc/${MYLIB}/${MYLIB}.qahist \
mhist=mapqc/${MYLIB}/${MYLIB}.mhist \
idhist=mapqc/${MYLIB}/${MYLIB}.idhist \
scafstats=mapqc/${MYLIB}/${MYLIB}.scafstats \
statsfile=mapqc/${MYLIB}/${MYLIB}.mapstats \
in=processedreads/${MYLIB}.06-maxee.fastq.gz \
out=mappedreads/${MYLIB}.sam

samtools stats mappedreads/${MYLIB}.sam > mapqc/${MYLIB}/${MYLIB}.samtoolsstats

# get number of reads mapped to each 16S sequence
pileup.sh ow=t in=mappedreads/${MYLIB}.sam \
out=mappedreads/${MYLIB}.coverage \
rpkm=mappedreads/${MYLIB}.rpkm

done