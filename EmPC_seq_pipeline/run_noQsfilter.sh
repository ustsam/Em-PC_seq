#!/bin/bash

#$ -l h_rt=336:0:0
#$ -cwd
#$ -l arch=linux-x64
#$ -S /bin/bash
#$ -N cirseq

set -e

WORKDIR=$1
REFFILE=$2
SCRIPTDIR=$3
MAPPER=$4
REPEATS=$5	#JHL, 20160305
READLENGTH=$6	#JHL, 20160305

if [ $MAPPER == 'bwa' ]; then

	FASTQS="${@:7}"
elif [ $MAPPER == 'STAR' ]; then

	FASTQS="${@:8}"
	GFF=$7
else

	echo "$MAPPER not detected, aborted"
	exit
fi


#save command line arguments, JHL, 20160305
echo "$WORKDIR $REFFILE $SCRIPTDIR $MAPPER $REPEATS $READLENGTH $FASTQS" > cmd

#clean up
mkdir -p ${WORKDIR}/index

# clean up reference file path, just get the filename
REFNAME=`echo $REFFILE | sed -e 's/.*\///g' | awk -F '.' '{print $1}'`
REFIDX=${WORKDIR}/index/${REFNAME}

# generate consensus
# JHL, 20160303, repeatcutoff (3) & readlength (300), to be put up to sys argv
python ${SCRIPTDIR}/ConsensusGeneration.py $WORKDIR $REPEATS $READLENGTH $FASTQS

## MAPPER CHEKCER

#check if mapper exist
which $MAPPER > /dev/null 2>&1
if [ $? -eq 0 ] ; then

        echo "$MAPPER is detected, relocalization will be done using $MAPPER"

else

        echo "$MAPPER not detected, aborted"
        exit

fi

if [ $MAPPER == 'bwa' ]; then
# build bwa index #CTH, 20181217
bwa index ${REFFILE}

# align consensus, change to bwa and add screening of supplementary and secondary alignment using samtools
# CTH, 20181217
bwa mem ${REFFILE} ${WORKDIR}/1_consensus.fastq.gz |samtools view -h -F 0x100|samtools view -F 0x800|grep -E -v "^@"|gzip -c > ${WORKDIR}/2_alignment.sam.gz

# preprocess 1, change bwa version preprocessing
python ${SCRIPTDIR}/preprocessing_1_bwa.py ${WORKDIR}

gunzip 3_alignment.sam.gz
python ${SCRIPTDIR}/relocalize.py 3_alignment.sam
bwa mem ${REFFILE} relocalized.fastq | awk '$2!=4 && $5>=30 && $6 !~ /D/ && $6 !~ /S/' | grep -v fastq|  gzip -c  > ${WORKDIR}/relocalized.sam.gz
python ${SCRIPTDIR}/check_ambiguity.py ${WORKDIR} relocalized.sam.gz 1
mv 11_alignment.sam.gz 3_alignment_relocalized.sam.gz
rm relocalized.sam.gz relocalized.fastq 
# align again, change to bwa and add screening of supplementary and secondary alignment using samtools
# CTH, 20181217
bwa mem ${REFFILE} ${WORKDIR}/4_rearranged.fastq.gz |samtools view -h -F 0x100|samtools view -F 0x800|grep -E -v "^@"|gzip -c > ${WORKDIR}/6_alignment.sam.gz

# preprocess 2, change to bwa and add screening of segment unmapped using samtools
# CTH, 20181217
python ${SCRIPTDIR}/preprocessing_2_bwa.py ${WORKDIR}

gunzip 7_alignment.sam.gz
python ${SCRIPTDIR}/relocalize.py 7_alignment.sam
bwa mem ${REFFILE} relocalized.fastq | awk '$2!=4 && $5>=30 && $6 !~ /D/ && $6 !~ /S/' | grep -v fastq | gzip -c  > ${WORKDIR}/relocalized.sam.gz 
python ${SCRIPTDIR}/check_ambiguity.py ${WORKDIR} relocalized.sam.gz 2
mv 11_alignment.sam.gz 7_alignment_relocalized.sam.gz
#rm relocalized.sam.gz relocalized.fastq

bwa mem ${REFFILE} ${WORKDIR}/5_rotated.fastq.gz |samtools view -h -F 0x4|grep -E -v "^@"|gzip -c > ${WORKDIR}/9_alignment.sam.gz

bwa mem ${REFFILE} ${WORKDIR}/8_rotated.fastq.gz |samtools view -h -F 0x4|grep -E -v "^@"|gzip -c > ${WORKDIR}/10_alignment.sam.gz

# preprocess 3, change to Alignment and edit distance screening.
# CTH, 20181217
python ${SCRIPTDIR}/preprocessing_3_bruteforce_bwa.py ${WORKDIR}

# combine results
cat ${WORKDIR}/3_alignment_relocalized.sam.gz ${WORKDIR}/7_alignment_relocalized.sam.gz ${WORKDIR}/11_alignment.sam.gz > ${WORKDIR}/data.sam.gz

# cleanup intermediate files
rm -f ${WORKDIR}/3_alignment.sam ${WORKDIR}/7_alignment.sam ${WORKDIR}/2_alignment.sam.gz ${WORKDIR}/4_rearranged.sam.gz ${WORKDIR}/5_rotated.fastq.gz ${WORKDIR}/6_alignment.sam.gz ${WORKDIR}/8_rotated.fastq.gz ${WORKDIR}/9_alignment.sam.gz ${WORKDIR}/10_alignment.sam.gz relocalized.sam.gz relocalized.fastq  

fi
