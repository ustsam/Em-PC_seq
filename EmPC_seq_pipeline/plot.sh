##Author: Ilona Christy Unarta 02/09/2020 (2nd Sept 2020)
##Requirement:
## 1. samtools version 1.9
## 2. python 2
## 3. pysam version
## 4. Numpy
## For Plotting one will need the following python libraries
## 5. Matplotlib


WORKDIR=$1 # Full Path to working directory
REFFILE=$2 # Full Path to reference fasta file
SCRIPTDIR=$3 # Full Path to directory containing all the scripts (python and C scripts)
#AMB=$4 # Ambiguity Threshold, ambiguity is the number of ways a transcript can be mapped to a reference genome
#MAX_DEPTH=$5 #Maximum Depth per sites for mpileup
minBQ=$4
minMQ=$5
PLOT=$6

if [ $PLOT == 1 ]
then
python ${SCRIPTDIR}/plotting.py ${WORKDIR}/data_sorted_filtered_MQ${minMQ}_BQ0.pileup_pysam_count ${WORKDIR}/data_sorted_filtered_MQ${minMQ}_BQ${minBQ}.pileup_pysam_count_filtered ${WORKDIR}/muta_reads_MQ${minMQ}_BQ${minBQ}_PosReads.txt ${WORKDIR} ${REFFILE}

elif [ ${PLOT} == 0 ]
then
python ${SCRIPTDIR}/plotting_preprocess.py ${WORKDIR}/data_sorted_filtered_MQ${minMQ}_BQ0.pileup_pysam_count ${WORKDIR}/data_sorted_filtered_MQ${minMQ}_BQ${minBQ}.pileup_pysam_count_filtered ${WORKDIR}/muta_reads_MQ${minMQ}_BQ${minBQ}_PosReads.txt ${WORKDIR} ${REFFILE}

fi
