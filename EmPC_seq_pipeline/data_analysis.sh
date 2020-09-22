##Author: Ilona Christy Unarta 02/09/2020 (2nd Sept 2020)
##Requirement: 
## 1. samtools (version 1.9)
## 2. python (version 2.7)
## 3. pysam (version 0.15.0)
## 4. Numpy (version 1.11.0)
## For Plotting one will need the following python libraries
## 5. Matplotlib (version 2.2.2)


WORKDIR=$1 # Full Path to working directory
REFFILE=$2 # Full Path to reference fasta file 
SCRIPTDIR=$3 # Full Path to directory containing all the scripts (python and C scripts)
AMB=$4 # Ambiguity Threshold, ambiguity is the number of ways a transcript can be mapped to a reference genome
MAX_DEPTH=$5 #Maximum Depth per sites for mpileup 
minBQ=$6
minMQ=$7
NUMSIM=$8
##SAM FILE PROCESSING
samtools view -bT ${REFFILE} ${WORKDIR}/data.sam.gz | samtools sort > ${WORKDIR}/data_sorted.bam

samtools view ${WORKDIR}/data_sorted.bam | awk -F 'UN:i:' '{print $2}' | sort -n | uniq -c | awk '{print $2,$1}' >& ${WORKDIR}/Distribution_of_Ambiguity.txt 

samtools view ${WORKDIR}/data_sorted.bam | awk -F 'UN:i:' -v j=${AMB} '$2<=j' | awk '$6!~/D/ && $6!~/I/' | awk -v minMQ=$minMQ '$5>=minMQ' |  samtools view -bT ${REFFILE} >& ${WORKDIR}/data_sorted_filtered.bam

#this is specifically for ribosomal DNA, since most rDNA are transcribed in negative direction, thus we assume that the 3' end is in the beginning of the transcript. 
#To remove the biotin position, we add the start of the transcript by one (assuming rDNA and negative transcription). we also remove the length of each trancsript by one.
#To make the file containing the starting site and length for each transcript for input to simulation.py
samtools view ${WORKDIR}/data_sorted_filtered.bam | awk '{print $3,$4+1,length($10)-1}' >& ${WORKDIR}/startlengths.txt 

#use this command if direction of transcription is known
#samtools view ${WORKDIR}/data_sorted_filtered.bam | awk '{if($2==16) print $4+1,length($10)-1,$2 ; else if($2==0) print $4,length($10)-1,$2}' >& startlengths.txt 

samtools index ${WORKDIR}/data_sorted_filtered.bam
#Make two kindes of pileup file
#we need to use the pysam to remove the biotin site
#the pysam code is specific for rDNA, where we assume that all transcription is in the negative direction, thus the starting of the transcript is the 3' end 
python ${SCRIPTDIR}/pysam_make_pileup.py -b ${WORKDIR}/data_sorted_filtered.bam -f ${REFFILE} -d ${MAX_DEPTH} -q ${minMQ} -Q ${minBQ} -c 0 -w ${WORKDIR}
python ${SCRIPTDIR}/pysam_make_pileup.py -b ${WORKDIR}/data_sorted_filtered.bam -f ${REFFILE} -d ${MAX_DEPTH} -q ${minMQ} -Q 0 -c 0 -w ${WORKDIR}

#Filter the mutation to remove the mutation with error rate per site larger than 0.01
awk '$5/$4<=0.01 && $5/$4>0.0' ${WORKDIR}/data_sorted_filtered_MQ${minMQ}_BQ${minBQ}.pileup_pysam_count >${WORKDIR}/mutaToKeep.txt
temp=$(wc -l ${WORKDIR}/mutaToKeep.txt | awk '{print $1}')
declare -i filelength=$temp
for ((i=1; i<=$filelength; i+=1))  ; do chrom=$(sed -n ""$i"p" ${WORKDIR}/mutaToKeep.txt | awk '{print $1}' ) ; pos=$(sed -n ""$i"p" ${WORKDIR}/mutaToKeep.txt | awk '{print $2}') ; awk -v chrom=$chrom -v pos=$pos '$1==chrom && $2==pos' ${WORKDIR}/muta_reads_MQ${minMQ}_BQ${minBQ}.txt  ; done &> ${WORKDIR}/muta_reads_MQ${minMQ}_BQ${minBQ}_filtered.txt
awk '{if($5/$4>0.01) $5=0 ; print $0}' ${WORKDIR}/data_sorted_filtered_MQ${minMQ}_BQ${minBQ}.pileup_pysam_count | sed 's/\ /\t/g' >& ${WORKDIR}/data_sorted_filtered_MQ${minMQ}_BQ${minBQ}.pileup_pysam_count_filtered

#Count the number of mutations at each position in a transcript
awk '{print $6}' ${WORKDIR}/muta_reads_MQ${minMQ}_BQ${minBQ}_filtered.txt | sort -n | uniq -c | awk '{print $2-1,$1}' >& ${WORKDIR}/muta_reads_MQ${minMQ}_BQ${minBQ}_PosReads.txt

#Calculate the error rate for each mutation type in the RNA transcript
awk '{print $3">"$5}' ${WORKDIR}/muta_reads_MQ${minMQ}_BQ${minBQ}_filtered.txt | sort -n | uniq -c | awk '{print $2,$1}' | sed 's/>/\ /g' >& ${WORKDIR}/MutaCountType ; for i in "A" "U" "G" "C" ; do COV=$(awk -v k=$i '{if($1==k) print $2}' ${WORKDIR}/BaseTypeCount_MQ${minMQ}_BQ${minBQ}.txt) ; awk -v k=$i -v cov=$COV '{if($1==k) print $1">"$2,$3/cov}' ${WORKDIR}/MutaCountType ; done >& ${WORKDIR}/MutationTypeSpectrum.txt

#Generate simulation files
#For this example, since we use the most strict requirement for the ambiguity (ambiguity is 1). 
#We do not need to do circularization for the simulation data (without -r tag for the simulation.py)
COVERAGE=$(awk '{sum+=$4} END {print sum}' ${WORKDIR}/data_sorted_filtered_MQ${minMQ}_BQ0.pileup_pysam_count)
MUTA=$(wc -l ${WORKDIR}/muta_reads_MQ${minMQ}_BQ${minBQ}_filtered.txt | awk '{print $1}')

#Create simulation data 
echo ${COVERAGE},${MUTA}
if [ $AMB == 1 ]
then

python ${SCRIPTDIR}/simulation.py -f ${REFFILE} -c ${COVERAGE} -m ${MUTA} -s ${NUMSIM} -w ${WORKDIR}
mkdir ${WORKDIR}/Simulation_Data
mv ${WORKDIR}/sim_*fastq.gz ${WORKDIR}/Simulation_Data

#Calculate the p-values for all the genomic sites with mutation and the site in the transcript
python ${SCRIPTDIR}/binomial_distribution.py -f ${REFFILE} -w ${WORKDIR} -p ${WORKDIR}/data_sorted_filtered_MQ${minMQ}_BQ${minBQ}.pileup_pysam_count_filtered  -m ${WORKDIR}/muta_reads_MQ${minMQ}_BQ${minBQ}_PosReads.txt 
fi
##PLOTTING
##This command can be commented if prefer to plot with other plotting tools
#python plotting.py "data_sorted_filtered_MQ${minMQ}_BQ0.pileup_pysam_count" "data_sorted_filtered_MQ${minMQ}_BQ${minBQ}.pileup_pysam_count_filtered" "muta_reads_MQ${minMQ}_BQ${minBQ}_PosReads.txt" ${WORKDIR}

#Extra commands
#samtools mpileup -d ${MAX_DEPTH} -B -q 30 -Q 30 -f ${REFFILE} ${WORKDIR}/data_sorted_filtered.bam | awk '{mystr=$5 ; A=gsub(/[Aa]/,"",mystr); T=gsub(/[Tt]/,"",mystr) ; C=gsub(/[Cc]/,"",mystr) ; G=gsub(/[Gg]/,"",mystr) ; print $1"\t"$2"\t"$3"\t"$4"\t"A+T+G+C"\t"A"\t"T"\t"G"\t"C;}' > ${WORKDIR}/data_sorted_filtered.mpileup

#awk '{coverage+=$4;error+=$5} END {printf "total coverage: %d\ntotal errors: %d\nerror rate: %f\n", coverage,error,error/coverage}' ${WORKDIR}/data_sorted_filtered.mpileup >& ${WORKDIR}/Coverage_Error_summary.txt

