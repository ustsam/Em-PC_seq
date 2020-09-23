#Author: Ilona Christy Unarta (e-mail: icunarta@connect.ust.hk)
#Institution: HKUST (September 2020)
#Objective: To create pileup file with the biotin site neglected.


import numpy as np
import os 
import pysam
import optparse
import sys

p = optparse.OptionParser()
p.add_option('--bamfile', '-b', default="data_sorted_filtered.bam", help="input bam file")
p.add_option('--reffasta', '-f', default="rDNA1.fasta",  help="reference fasta file")
p.add_option('--max_depth', '-d', default="100000",  help="")
p.add_option('--min_mapq', '-q', default="30") 
p.add_option('--min_baseq','-Q',default="30")
p.add_option('--chop','-c',default="0")
p.add_option('--workdir','-w',default="./")

options, arguments = p.parse_args()
filename=str(options.bamfile)
ref_filename=str(options.reffasta)
Max_Depth=int(options.max_depth)
min_mapping_qual=int(options.min_mapq)
min_base_qual=int(options.min_baseq)
workdir=str(options.workdir)
chop=int(options.chop)

samfile = pysam.AlignmentFile(workdir+filename, "rb" )
m=open(ref_filename)

DNA_complement_dict={"A":"T","T":"A","C":"G","G":"C"}
RNA_complement_dict={"A":"U","U":"A","C":"G","G":"C"}
DNAtoRNA={"A":"A","T":"U","G":"G","C":"C"}

#Read reference fasta file and turn it into dictionary. we choose to read line by line to avoid running out of memory for large reference fasta file.
#dictionary keys are the chromosome name
#dictionary values are the fasta sequence of the corresponding chromosome
ref_dict={}
for line in m :
    temp=line.strip().split()
    if temp[0][0]==">":
        chrs=temp[0][1:]
        ref_dict[temp[0][1:]]=""
    else:
        ref_dict[chrs]+=temp[0]
m.close()

BaseCount_dict={"A":0,"U":0,"G":0,"C":0,"N":0}

pout=open(workdir+filename[:-4]+"_MQ"+str(min_mapping_qual)+"_BQ"+str(min_base_qual)+".pileup_pysam_count","wt")
m=open(workdir+"muta_reads_MQ"+str(min_mapping_qual)+"_BQ"+str(min_base_qual)+".txt","wt")
ref=pysam.FastaFile(ref_filename)

#Process the sam file per chromosome
for i in ref_dict.keys() :

    #compute_baq is false to avoid recalculation of the base quality score and preserve the errors at the end of the transcript which is especially important for nascent RNA transcripts.
    for pileupcolumn in samfile.pileup(i,compute_baq=False,max_depth=Max_Depth,min_mapping_quality=min_mapping_qual,fastafile=ref,min_base_quality=min_base_qual):
        refBase=ref_dict[i][pileupcolumn.pos]
        count_pos={"A":0,"T":0,"G":0,"C":0,"N":0}
        count=0
        mutaCount=0
        for pileupread in pileupcolumn.pileups :
            if not pileupread.is_del  and not pileupread.is_refskip and pileupread.indel==0 : 
                 direction=int(pileupread.alignment.flag)
                 #Ideally the direction of transcription for these nascent RNA is preserved, so that the flag in sam file would really indicate the direction of transcription
                 #For the current implementation of emPC-seq where the direction of the RNA is unknown, we assume that all transcripts are negatively transcribed for nascent rRNA, thus direction==16
                 direction=16 
                 if direction==0 :
                      readBase=pileupread.alignment.query_sequence[pileupread.query_position]
                      init=pileupread.alignment.reference_start
                      if pileupread.query_position >=chop and pileupread.alignment.qlen - pileupread.query_position >chop+1 : # plus one is to remove the biotin position
                           count+=1
                           count_pos[readBase]+=1
                           BaseCount_dict[DNAtoRNA[refBase]]+=1
                           if refBase != readBase:
                               mutaCount+=1
                               pos_from_3pr=pileupread.alignment.qlen-pileupread.query_position-1                               
                               
                               m.write(i+"\t"+str(pileupcolumn.pos+1)+"\t"+refBase+"\t"+str(pileupread.alignment.query_name)+"\t"+readBase+"\t"+str(pos_from_3pr)+"\t"+str(pileupread.alignment.qlen)+"\t"+str(direction)+"\t"+str(pileupread.alignment.query_qualities[pileupread.query_position])+"\t"+str(pileupread.alignment.mapping_quality)+"\t"+str(init)+"\n")

                 elif direction==16 :
                      readBase=pileupread.alignment.query_sequence[pileupread.query_position]
                      refBase_TN=DNA_complement_dict[refBase]
                      readBaseRNA=RNA_complement_dict[DNAtoRNA[readBase]]
                      init=pileupread.alignment.reference_start
                      if pileupread.query_position >=chop+1 and pileupread.alignment.qlen - pileupread.query_position >chop : #plus one is to remove the biotin position
                           count+=1
                           count_pos[readBase]+=1
                           BaseCount_dict[DNAtoRNA[refBase_TN]]+=1
                           if refBase != readBase:
                               mutaCount+=1
                               pos_from_3pr=pileupread.query_position

                               m.write(i+"\t"+str(pileupcolumn.pos+1)+"\t"+DNAtoRNA[refBase_TN]+"\t"+str(pileupread.alignment.query_name)+"\t"+readBaseRNA+"\t"+str(pos_from_3pr)+"\t"+str(pileupread.alignment.qlen)+"\t"+str(direction)+"\t"+str(pileupread.alignment.query_qualities[pileupread.query_position])+"\t"+str(pileupread.alignment.mapping_quality)+"\t"+str(init)+"\n")

        pout.write(str(i)+"\t"+str(pileupcolumn.pos+1)+"\t"+str(refBase)+"\t"+str(count)+"\t"+str(mutaCount)+"\t"+str(count_pos["A"])+"\t"+str(count_pos["T"])+"\t"+str(count_pos["G"])+"\t"+str(count_pos["C"])+"\t"+str(count_pos["N"])+"\n")

out=open(workdir+"BaseTypeCount_MQ"+str(min_mapping_qual)+"_BQ"+str(min_base_qual)+".txt","wt")
for key in BaseCount_dict.keys():
     out.write(str(key)+"\t"+str(BaseCount_dict[key])+"\n")
out.close()
pout.close()
m.close()
samfile.close()
