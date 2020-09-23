#Author: Ilona Christy Unarta (e-mail: icunarta@connect.ust.hk)
#Institution: HKUST (September 2020)
#Objective: Generate simulation file with equally distribution mutation along a transcript


import random as rd
import optparse
import numpy as np
#import statistics as stat
import math 
import gzip

p = optparse.OptionParser()
p.add_option('--reffasta', '-f', default="rDNA1.fasta",  help="reference fasta file")
p.add_option('--nummuta', '-m', default=1, help="number of mutation in experiment")
p.add_option('--numcov', '-c', default=1, help="number of total coverage in experiment")
p.add_option('--numSim', '-s',default=100, help="the number of simulations performed")
p.add_option('--rotation','-r', action="store_true",dest="rotation")
p.add_option('--workdir','-w',default="./", help="working directory")
p.add_option('--pileup','-p',default="data_sorted_filtered_MQ30_BQ0.pileup_pysam_count")

options, arguments = p.parse_args()
ref_filename=str(options.reffasta)
MUTA=float(options.nummuta)
COV=float(options.numcov)
rotation=bool(options.rotation)
numSim=int(options.numSim)
workdir=str(options.workdir)
if rotation:
    print("Rotation is turned on. it will do circularization and two copies of the read are printed in the new read.")
else :
    print("Rotation is turned off. it will not do circularization and only one copy of the read is printed.")

def mut1base(bw):
    N = ['A','T','G','C']
    if bw in N :
         N.remove(bw)#[N.index(bw)]=[]
    else :
         error('unknown input base!')
    bm = N[rd.randint(0,len(N)-1)]
    return bm


mutrate=float(MUTA)/float(COV)
print("mutation rate:",mutrate)
#Read reference fasta file 
m=open(ref_filename,"r")
ref_dict={}
for line in m :
    temp=line.strip().split()
    if temp[0][0]==">":
        chrs=temp[0][1:]
        ref_dict[temp[0][1:]]=""
    else:
        ref_dict[chrs]+=temp[0]
m.close()

muta_sim_dict={}
for keys in ref_dict:
    muta_sim_dict[keys]=[[]]*len(ref_dict[keys])
    for i in range(len(muta_sim_dict[keys])):
         muta_sim_dict[keys][i]=[0]*numSim

muta_pos_dict={}
for i in range(1000):
    muta_pos_dict[i]=[0]*numSim

#Read start and length table file
f=open(workdir+"startlengths.txt","r")
for line in f:
    temp=line.strip().split()
    chrom=str(temp[0])
    start=int(temp[1])
    length=int(temp[2])
    read=ref_dict[chrom][start-1:start+length-1]

    for sim in range(numSim):
            strsim=str(sim)
            outname=strsim.zfill(5)
            fout=gzip.open(workdir+"sim_"+outname+".fastq.gz","a+")
            #fout2=open(workdir+"mutaPos_"+outname+".txt","a+")

            #mutating
            random_array=np.random.rand(length)#[rd.random() for _ in range(length)]
            temp2=np.where(random_array<mutrate)
            nmut=len(temp2[0])
            mut_pos=temp2[0]
            readname="@"+str(chrom)+":"+str(start)+"-"+str(start+length-1)
            mutname=""
            if nmut>0:
                 for j in range(nmut):
                     refbase=ref_dict[chrom][start+mut_pos[j]]
                     mutbase=mut1base(refbase)
                     newread=read[:mut_pos[j]+1]+mutbase+read[mut_pos[j]+2:]
                     read=newread
                     mutname+="$mut:"+str(start+mut_pos[j])+":"+str(refbase)+">"+str(mutbase) 
                     muta_sim_dict[chrom][start+mut_pos[j]-1][sim]+=1
                     muta_pos_dict[mut_pos[j]][sim]+=1
             #        fout2.write(str(chrom)+" "+str(start+mut_pos[j])+" "+str(mut_pos[j])+" "+str(refbase)+">"+str(mutbase)+"\n")

            #rotating 
            if rotation:
                cirpos=rd.randint(0,length-1)
                newread=read[cirpos:]+read[:cirpos] 
                read=newread
                triplereads=read*2
                readname+="#cir:"+str(start+cirpos)

            readname+=mutname
            fout.write(readname+"\n")
            if rotation:
                fout.write(triplereads+"\n")
                fout.write("+\n")
                qual="H"*len(triplereads)
                fout.write(qual+"\n")
            else:
                fout.write(read+"\n")
                fout.write("+\n")
                qual="H"*len(read)
                fout.write(qual+"\n")

            fout.close()
#            fout2.close()

#fout3=open(workdir+"simulation_mutaGenomePos_average.txt","wt")
#fout4=open(workdir+"simulation_mutaReadPos_average.txt","wt")
for sim in range(numSim):
    strsim=str(sim)
    outname=strsim.zfill(5)
    fout3=open(workdir+"simulation_mutaGenomePos_"+outname+".txt","a+")
    fout4=open(workdir+"simulation_mutaReadPos_"+outname+".txt","a+")
    for keys in muta_sim_dict:
         for j in range(len(muta_sim_dict[keys])):
              fout3.write(str(keys)+" "+str(j+1)+" "+str(muta_sim_dict[keys][j][sim])+"\n")
    for keys in muta_pos_dict:
        fout4.write(str(keys)+" "+str(muta_pos_dict[keys][sim])+"\n")
    fout3.close()
    fout4.close()

