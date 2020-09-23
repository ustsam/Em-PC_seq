#Author: Ilona Christy Unarta (e-mail: icunarta@connect.ust.hk)
#Institution: HKUST (September 2020)
#Objective: Create a set of data based on the binomial distribution fitted to experimental data by maximum likelihood. 
#Then calculate the p-value of this set of data compared with simulation data

import os
import numpy as np
import optparse 

p = optparse.OptionParser()
p.add_option('--reffasta', '-f', default="rDNA1.fasta",  help="reference fasta file")
p.add_option('--workdir','-w',default="./", help="working directory")
p.add_option('--pileup','-p',default="data_sorted_filtered_MQ30_BQ0.pileup_pysam_count")
p.add_option('--mutafile','-m',default="muta_reads_MQ30_BQ30_PosReads.txt")
options, arguments = p.parse_args()
workdir=str(options.workdir)
pileup=str(options.pileup)
ref_filename=str(options.reffasta)
muta_file=str(options.mutafile)

def diffvaf(sim,binom):
	if len(sim)!=len(binom):
		print("length of arrays are not the same.",len(sim),len(binom))
		return 
	s=0
	w=0
	pv=1
	for i in range(len(sim)):
		for j in range(len(sim)):
			w+=1
			if sim[i]>=binom[j]: #H0 is true
				s+=1
	if w>0:
		pv=float(s)/float(w)
	
	return pv


os.system("find "+workdir+"simulation_mutaGenomePos_*.txt >& "+workdir+"list")
file_lists=open(workdir+"list","r").read().split("\n")[:-1]
numSim=len(file_lists)

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
	
for j in range(numSim):
	FF=np.loadtxt(file_lists[j],dtype=str)
	for i in range(len(FF)):
		muta_sim_dict[FF[i][0]][int(FF[i][1])-1][j]=float(FF[i][2])


pvalue=[]
a=np.loadtxt(workdir+pileup,dtype=str)
out=open(workdir+"GenomePos_pvalues.txt","wt")
out.write("Chrom Position Coverage NoOfErrors AverageErrorRate_Simulation StdErrorRate_Simulation AverageErrorRate_Binom StdErrorRate_Binom P_values\n")
for keys in ref_dict:
	for i in range(len(a)):
		line=[]
		cov=float(a[i][3])
		muta=int(a[i][4])
		pos=int(a[i][1])
		if cov>0:
			random_array=np.random.binomial(int(a[i][3]),float(a[i][4])/cov,size=numSim)
			simulation=np.array(muta_sim_dict[keys][pos-1])/cov
			binom=random_array/cov
			average_simulation=np.mean(simulation)
			std_simulation=np.std(simulation)
			average_binom=np.mean(binom)
			std_binom=np.std(binom)
			pv=diffvaf(simulation,binom)
			out.write(str(keys)+" "+str(pos)+" "+str(cov)+" "+str(muta)+" "+str(average_simulation)+" "+str(std_simulation)+" "+str(average_binom)+" "+str(std_binom)+" "+str(pv)+"\n")


os.system("find "+workdir+"simulation_mutaReadPos_*.txt >& "+workdir+"list")
file_lists=open(workdir+"list","r").read().split("\n")[:-1]
numSim=len(file_lists)

muta_pos_dict={}
for i in range(1000):
	muta_pos_dict[i]=[]

for j in range(numSim):
	FF=np.loadtxt(file_lists[j],dtype=str)
	for i in range(len(FF)):
		muta_pos_dict[int(FF[i][0])].append(float(FF[i][1]))

SL=np.loadtxt(workdir+"startlengths.txt",dtype=str)[:,2].astype(int)
histo,bins=np.histogram(SL,bins=np.arange(min(SL),max(SL)+1))
depth_per_readpos=[0]*(max(bins)+1)
bins=bins[:-1]-1
depth_per_readpos[bins[len(bins)-1]]=histo[len(histo)-1]
for i in range(len(histo)-2,-1,-1):
	depth_per_readpos[bins[i]]=histo[i]+depth_per_readpos[bins[i+1]]
for i in range(bins[0]+1):
	depth_per_readpos[i]=len(SL)

FF=np.loadtxt(muta_file)#"muta_reads_MQ30_BQ30_PosReads.txt")
Exp_MutaPos=[0]*len(depth_per_readpos)
for i in range(len(FF)):
	Exp_MutaPos[int(FF[i][0])]=float(FF[i][1])

out=open(workdir+"ReadPos_pvalues.txt","wt")
out.write("PositionInTranscript CoverageInPosition NoOfErrors AverageErrorRate_Simulation StdErrorRate_Simulation AverageErrorRate_Binom StdErrorRate_Binom P_values\n")
for i in range(len(depth_per_readpos)):
	if depth_per_readpos[i]>0 :
		random_array=np.random.binomial(depth_per_readpos[i],float(Exp_MutaPos[i])/depth_per_readpos[i],size=numSim)
		simulation=np.array(muta_pos_dict[i])/float(depth_per_readpos[i])
		binom=random_array/float(depth_per_readpos[i])
		pv=diffvaf(simulation,binom)
		average_simulation=np.mean(simulation)
		std_simulation=np.std(simulation)
		average_binom=np.mean(binom)
		std_binom=np.std(binom)
		out.write(str(i)+" "+str(depth_per_readpos[i])+" "+str(Exp_MutaPos[i])+" "+str(average_simulation)+" "+str(std_simulation)+" "+str(average_binom)+" "+str(std_binom)+" "+str(pv)+"\n")
