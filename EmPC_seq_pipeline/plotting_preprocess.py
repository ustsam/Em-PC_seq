import numpy as np
import sys
import os 


file1=str(sys.argv[1])
file2=str(sys.argv[2])
file3=str(sys.argv[3])
workdir=str(sys.argv[4])
ref_filename=str(sys.argv[5])
GG=open(workdir+file1,"r")
HH=open(workdir+file2,"r")

os.system("awk '{print $1}' "+workdir+file2+" | sort | uniq > chrom.list")
chroms=open(workdir+"chrom.list","r").read().split("\n")[:-1]

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

#Plot the number of errors and coverage along the rDNA
Exp_Genome_Muta={}
Exp_Genome_Cov={}
for keys in chroms:
    Exp_Genome_Muta[keys]=[0]*len(ref_dict[keys])
    Exp_Genome_Cov[keys]=[0]*len(ref_dict[keys])

for lines in GG:
    temp=lines.strip().split()
    Exp_Genome_Cov[temp[0]][int(temp[1])]=float(temp[3])
    
for lines in HH:
    temp=lines.strip().split()
    Exp_Genome_Muta[temp[0]][int(temp[1])]=float(temp[4])
    
for keys in Exp_Genome_Muta:
    Exp_Genome_Muta[keys]=np.array(Exp_Genome_Muta[keys])
    Exp_Genome_Cov[keys]=np.array(Exp_Genome_Cov[keys])

GenomePos=np.loadtxt(workdir+"GenomePos_pvalues.txt",dtype=str)[1:]
Hotspot=GenomePos[GenomePos[:,-1].astype(float)<=0.05]
Sim_GenomePos=GenomePos[:,[0,1,4,5]]

for keys in Exp_Genome_Muta:
    Exp_Genome_Cov[keys][Exp_Genome_Cov[keys]==0]=1
    Sim_GenomePos_key=Sim_GenomePos[Sim_GenomePos[:,0]==keys][:,[1,2,3]].astype(float)
    Sim_GenomePos_key=Sim_GenomePos_key[Sim_GenomePos_key[:,1]>0.0]
    muta_frequency=Exp_Genome_Muta[keys]/Exp_Genome_Cov[keys]
    muta2=muta_frequency[muta_frequency>0]    
    Sim_toplot=[]
    for i in range(len(Sim_GenomePos_key)):
         if Exp_Genome_Muta[keys][int(Sim_GenomePos_key[i][0])]>0:
             Sim_toplot.append(Sim_GenomePos_key[i])
    Sim_toplot=np.array(Sim_toplot)

    label=np.arange(len(muta_frequency))
    label=label[muta_frequency>0]
    exp_plot=np.stack((label,muta2)).T
    np.savetxt("MutationalFrequency_Exp_chrom_"+str(keys)+".txt",exp_plot)
    np.savetxt("MutationalFrequency_Sim_chrom_"+str(keys)+".txt",Sim_toplot)
    HS=Hotspot[Hotspot[:,0]==keys]
    HS_label=HS[:,1].astype(float)
    hs_rate=[]
    for i in range(len(HS)):
         hs_rate.append(muta2[label==int(HS[i][1])])
    hs_rate=np.array(hs_rate)[:,0]
    HS2=np.stack((HS_label,hs_rate)).T
    np.savetxt("MutationalFrequency_TEEL_chrom_"+str(keys)+".txt",HS2)

#Plot the error frequency at each position in the transcripts (from the 3prime end of the transcript).
ReadPos=np.loadtxt(workdir+"ReadPos_pvalues.txt",dtype=str)[1:].astype(float)

cutoff=30
SimY_axis=ReadPos[:,3][:cutoff]
SimY_axis_error=ReadPos[:,4][:cutoff]
ExpY_axis=ReadPos[:,5][:cutoff]
ExpY_axis_error=ReadPos[:,6][:cutoff]
X_axis=ReadPos[:,0][:cutoff]
ALL=np.stack((X_axis,SimY_axis,SimY_axis_error,ExpY_axis,ExpY_axis_error)).T
np.savetxt("MutationalFrequency_PerPositionInTranscript.txt",ALL)

