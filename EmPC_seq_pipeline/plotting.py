import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.switch_backend('Agg')
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

#Plot the distribution of number of ways to map a transcript in the reference genome
a=np.loadtxt(workdir+"Distribution_of_Ambiguity.txt")
f,ax=plt.subplots(figsize=(12,8))
width=0.6
ax.bar(a[:,0],a[:,1]/sum(a[:,1])*100,width,color="gray")
fontsize=20
ax.set_xticks(np.arange(len(a)))
plt.tick_params("y",labelsize=fontsize)
plt.tick_params("x",labelsize=fontsize)
ax.set_xlim(0.5,max(a[:,0])+0.5)
ax.set_xticks(np.arange(len(a))+1)
ax.set_xticklabels(a[:,0].astype(int))
ax.set_xlabel("Number of ways to map to reference genome",fontsize=fontsize+5)
ax.set_ylabel("Percentage of transcripts",fontsize=fontsize+5)
f.savefig(workdir+"Distribution_NumberOfWaysToMap.png")

#Plot the error rate for each mutation type
MutaType=np.loadtxt(workdir+"MutationTypeSpectrum.txt",dtype=str)
f,ax=plt.subplots(figsize=(12,8))
label=MutaType[:,0]
muta=MutaType[:,1].astype(float)
width=0.6
ax.bar(np.arange(len(MutaType)),muta,width,color="gray")
fontsize=20
ax.set_xticks(np.arange(len(MutaType)))
ax.set_xticklabels(label,fontsize=fontsize)
ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
ax.yaxis.offsetText.set_fontsize(fontsize)
plt.tick_params("y",labelsize=fontsize)
ax.set_ylabel("Mutational Frequency",fontsize=fontsize+5)
f.savefig(workdir+"MutationTypeSpectrum.png")

#Plot the number of errors and coverage along the rDNA
Exp_Genome_Muta={}
Exp_Genome_Cov={}
for keys in chroms:
    Exp_Genome_Muta[keys]=[0]*len(ref_dict[keys])
    Exp_Genome_Cov[keys]=[0]*len(ref_dict[keys])
    Sim_GenomePos=[0]*len(ref_dict[keys])

'''
for lines in FF:
    temp=lines.strip().split()
    line=[]
    if float(temp[2])>0.0:
        line=[float(temp[1]),float(temp[2]),float(temp[3])]
        Sim_Genome_Muta[temp[0]].append(line)
'''
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
    f,ax=plt.subplots(figsize=(20,8))
    fontsize=15
    width=1
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
    ax.plot(Sim_toplot[:,0],Sim_toplot[:,1],color="k",zorder=0,label="Simulation",lw=1)
    ax.fill_between(Sim_toplot[:,0],[0]*len(Sim_toplot[:,0]),Sim_toplot[:,1],color='gray',alpha=0.5)
    for i in range(len(label)):
          ax.plot([label[i],label[i]],[0,muta2[i]],color="red",lw=2,alpha=0.8,zorder=1)
    ax.plot([label[i],label[i]],[0,muta2[i]],color="red",lw=1,zorder=1,label="Experiment")
    exp_plot=np.stack((label,muta2)).T
    np.savetxt("MutationalFrequency_Exp_chrom_"+str(keys)+".txt",exp_plot)
    np.savetxt("MutationalFrequency_Sim_chrom_"+str(keys)+".txt",Sim_toplot)
    hs_rate=[]
    HS=Hotspot[Hotspot[:,0]==keys]
    HS_label=HS[:,1].astype(float)

    for i in range(len(HS)):
         hs_rate.append(muta2[label==int(HS[i][1])])
    ax.scatter(HS[:,1].astype(float),np.array(hs_rate),s=25,marker="o",color="red",edgecolor="k",zorder=2,label="TEEL")
    hs_rate=np.array(hs_rate)[:,0]
    HS2=np.stack((HS_label,hs_rate)).T
    np.savetxt("MutationalFrequency_TEEL_chrom_"+str(keys)+".txt",HS2)

    MinX=[]
    ax.legend(loc=0,fontsize=fontsize)
    ax.set_ylim(0,max(muta_frequency)*1.1)
    ax.set_xlim(min(label)*0.9,max(label)*1.1)
    ax.set_ylabel("Mutational Frequency",fontsize=fontsize)
    ax.set_xlabel(keys,fontsize=fontsize)
    f.savefig(workdir+"Muta_Frequency_inChrom_"+str(keys)+".png")

#Plot the error frequency at each position in the transcripts (from the 3prime end of the transcript).
f,ax=plt.subplots(figsize=(15,8))
ReadPos=np.loadtxt(workdir+"ReadPos_pvalues.txt",dtype=str)[1:].astype(float)

cutoff=30
SimY_axis=ReadPos[:,3][:cutoff]
SimY_axis_error=ReadPos[:,4][:cutoff]
ExpY_axis=ReadPos[:,5][:cutoff]
ExpY_axis_error=ReadPos[:,6][:cutoff]
X_axis=ReadPos[:,0][:cutoff]
ax.errorbar(X_axis, SimY_axis, yerr=SimY_axis_error,color="blue",elinewidth=2,markersize=10,fmt='o',capthick=2,capsize=5,label="Simulation")
ax.errorbar(X_axis, ExpY_axis, yerr=ExpY_axis_error,color="red",elinewidth=2,markersize=10,fmt='o',capthick=2,capsize=5,label="Experiment")
ax.set_xlim(-1.5,30)
MAXs=[max(SimY_axis+SimY_axis_error),max(ExpY_axis+ExpY_axis_error)]
mins1=np.unique(np.sort(SimY_axis-SimY_axis_error))[0]
mins2=np.unique(np.sort(ExpY_axis+ExpY_axis_error))[0]
mins=min(mins1,mins2)
ax.set_ylim(mins-2*mins,max(MAXs)*1.1)
fontsize=25
ax.legend(loc=0,fontsize=fontsize)
plt.tick_params("x",labelsize=fontsize)
plt.tick_params("y",labelsize=fontsize)
ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
ax.yaxis.offsetText.set_fontsize(fontsize)
ax.set_ylabel("Error Rate",fontsize=fontsize+5)
ax.set_xlabel("Distance to 3-Prime end",fontsize=fontsize+5)
ALL=np.stack((X_axis,SimY_axis,SimY_axis_error,ExpY_axis,ExpY_axis_error)).T
np.savetxt("MutationalFrequency_PerPositionInTranscript.txt",ALL)

f.savefig(workdir+"ErrorRate_per_PositionInTranscripts.png")
