from numpy import * 
import os 
import sys
filename=sys.argv[1]
a=open(filename,"r")
out=open("relocalized.fastq","wt")
for line in a:
	temp=line.strip().split()
	seq=temp[9]
	seq_array=[]
	qual_array=[]
	seq_array.append(seq)
	qual_array.append(temp[10])
	qual=temp[10]
	for i in range(10):
		seq2=seq[i+1:len(seq)]+seq[:i+1]
		seq_array.append(seq2)
		qual_array.append(qual[i+1:len(seq)]+qual[:i+1])
	for i in range(10):
		seq2=seq[len(seq)-i-1:]+seq[:len(seq)-i-1]
		seq_array.append(seq2)
		qual_array.append(qual[len(qual)-i-1:]+qual[:len(qual)-i-1])
	for i in range(len(seq_array)):
		out.write("@"+str(temp[0])+"\n")#+"_"+str(i)+"\n")	
		out.write(str(seq_array[i])+"\n")
		out.write("+"+str(temp[0])+"\n")#"_"+str(i)+"\n")
		out.write(str(qual_array[i])+"\n")
			
