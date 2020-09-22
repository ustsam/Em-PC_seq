import sys,gzip
##One change, Illon correct that I and D would not be count in XM(bowtie2) or nM(STAR) field, so i should consider it as well.
##Sorting problems for all the same AS score, cause sorting defect,must be fixed, SAM, 20180926
workdir = sys.argv[1]
filename = sys.argv[2]
PRE=sys.argv[3]
infiles = ["/"+filename]
outfile = gzip.open(workdir + "/11_alignment.sam.gz","wb")
outfile1 = gzip.open(workdir + "/13_Enderror.sam.gz","wb")
outfile2 = gzip.open(workdir + "/14_hardclipped.sam.gz","wb")

for file in infiles:
	infile = gzip.open(workdir + file,"rb")
	
##running lists
	Alignment = []
	AlignmentScore = []
	nM = []	#For if condition to screen Spliced read with difference AS
##running lists	
#current read
	SequenceID = ""
	MDtag=[]
	
	while True : #for line in infile:		
		RawLine = infile.readline().strip()
		line = RawLine.split()
		if RawLine=="":
                        if len(Alignment) > 0:
                                BestIndex = []
                                g = 0
                                subsetNM=[]
                                BestIndex_temp=[]
                                while g < len(AlignmentScore):
                                        if AlignmentScore[g]==max(AlignmentScore):
                                                subsetNM.append(nM[g])
                                                BestIndex_temp.append(g)
                                        g+=1
                                for i in range(len(subsetNM)):
                                        if subsetNM[i]==min(subsetNM):
                                                BestIndex.append(BestIndex_temp[i])
                                uniqueness = 1
                                if len(BestIndex)==2:
                                        for i in range(len(BestIndex)):
                                                if MDtag[BestIndex[i]][0]=="0":
                                                        uniqueness=1
                                                        newBest=BestIndex[i]
                                                        BestIndex=[]
                                                        BestIndex.append(newBest)
                                                        break
                                        if len(BestIndex)==2:
                                                uniqueness=len(BestIndex)
                                elif len(BestIndex)>2:
                                        uniqueness=len(BestIndex)
                                if len(BestIndex) > 0 :
                                        Alignment[BestIndex[0]] = Alignment[BestIndex[0]].split("\n")[0]
                                        Alignment[BestIndex[0]] = Alignment[BestIndex[0]]+"\tPR:i:"+str(PRE)+"\tUN:i:"+str(uniqueness)+'\n'
                                        Splitchecker = Alignment[BestIndex[0]]
                                        Splitchecker = Splitchecker.split()
                                        if Splitchecker[5].count("S") == 0 and Splitchecker[5].count("H") == 0:
                                                outfile.write(Alignment[BestIndex[0]])
                                        elif Splitchecker[5].count("S") > 0 and Splitchecker[5].count("H") == 0:
                                                outfile1.write(Alignment[BestIndex[0]])
                                        elif Splitchecker[5].count("H") > 0:
                                                outfile2.write(Alignment[BestIndex[0]])
			break

		if line[0] == SequenceID : 
			for i in range(11,len(line)):
				temp2=line[i].split(":")
				if temp2[0]=="MD":
					MDtag.append(str(temp2[2]))
				elif temp2[0]=="NM":
					MuTag=temp2
				elif temp2[0]=="AS":
					AlignmentScore.append(int(temp2[2]))
			numbers = ""
			letters = []
			for character in line[5]:       #Chop according to the number field 6
				if character == "M" or character == "S" or character == "N" or character == "I" or character == "D" or character == "H" or character == "P":
					letters.append(character)
					numbers += "X"
				else:
					numbers += character
			numbers = numbers.split("X")
			i = 0
			while i < len(letters):
				if letters[i] == "S" or letters[i] == "I" or letters[i] == "D" or letters[i] == "H":
					MuTag[2] = int(MuTag[2]) + int(numbers[i])    #!int! is needed since what append in above is string type.
				i+=1
			nM.append(int(MuTag[2]))
			Alignment.append(RawLine)

		else:
			if len(Alignment) > 0:
				BestIndex = []
				g = 0
				subsetNM=[]
				BestIndex_temp=[]
				while g < len(AlignmentScore):
					if AlignmentScore[g]==max(AlignmentScore):
						subsetNM.append(nM[g])
						BestIndex_temp.append(g)
					g+=1
				for i in range(len(subsetNM)):
					if subsetNM[i]==min(subsetNM):
						BestIndex.append(BestIndex_temp[i])

##Activate unique read only or not, hard-code for activate first
				uniqueness = 1

##Ingore all read that are not unique, less strigent way would be accept the read having error spot unmodified, the following now is the most strigent one.
##Suggestion from illona, print field and can be screen after, no need to relocalize every time.

				if len(BestIndex)==2:
					for i in range(len(BestIndex)):
						if MDtag[BestIndex[i]][0]=="0":
							uniqueness=1
							newBest=BestIndex[i]
							BestIndex=[]
							BestIndex.append(newBest)
							break
					if len(BestIndex)==2:
						uniqueness=len(BestIndex)
				elif len(BestIndex)>2:
					uniqueness=len(BestIndex)
				if len(BestIndex) > 0 :
					Alignment[BestIndex[0]] = Alignment[BestIndex[0]].split("\n")[0]
                                        Alignment[BestIndex[0]] = Alignment[BestIndex[0]]+"\tPR:i:"+str(PRE)+"\tUN:i:"+str(uniqueness)+'\n'
					Splitchecker = Alignment[BestIndex[0]]
					Splitchecker = Splitchecker.split()
					if Splitchecker[5].count("S") == 0 and Splitchecker[5].count("H") == 0:
						outfile.write(Alignment[BestIndex[0]])
					elif Splitchecker[5].count("S") > 0 and Splitchecker[5].count("H") == 0:
						outfile1.write(Alignment[BestIndex[0]])
					elif Splitchecker[5].count("H") > 0:
						outfile2.write(Alignment[BestIndex[0]])

			SequenceID = line[0] 
			Alignment = []
			AlignmentScore = []
			nM = []
			MDtag=[]
			for i in range(11,len(line)):
				temp2=line[i].split(":")
				if temp2[0]=="MD":
					MDtag.append(str(temp2[2]))
				elif temp2[0]=="NM":
					MuTag=temp2
				elif temp2[0]=="AS":
					AlignmentScore.append(int(temp2[2]))

			numbers = ""
			letters = []
			for character in line[5]:       #Chop according to the number field 6
				if character == "M" or character == "S" or character == "N" or character == "I" or character == "D" or character == "H" or character == "P":
					letters.append(character)
					numbers += "X"
				else:
					numbers += character
			numbers = numbers.split("X")

			i = 0
			while i < len(letters):
				if letters[i] == "S" or letters[i] == "I" or letters[i] == "D":
					MuTag[2] = int(MuTag[2]) + int(numbers[i])    #!int! is needed since what append in above is string type.
				i+=1
			nM.append(int(MuTag[2]))
			Alignment.append(RawLine)

