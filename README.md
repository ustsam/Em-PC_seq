# Em-PC_seq
Scripts for treating Em-PC-seq data

## Ownership
[Xuhui HUANG Lab at HKUST]  
[Jihuang WANG Lab at HKUST]

## Status
Active Development

## Introduction
Circular sequencing was initially developed by [Andino Group at UCSF](https://andino.ucsf.edu/) in 2014 for detecting low-frequency variants in the study of RNA virus by reducing sequencing errors [1]. The authors also released their data analysis software, namely CirSeq, which is open-access at https://andino.ucsf.edu/CirSeq.

However, the original was written in hard-code form for repeat and read length. The code is rewritten to perform >=2 repeat identification and variable readlength. The relocalization scripts are rewritten to fit the use of bwa-mem as mapper. Selection of the relocalized read is also changed to be more stringent parameter.

## Demo data

Script for generating simulation data is called simulation.py. It is a python script based on the matlab script written by Dr. Biaobin Jiang from Prof. Jiguang Wang's group at HKUST. original script can be obtained https://github.com/ustsam/emProc-seq in the simulation_script folder. A file called, sim_00000.fastq.gz, is provided as one set of demo data, typical run time would be around 20-30 minutes.

Demo data can be run as following:

1. Download the zipped script file.
2. Unzip the file.
3. Enter script directory
3. Compile the codes by typing `python setup_newreloc.py build_ext --inplace`.
4. `mkdir demo_data_real_time`
5. Call the function using `./run_noQsfilter.sh ./demo_data_real_time ./reference/rDNA1.fasta ./ DUMMY 2 600 ./demo_data/sim_00000.fastq.gz `

## Usage
1. Download the zipped script file.
2. Unzip the file.
3. Enter script directory
3. Compile the codes by typing `python setup_newreloc.py build_ext --inplace`.
4. Call the function using 
`bash run_noQsfilter.sh {PATH of the output directory} {PATH of the reference file} {PATH of the script directory} DUMMY 2 ${twice of the max readlength} ${PATH of the data file in gzipped form}`
notes:
1. The data file are suggested to process with triming software for adaptor and low qulaity base before running this script.

## Data analysis
1. The run_noQsfilter.sh script will output data.sam.gz
2. To analyze the data, type the following command:
`bash data_analysis.sh {PATH to the working directory} {PATH to the reference file} {PATH to the script directory} 1 {maximum depth per site} {minimum base quality} {minimum mapping quality} {number of simulation fastq files generated}`

notes:
1. ambiguity is defined as the number of ways a transcript can be mapped to the reference genome. Ambiguity occurs due to the circularization method, which remove the information about the starting point of the transcript. if ambiguity threshold is set to one, only the transcripts having only one way to be mapped are considered in analysis. In the originial publication (https://doi.org/10.1016/j.jmb.2020.04.011), ambiguity threshold is 1. 

## Figure plotting

To obtain figures in the manuscript, please run the following three main functions:
`bash plot.sh {PATH to the working directory} {PATH to the reference file} {PATH to the script directory} {ambiguity threshold} {maximum depth per site} {minimum base quality} {minimum mapping quality} 1 `

matplotlib should be installed to plot the figures. if it is not installed, please run the following script:
`bash plot.sh {PATH to the working directory} {PATH to the reference file} {PATH to the script directory} {ambiguity threshold} {maximum depth per site} {minimum base quality} {minimum mapping quality} 0 `

## Output files:
The output of plot.sh files figures (.png format) are the following:
a)	“Distribution_NumberOfWaysToMap.png” : the distribution of the number of ways the transcripts can be mapped to the reference genome (ambiguity). The y-axis is the number of transcripts and the x-axis is ambiguity.

b)	“MutationTypeSpectrum.png” : The mutational frequency for each type of mutation  in the RNA transcript. The mutational frequency is the number of errors divided by the coverage of the corresponding reference base. For example: number of A>C errors divided by coverage of base A.

c)	“Muta_Frequency_inChrom_###.png”: The mutational frequency along the sites in the chromosome for the experimental and simulation data. “###" is the name of the chromosome. The transcription error-enriched genomic loci (TEEL) is shown as red dots.

d)	“ErrorRate_per_PositionInTranscripts.png”: The error rate at each position in the transcript. The Position 0 corresponds to the 3’ end of the transcript.

If one would prefer to use othe software tools, the output text files to plot the figures are:
a)	“Distribution_of_Ambiguity.txt” can be used to plot Figure a).

b)	“MutationTypeSpectrum.txt” can be used plot Figure b). 

c)	For Figure c), “MutationalFrequency_Exp_chrom_###.txt” and “MutationalFrequency_Sim_chrom_###.txt” are the files containing the Mutational Frequencies per site in a chromosome for experimental and simulation,respectively. “MutationalFrequency_TEEL_chrom_rDNA1.txt” contains the mutational frequency of the sites considered as TEEL. The first column is the position in the chromosome and the second column is the mutational frequency. 

d)	“MutationalFrequency_PerPositionInTranscript.txt” can be used to plot Figure d). The first column is the position in the transcript. The second and third column is the average and standard deviation of the simulation data. The fourth and fifth column is the average and standard deviation of the experimental data. The average and standard deviation of experimental data is calculated from error rate binomial distribution estimated by maximum likelihood.

## Methodology:

Please check methodology.pdf for detail of the script and the meaning of all generated files. Paper is not yet publish, link will be updated as soon as possible.

## System requirements

The following OS had been used for running the script.

Ubuntu 16.04 LTS (Xenial Xerus)

CentOS release 6.6 (Final)


The following packages are prerequisites for using emProc-seq

1. Python (version 2.7.12)    
2. Cython (version 0.23.4)   
3. NumPy (version 1.11.0)     
4. SciPy (version 0.17.0)    
5. bwa (version 0.7.17-r1188)   
6. samtools (version 1.4.1)
7. pysam (version 0.15.0)
8. matplotlib (version 2.2.2)

NOTE 1: Cython requires a compiler. For OSX this may require installation of Xcode.

NOTE 2: bwa and samtools binaries must be in the PATH.

NOTE 3: Install time are estimated to be no more than 1.5 hours.


## Reference
[1] Acevedo, A., Brodsky, L., & Andino, R. (2014). Mutational and fitness landscapes of an RNA virus revealed through population sequencing. Nature, 505, pp.686-690.

## Contact
For technical questions, please contact TinHang via email: thchong@connect.ust.hk, Ilona via email: icunarta@connect.ust.hk


