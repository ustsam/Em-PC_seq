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

The data file are suggested to process with triming software for adaptor and low qulaity base before running this script.
## Data analysis
1. The run_noQsfilter.sh script will output data.sam.gz
2. To analyze the data, type the following command:
`bash data_analysis.sh {PATH to the working directory} {PATH to the reference file} {PATH to the script directory} {ambiguity threshold} {maximum depth per site} {minimum base quality} {minimum mapping quality} {number of simulation data generated}`

## Figure plotting

To obtain figures in the manuscript, please run the following three main functions:
`bash plot.sh {PATH to the working directory} {PATH to the reference file} {PATH to the script directory} {ambiguity threshold} {maximum depth per site} {minimum base quality} {minimum mapping quality} 1 `

matplotlib should be installed to plot the figures. if it is not installed, please run the following script:
`bash plot.sh {PATH to the working directory} {PATH to the reference file} {PATH to the script directory} {ambiguity threshold} {maximum depth per site} {minimum base quality} {minimum mapping quality} 0 `

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


