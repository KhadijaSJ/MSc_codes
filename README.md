## Table of contents
* [Introduction](#Introduction)
* [Technologies](#technologies)
* [Setup](#setup)

## Introduction
This project is based on identifying potential biomarkers for a disease, in this case it is pancreatic cancer. The required inputs are raw fastq files that were pre-processed, processed and analysed accordingly. The output from DESeq2 and/or ClusterProfiler consists of potential genes that needs further in-silico and wet laboratory validation.  

### Pre-processing:
A python pipeline consisting of different packages, see preprocessing.py
	
## Technologies
Project used the following languages:
* Python
* R studio
	
## Setup
To run this project, the following packages will be used:
### Python
* FastQC
* Trimmomatic
* HISAT2
* HTSeq
### R
* DESeq2
* ClusterProfiler
* Random forest classifier


