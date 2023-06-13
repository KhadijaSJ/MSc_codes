import os
import glob
dir = '/media/gh527/Seagate/Khadija/masters/dataset3/'
arr = os.listdir(dir)

#fastqc: analysing raw fastq files
def qc():
    for i in arr:
        ar = os.listdir(dir + i)
        for a in ar:
            print(a)
            os.system('fastqc /media/gh527/Seagate/Khadija/masters/dataset3/'+ i +'/'+ a)

#trimmomatic: trimming and removing low quality reads
def trim():
    for i in arr:
        os.chdir(dir + i)
        a = glob.glob('./*.fastq.gz')
        os.chdir('/home/gh527/Downloads/Trimmomatic-0.39')
        os.system('java -jar trimmomatic-0.39.jar PE ' + dir + i +'/'+ a[0] + ' ' + dir + i + '/' +a[1] + ' /media/gh527/Seagate/Khadija/masters/dataset3/trim/' + i +'/forward_'+ i + '.fastq.gz'+ ' /media/gh527/Seagate/Khadija/masters/dataset3/trim/' + i +'/reverse_' + i +'.fastq.gz ILLUMINACLIP:/home/gh527/Downloads/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36')

#hisat2 alignment: alignment to the human reference genome and save alignment rate
def hisat1():
    for i in arr:
        os.chdir(dir + i)
        a = glob.glob('*.fastq.gz')
        os.chdir('/home/gh527/Downloads/hisat2-2.2.1')
        os.system( 'hisat2 -p 12 -x /media/gh527/Seagate/Khadija/masters/genome/genome -1 /media/gh527/Seagate/Khadija/masters/dataset3/'+ i +'/' + a[0] + ' -2 /media/gh527/Seagate/Khadija/masters/dataset3/'+ i +'/' + a[1] + ' -S /media/gh527/Seagate/Khadija/masters/dataset3/sam/' + i + '.sam --summary-file /media/gh527/Seagate/Khadija/masters/dataset3/sam_alignment/' + i)

#SAM to BAM: conversion of SAM to BAM files
def sam_bam():
    for i in arr:
        os.mkdir('/media/gh527/Seagate/Khadija/masters/dataset3/sort_bam/' + i)
        os.system('samtools view -@ 12 -S -b /media/gh527/Seagate/Khadija/masters/dataset3/sam/' + i + '.sam > /media/gh527/Seagate/Khadija/masters/dataset3/bam/' + i + '.bam')
        os.system('samtools sort  -@ 12 /media/gh527/Seagate/Khadija/masters/dataset3/bam/' + i + '.bam -o /media/gh527/Seagate/Khadija/masters/dataset3/sort_bam/' + i + '/' + i + '.bam')
        os.system('samtools index -@ 12 /media/gh527/Seagate/Khadija/masters/dataset3/sort_bam/' + i + '/' + i + '.bam /media/gh527/Seagate/Khadija/masters/dataset3/sort_bam/' + i + '/' + i + '.bai')

#bam quality control
def qbam():
    for i in arr:
        os.system('samtools stats -@ 12 /media/gh527/Seagate/Khadija/masters/dataset3/sort_bam/' + i + '/' + i + '.bam > /media/gh527/Seagate/Khadija/masters/dataset3/bam_stats/' + i + '.stats')
        os.system('plot-bamstats -p /media/gh527/Seagate/Khadija/masters/dataset3/bam_stats/' + i + '.stats > /media/gh527/Seagate/Khadija/masters/dataset3/bam_quality/' + i)

#HTSeq: producing a count matrix from the fastq files
def htseq():
    for i in arr:
        os.system('htseq-count -f bam -s yes /media/gh527/Seagate/Khadija/masters/dataset3/sort_bam/' + i + '/' + i + '.bam /media/gh527/Seagate/Khadija/masters/gencode.v40.chr_patch_hapl_scaff.annotation.gff3 > /media/gh527/Seagate/Khadija/masters/dataset3/count/' + i + '.txt')

if __name__ == '__main__':
    qc()
    trim()
    hisat()
    sam_bam()
    #qbam()
    bam()
    htseq()
