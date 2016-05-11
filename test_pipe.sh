#!/bin/sh

#BSUB -q priority
#BSUB -W 12:00  # hours:minutes runlimit after which job wille be killed.
#BSUB -e errors_cnv_2.5M_100x.%J  # send errors (stderr) to file errfile.  If file exists, add to it.
#BSUB -n 4
#BSUB -N
#BSUB -o output_cnv_2.5M_100x.%J

export PATH=$PATH:/groups/nabavi/David/bin
export PATH=$PATH:/opt/python-2.7.6/bin
export PATH=$PATH:/opt/lsf/7.0/linux2.6-glibc2.3-x86_64/bin

module load seq/samtools/1.1 seq/bowtie/2.2.4 dev/python/2.7.6 seq/bwa/0.7.8

while read line
do
	samtools faidx $line.fa
	/groups/nabavi/David/bin/faToTwoBit $line.fa $line.2bit
done <list_A1.txt

while read line
do	
	samtools faidx $line.fa 
	/groups/nabavi/David/bin/faToTwoBit $line.fa $line.2bit
done <list_A2.txt


TARGET=$1
MODEL=/groups/nabavi/David/tools/GemSIM_v1.6/models/ill100v5_p.gzip

### WES simulation ### 

while read line
do
	bsub -q priority -W 12:00 -e errors_cnv_2.5Mx_100.%J -n 4 -N -o output_cnv_2.5M_100x.%J python2.7 Wessim1.py -R $line.fa -B $TARGET -n 2500000 -l 100 -M $MODEL -p -v -t 4 -o WESim_$line
done <list_A1.txt

while read line
do	
	bsub -q priority -W 12:00 -e errors_cnv_2.5Mx_100.%J -n 4 -N -o output_cnv_2.5M_100x.%J python2.7 Wessim1.py -R $line.fa -B $TARGET -n 2500000 -l 100 -M $MODEL -p -v -t 4 -o WESim_$line
done <list_A2.txt

#### Alignment ###########
#######To index your reference genome use:
#bsub -q priority -W 12:00 bwa index -a bwtsw -p hg19_chr1 chr1.fa

#####(Make sure to bwa index your reference genome sequence prior to this step and specify the prefix) ###########

while read line
do
       bsub -n 4 -q priority -W 12:00 -e errors_BWA.%J -o output_BWA.%J bwa aln -t 4 -f $line.sai /groups/nabavi/David/data/chr1_hg19/Hg19_chr1 $line.fastq
done <fastq_1.txt

while read line
do
       bsub -n 4 -q priority -W 12:00 -e errors_BWA.%J -o output_BWA.%J bwa aln -t 4 -f $line.sai /groups/nabavi/David/data/chr1_hg19/Hg19_chr1 $line.fastq
done <fastq_2.txt
