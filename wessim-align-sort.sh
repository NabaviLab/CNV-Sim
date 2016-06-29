#!/usr/bin/env bash

ORIGINAL_GENOME=$1      # used for alignment
GENOME_FILE=$2          # used as a reference for reads by wessim
TARGET_FILE=$3          # used as a target for reads by wessim
EXPERIMENT_NAME=$4      # used to name output files

CNV_HOME=$(pwd)
MODEL_FILE=$CNV_HOME'/input/models/ill100v4_p.gzip'

# Step 2 (WESSIM)
READ_LENGTH=100
NUMBER_OF_THREADS=8
OUTPUT_FILE=$CNV_HOME'/output/wessim_output/'$EXPERIMENT_NAME
NUMBER_OF_TARGETS=$(wc -l $TARGET_FILE | awk '{print($1)}')
NUMBER_OF_READS=$(($NUMBER_OF_TARGETS*10))

cd $CNV_HOME
cd 2_wessim/Wessim/
python Wessim1.py -R $GENOME_FILE -B $TARGET_FILE -n $NUMBER_OF_READS -l $READ_LENGTH -M $MODEL_FILE -z -o $OUTPUT_FILE -t $NUMBER_OF_THREADS -p

# Step 3 (BWA)
NUMBER_OF_THREADS=4
INPUT_FILE_1=$CNV_HOME'/output/wessim_output/'$EXPERIMENT_NAME'_1.fastq.gz'
INPUT_FILE_2=$CNV_HOME'/output/wessim_output/'$EXPERIMENT_NAME'_2.fastq.gz'
OUTPUT_FILE=$CNV_HOME'/output/bwa_output/'$EXPERIMENT_NAME'.sam'

cd $CNV_HOME
cd 3_bwa_aligner/bwa/
# ./bwa index $ORIGINAL_GENOME
./bwa mem -t $NUMBER_OF_THREADS $ORIGINAL_GENOME $INPUT_FILE_1 $INPUT_FILE_2 > $OUTPUT_FILE

# Step 4 (SAMTOOLS)
INPUT_FILE=$CNV_HOME'/output/bwa_output/'$EXPERIMENT_NAME'.sam'
OUTPUT_FILE=$CNV_HOME'/output/samtools_output/'$EXPERIMENT_NAME'.sorted.sam'
TEMP_FILE=$CNV_HOME'/output/samtools_output/'$EXPERIMENT_NAME'.sort'

cd $CNV_HOME
samtools sort -O sam -T $TEMP_FILE -o $OUTPUT_FILE $INPUT_FILE
