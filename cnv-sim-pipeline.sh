#!/bin/bash

CNV_HOME=$(pwd)
EXPERIMENT_NAME='my_test'

GENOME_FILE=$CNV_HOME'/input/test/chr1.fa'
TARGET_FILE=$CNV_HOME'/input/test/test.bed'
MODEL_FILE=$CNV_HOME'/input/models/ill100v4_p.gzip'

# Step 1 (CNV SIM)
sort -k1,1 -k2,2n $TARGET_FILE > $TARGET_FILE'.sorted'
bedtools merge -i $TARGET_FILE'.sorted' > $TARGET_FILE'.sorted.merged'
echo 'sorted and merged exons target'

cd $CNV_HOME
cd 1.\ cnv_simulator
python cnv-sim.py $GENOME_FILE $TARGET_FILE'.sorted.merged' $EXPERIMENT_NAME

<<COMMENT
# Step 2 (WESSIM)
NUMBER_OF_READS=2000
READ_LENGTH=100
NUMBER_OF_THREADS=8
OUTPUT_FILE=$CNV_HOME'/output/wessim_output/'$EXPERIMENT_NAME

cd $CNV_HOME
cd 2.\ wessim/Wessim/
python Wessim1.py -R $GENOME_FILE -B $TARGET_FILE -n $NUMBER_OF_READS -l $READ_LENGTH -M $MODEL_FILE -z -o $OUTPUT_FILE -t $NUMBER_OF_THREADS -p

# Step 3 (BWA)
NUMBER_OF_THREADS=8
INPUT_FILE_1=$CNV_HOME'/output/wessim_output/'$EXPERIMENT_NAME'_1.fastq.gz'
INPUT_FILE_2=$CNV_HOME'/output/wessim_output/'$EXPERIMENT_NAME'_2.fastq.gz'
OUTPUT_FILE=$CNV_HOME'/output/bwa_output/'$EXPERIMENT_NAME'.sam'

cd $CNV_HOME
cd 3.\ bwa_aligner/bwa/
./bwa index $GENOME_FILE
./bwa mem -t $NUMBER_OF_THREADS $GENOME_FILE $INPUT_FILE_1 $INPUT_FILE_2 > $OUTPUT_FILE

# Step 4 (SAMTOOLS)
INPUT_FILE=$CNV_HOME'/output/bwa_output/'$EXPERIMENT_NAME'.sam'
OUTPUT_FILE=$CNV_HOME'/output/samtools_output/'$EXPERIMENT_NAME'.sorted.sam'
TEMP_FILE=$CNV_HOME'/output/samtools_output/'$EXPERIMENT_NAME'.sort'

cd $CNV_HOME
samtools sort -O sam -T $TEMP_FILE -o $OUTPUT_FILE $INPUT_FILE

# Step 5 (IGV)
COMMENT