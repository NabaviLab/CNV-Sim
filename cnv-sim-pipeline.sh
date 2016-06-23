#!/bin/bash

CNV_HOME=$(pwd)
EXPERIMENT_NAME='my_test_optimized3'

ORIGINAL_GENOME_FILE=$CNV_HOME'/input/test/chr1.fa'
ORIGINAL_TARGET_FILE=$CNV_HOME'/input/test/chr1-target.bed'

# Step 1 (CNV SIM)
sort -k1,1 -k2,2n $ORIGINAL_TARGET_FILE > $ORIGINAL_TARGET_FILE'.sorted'
bedtools merge -i $ORIGINAL_TARGET_FILE'.sorted' > $ORIGINAL_TARGET_FILE'.sorted.merged'
echo 'sorted and merged exons target'

cd $CNV_HOME
cd 1.\ cnv_simulator
python cnv-sim.py $ORIGINAL_GENOME_FILE $ORIGINAL_TARGET_FILE'.sorted.merged' $EXPERIMENT_NAME

cd $CNV_HOME
GENOME=$CNV_HOME'/output/cnvsim_output/'$EXPERIMENT_NAME'-ControlGenome.fa'
TARGET=$CNV_HOME'/output/cnvsim_output/'$EXPERIMENT_NAME'-ControlTarget.bed'
NUMBER_OF_READS=200000
sh wessim-align-sort-vis.sh $ORIGINAL_GENOME_FILE $GENOME $TARGET $EXPERIMENT_NAME'.normal' $NUMBER_OF_READS

GENOME=$CNV_HOME'/output/cnvsim_output/'$EXPERIMENT_NAME'-CNVGenome.fa'
TARGET=$CNV_HOME'/output/cnvsim_output/'$EXPERIMENT_NAME'-CNVTarget.bed'
NUMBER_OF_READS=400000
sh wessim-align-sort-vis.sh $ORIGINAL_GENOME_FILE $GENOME $TARGET $EXPERIMENT_NAME'.cnv' $NUMBER_OF_READS

# Step 5 (IGV)
CNV_FILE=$CNV_HOME'/output/cnvsim_output/'$EXPERIMENT_NAME'-CNVList.bed'
SAMPLE_AMPLIFICATION_FILE=$CNV_HOME'/output/igv_output/'$EXPERIMENT_NAME'-CNVList.amplifications.bed'
SAMPLE_DELETION_FILE=$CNV_HOME'/output/igv_output/'$EXPERIMENT_NAME'-CNVList.deletions.bed'
awk '($4 == 10)' CNV_FILE | head -20 > SAMPLE_AMPLIFICATION_FILE
awk '($4 == -1)' CNV_FILE | tail -20 > SAMPLE_DELETION_FILE

OUTPUT_PATH=$CNV_HOME'/output/igv_output/'$EXPERIMENT_NAME'-amplifications-visualizations/'
bedtools igv -path $OUTPUT_PATH -i $SAMPLE_AMPLIFICATION_FILE > $SAMPLE_AMPLIFICATION_FILE'.sh'

OUTPUT_PATH=$CNV_HOME'/output/igv_output/'$EXPERIMENT_NAME'-deletions-visualizations/'
bedtools igv -path $OUTPUT_PATH -i $SAMPLE_DELETION_FILE > $SAMPLE_DELETION_FILE'.sh'
