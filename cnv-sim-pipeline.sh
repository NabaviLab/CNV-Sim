#!/bin/bash

CNV_HOME=$(pwd)
EXPERIMENT_NAME='my_test_del'

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
sh wessim-align-sort-vis.sh $GENOME_FILE $GENOME $TARGET $EXPERIMENT_NAME'.normal'

GENOME=$CNV_HOME'/output/cnvsim_output/'$EXPERIMENT_NAME'-CNVGenome.fa'
TARGET=$CNV_HOME'/output/cnvsim_output/'$EXPERIMENT_NAME'-CNVTarget.bed'
sh wessim-align-sort-vis.sh $GENOME_FILE $GENOME $TARGET $EXPERIMENT_NAME'.cnv'