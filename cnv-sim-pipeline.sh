#!/bin/bash

CNV_HOME=$(pwd)
EXPERIMENT_NAME='my_test'

GENOME_FILE=$CNV_HOME'/input/test/chr1.fa'
TARGET_FILE=$CNV_HOME'/input/test/chr1-target.bed'

# Step 1 (CNV SIM)
sort -k1,1 -k2,2n $TARGET_FILE > $TARGET_FILE'.sorted'
bedtools merge -i $TARGET_FILE'.sorted' > $TARGET_FILE'.sorted.merged'
echo 'sorted and merged exons target'

cd $CNV_HOME
cd 1.\ cnv_simulator
python cnv-sim.py $GENOME_FILE $TARGET_FILE'.sorted.merged' $EXPERIMENT_NAME

cd $CNV_HOME
CONTROL_GENOME=$CNV_HOME'/output/cnvsim_output/'$EXPERIMENT_NAME'-ControlGenome.fa'
CONTROL_TARGET=$CNV_HOME'/output/cnvsim_output/'$EXPERIMENT_NAME'-ControlTarget.bed'
sh wessim-align-sort-vis.sh $CONTROL_GENOME $CONTROL_TARGET 'normal'

CONTROL_GENOME=$CNV_HOME'/output/cnvsim_output/'$EXPERIMENT_NAME'-CNVGenome.fa'
CONTROL_TARGET=$CNV_HOME'/output/cnvsim_output/'$EXPERIMENT_NAME'-CNVTarget.bed'
sh wessim-align-sort-vis.sh $CONTROL_GENOME $CONTROL_TARGET 'cnv'