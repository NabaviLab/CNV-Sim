#!/bin/bash

CNV_HOME=$(pwd)
EXPERIMENT_NAME='dev_test_adjusted'

ORIGINAL_GENOME_FILE=$CNV_HOME/input/test/chr1.fa
ORIGINAL_TARGET_FILE=$CNV_HOME/input/test/chr1-target.bed

echo '---------------------------------------------------------------------------------------------------------'
echo '-------------------------------- Copy Number Variation Simulator v0.9 -----------------------------------'
echo '---------------------------------------------------------------------------------------------------------'
echo ''

# Step 1 (CNV SIM)
NOW=$(date +"%Y-%m-%d %H:%M:%S")
echo "[CNV SIM $NOW] sorting target file .."
sort -k1,1 -k2,2n $ORIGINAL_TARGET_FILE > $ORIGINAL_TARGET_FILE'.sorted'
NOW=$(date +"%Y-%m-%d %H:%M:%S")
echo "[CNV SIM $NOW] merging target overlaps to obrain larger regions .."
bedtools merge -i $ORIGINAL_TARGET_FILE'.sorted' > $ORIGINAL_TARGET_FILE'.sorted.merged'
NOW=$(date +"%Y-%m-%d %H:%M:%S")
echo "[CNV SIM $NOW] sorted and merged exons target."

cd $CNV_HOME
cd 1_cnv_simulator
python cnv-sim.py -n $EXPERIMENT_NAME $ORIGINAL_GENOME_FILE $ORIGINAL_TARGET_FILE'.sorted.merged'
CONTROL_TARGETS_LENGTH=$(cat $CNV_HOME/output/cnvsim_output/$EXPERIMENT_NAME-lengths.tmp | awk 'NR==1')
CNV_TARGETS_LENGTH=$(cat $CNV_HOME/output/cnvsim_output/$EXPERIMENT_NAME-lengths.tmp | awk 'NR==2')
rm $CNV_HOME/output/cnvsim_output/$EXPERIMENT_NAME-lengths.tmp
READ_LENGTH=100
CONTROL_NUMBER_OF_READS=1000000
COVERAGE=$(($CONTROL_NUMBER_OF_READS * $READ_LENGTH / $CONTROL_TARGETS_LENGTH ))
CNV_NUMBER_OF_READS=$(( $COVERAGE * $CNV_TARGETS_LENGTH / $READ_LENGTH))
echo "[CNV SIM $NOW] adjusted number of reads for WESSIM. Control: $CONTROL_NUMBER_OF_READS, CNV: $CNV_NUMBER_OF_READS"
echo "[CNV SIM $NOW] transferring job to WESSIM .."
echo ""

cd $CNV_HOME
GENOME=$CNV_HOME/output/cnvsim_output/$EXPERIMENT_NAME-ControlGenome.fa
TARGET=$CNV_HOME/output/cnvsim_output/$EXPERIMENT_NAME-ControlTarget.bed
sh wessim-align-sort.sh $ORIGINAL_GENOME_FILE $GENOME $TARGET $CONTROL_NUMBER_OF_READS $EXPERIMENT_NAME.normal

GENOME=$CNV_HOME/output/cnvsim_output/$EXPERIMENT_NAME-CNVGenome.fa
TARGET=$CNV_HOME/output/cnvsim_output/$EXPERIMENT_NAME-CNVTarget.bed
sh wessim-align-sort.sh $ORIGINAL_GENOME_FILE $GENOME $TARGET $CNV_NUMBER_OF_READS $EXPERIMENT_NAME.cnv

# Step 5 (IGV)
NOW=$(date +"%Y-%m-%d %H:%M:%S")
echo "[CNV SIM $NOW] generating IGV visualization scripts .."
CNV_FILE=$CNV_HOME/output/cnvsim_output/$EXPERIMENT_NAME-CNVList.bed
mkdir $CNV_HOME/output/igv_output/$EXPERIMENT_NAME
SAMPLE_AMPLIFICATION_FILE=$CNV_HOME/output/igv_output/$EXPERIMENT_NAME/$EXPERIMENT_NAME-CNVList.amplifications.bed
SAMPLE_DELETION_FILE=$CNV_HOME/output/igv_output/$EXPERIMENT_NAME/$EXPERIMENT_NAME-CNVList.deletions.bed
awk '($4 > 0 && ($3 - $2) > 1000)' $CNV_FILE | head -20 > $SAMPLE_AMPLIFICATION_FILE
awk '($4 < 0 && ($3 - $2) > 1000)' $CNV_FILE | tail -20 > $SAMPLE_DELETION_FILE

OUTPUT_PATH=$CNV_HOME/output/igv_output/$EXPERIMENT_NAME/$EXPERIMENT_NAME-amplifications-visualizations/
mkdir $OUTPUT_PATH
bedtools igv -path $OUTPUT_PATH -i $SAMPLE_AMPLIFICATION_FILE > $SAMPLE_AMPLIFICATION_FILE'.sh'

OUTPUT_PATH=$CNV_HOME/output/igv_output/$EXPERIMENT_NAME/$EXPERIMENT_NAME-deletions-visualizations/
mkdir $OUTPUT_PATH
bedtools igv -path $OUTPUT_PATH -i $SAMPLE_DELETION_FILE > $SAMPLE_DELETION_FILE'.sh'

rm $SAMPLE_AMPLIFICATION_FILE
rm $SAMPLE_DELETION_FILE
echo ''
NOW=$(date +"%Y-%m-%d %H:%M:%S")
echo "[CNV SIM $NOW] visualization scripts saved .."
echo '---------------------------------------------------------------------------------------------------------'
echo '---------------------------------------- Visualize on IGV -----------------------------------------------'
echo '---------------------------------------------------------------------------------------------------------'
echo '1. Run IGV browser'
echo '2. Load the normal sample reads from output/samtools'
echo '3. Load the CNV sample reads from output/samtools'
echo '4. From IGV Tools menu, select Run Batch Script ...'
echo '5. Select the batch script '$SAMPLE_AMPLIFICATION_FILE'.sh to generate images that show amplified targets'
echo '6. From IGV Tools menu, select Run Batch Script ...'
echo '7. Select the batch script '$SAMPLE_DELETION_FILE'.sh to generate images that show deleted targets'
echo '---------------------------------------------------------------------------------------------------------'