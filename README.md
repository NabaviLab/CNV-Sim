# Copy Number Variation Simulator
The Copy Number Variation Simulator (CNV SIM) is a tool used to simulate random copy number variations (amplifications/deletions) in whole exome sequencing.

_more detailed description will be available soon_

The pipeline goes as follows:

vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

### Step 1:

1. Inputs: ```control genome``` *and* ```target file```

 --> **CNV SIM** -->
2. Outputs:

 2.1 ```control genome``` *and* ```target file```

 2.2 ```simulated genome``` *and* ```simulated target file```

## Simulate the control genome and target
## Step 2:
3. Inputs: ```control genome``` *and* ```target file```

 -->> **WESSIM** -->

4. Output: short reads in ```fastq``` format

vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

## Step 3:
5. Inputs: short reads in ```fastq``` format

 --> **BWA Aligner** -->

6. Output: aligned reads in ```sam``` format to the ```control genome``` passed in the first step

vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

## Step 4:
7. Input: aligned reads in ```sam``` format

 --> **Samtools** -->
8. Outputs: ```Sorted aligned reads``` *and* ```sam index```

vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

## Step 5:
9. Inputs: ```sorted sam or bam``` *and* ```sam/bam index```

 --> **IGB Browser** -->
10. Output: visualization of the reads to the ```control genome``` passed in the first step.


## Simulate the simulated control genome and simulated target
Repeat steps **2** through **5** and compare the visualization for the two ..

# How it works?

All you have to do is to use ```cnv-sim-pipeline.sh``` with proper parameters.

details to be available soon