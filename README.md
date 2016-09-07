# Copy Number Variation Simulator (CNV-Sim)
In genomics, Copy Number Variations (CNVs) is a type of structural variation in a genome where sections of the genome are duplicated or deleted. 
The number of variations (duplications/deletions) varies between individuals in the human population.

The Copy Number Variation Simulator (CNV-Sim) is a simulation tool that extends the functionality of existing next-generation sequencing read simulators 
to introduce copy number variations in the generated reads. The resulting reads encompass amplifications as well as deletions according to a predefined list of variant regions.

CNV-Sim aids testing and benchmarking tools for copy number variation detection and analysis. The tool offers two types of simulation:

- CNV simulation in whole genome. CNV-Sim utilizes the functionality of [ART](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/) to introduce variations in the genome.
- CNV simulation in targeted exome. CNV-Sim utilizes the functionality of [Wessim](https://github.com/sak042/Wessim) to introduce variations in the targets.  

## Use CNV-Sim

### Download as a standalone application
Run the setup script appropriate for your operating system to install required dependencies. 

*Alternatively*, install the following dependencies on your machine (the setup script does all that for you):

- [Python 2.7](https://www.python.org/downloads/)
- [Python pip](https://pypi.python.org/pypi/pip)
- [Pysam Module](https://github.com/pysam-developers/pysam)
- [Biopython Module](http://biopython.org/)
- [bedtools](http://bedtools.readthedocs.io/en/latest/)
- Download [ART simulator](http://www.niehs.nih.gov/research/resources/software/art) and put the binary file `art_illumina` in `cnvsim/ART` directory.

Refer to the below [CNV-Sim options](#cnv-sim-options) section for more details on how to use all available command-line options


### Download as a Docker image
We prefer that you run CNV-Sim from the [Docker](http://www.docker.com) container as it has all the dependencies installed inside yet; no need to go through any of the setup scripts included here.

**New to Docker?** Read this [blog post](https://www.toptal.com/devops/getting-started-with-docker-simplifying-devops) to understand how it works
(yet you don't need to be a Docker expert to get CNV-Sim Docker container running!).

**Install Docker:** [Linux](https://docs.docker.com/engine/installation/#/on-linux), [Mac OS](https://docs.docker.com/docker-for-mac/) and [Windows](https://docs.docker.com/docker-for-windows/)

**How to use CNV-Sim Docker image:**


```shell
docker run -v <absolute_local_path_to_input_directory>:/my_data nabavilab/cnv-sim ./cnv-sim.py -o /my_data/<simulation_name> [OPTIONS] {genome, exome} /my_data/<genome_file> [/my_data/<target_file>]
```

where:

- `<absolute_local_path_to_input_directory>` is the directory where the required input files exist (genome.fa and target.bed if using exome simulation); *MUST* be absolute path 
- `<simulation_name>` is used as the name of the output folder in your data directory, where the output of CNV Sim will be saved
- `<genome_file>` is the FASTA file name for the genome reference
- `<target_file>` is the BED file name for the targets (only if using exome as the simulation type)

Refer to the below [CNV-Sim options](#cnv-sim-options) section for more details on how to use all available command-line options or 
simply type: ```docker run nabavilab/cnv-sim ./cnv-sim -h``` to read the help.


### Use from Galaxy

CNV-Sim is available as a [Galaxy](https://galaxyproject.org/) tool. 

Download it from the Galaxy Tool Shed: [https://toolshed.g2.bx.psu.edu/view/ahosny/cnvsim/](https://toolshed.g2.bx.psu.edu/view/ahosny/cnvsim/)


## Required Input
**IMPORTANT:** for the current version of CNV-Sim, the underlying *pysam* library used Wessim might not be able to handle large reference files (such as the whole human genome) when introducing amplifications. To avoid errors, simulate CNV in one chromosome at a time (reference and target for one chromosome in a single run). We are working on handling large files seamlessly and will be released in the next version.

### whole genome
- Reference genome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format.

### targeted exome
- Reference genome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format.
- Target exons file in [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format. 
The exome consists of introns and exons. The target file here should indicate the start and end positions of exons (not exomes).
 
## CNV-Sim options
```
Usage: cnv-sim.py [-h] [-v] [-o OUTPUT_DIR_NAME] [-n N_READS] [-l READ_LENGTH]
                  [--cnv_list CNV_LIST] [-g REGIONS_COUNT]
                  [-r_min REGION_MINIMUM_LENGTH]
                  [-r_max REGION_MAXIMUM_LENGTH] [-a AMPLIFICATIONS]
                  [-d DELETIONS] [-cn_min COPY_NUMBER_MINIMUM]
                  [-cn_max COPY_NUMBER_MAXIMUM]
                  {genome,exome} genome [target]

Generates NGS short reads that encompass copy number variations in whole
genome and targeted exome sequencing

Positional arguments:
  {genome,exome}        simulate copy number variations in whole genome or
                        exome regions
  genome                path to the referece genome file in FASTA format
  target                path to the target regions file in BED format (if
                        using exome)

Optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Show program's version number and exit.
  -o OUTPUT_DIR_NAME, --output_dir_name OUTPUT_DIR_NAME
                        a name to be used to create the output directory
                        (overrides existing directory with the same name).
  -n N_READS, --n_reads N_READS
                        total number of reads without variations
  -l READ_LENGTH, --read_length READ_LENGTH
                        read length (bp)
  --cnv_list CNV_LIST   path to a CNV list file in BED format chr | start |
                        end | variation. If not passed, it is randomly
                        generated using CNV list parameters below

CNV list parameters:
  parameters to be used if CNV list is not passed

  -g REGIONS_COUNT, --regions_count REGIONS_COUNT
                        number of CNV regions to be generated randomly
  -r_min REGION_MINIMUM_LENGTH, --region_minimum_length REGION_MINIMUM_LENGTH
                        minimum length of each CNV region
  -r_max REGION_MAXIMUM_LENGTH, --region_maximum_length REGION_MAXIMUM_LENGTH
                        maximum length of each CNV region
  -a AMPLIFICATIONS, --amplifications AMPLIFICATIONS
                        percentage of amplifications in range [0.0: 1.0].
  -d DELETIONS, --deletions DELETIONS
                        percentage of deletions in range [0.0: 1.0].
  -cn_min COPY_NUMBER_MINIMUM, --copy_number_minimum COPY_NUMBER_MINIMUM
                        minimum level of variations (copy number) introduced
  -cn_max COPY_NUMBER_MAXIMUM, --copy_number_maximum COPY_NUMBER_MAXIMUM
                        maximum level of variation (copy number) introduced
```

## Contributor
Abdelrahman Hosny (@abdelrahmanhosny)

## License
The MIT License (MIT)

Copyright (c) 2016 [NabaviLab](https://nabavilab.github.io/)

