# Copy Number Variation Simulator (CNV-Sim)
In genomics, Copy Number Variations (CNVs) is a type of structural variation in a genome where sections of the genome are repeated. 
The number if repetitions (duplications) varies between individuals in the human population.

The Copy Number Variation Simulator (CNV Sim) is a tool used to generate a set of artificial DNA fragments for Next Generation Sequencing (NGS) read simulation.
When aligned back to the reference genome, the artificial generated reads show variations in the CNV regions. Variations can be either amplifications of deletions.

CNV-Sim offers two types of simulation:

- CNV simulation in whole genome. CNV-Sim wraps the functionality of [ART](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/) to introduce variations in the genome.
- CNV simulation in whole exome. CNV-Sim wraps the functionality of [Wessim](https://github.com/sak042/Wessim) to introduce variations in the targets.  

## Use CNV-Sim

### Download as a standalone application
Run the setup script appropriate for your operating system to install required dependencies.

Refer to the below [CNV-Sim options](#cnv-sim-options) section for more details on how to use all available command-line options

### Download as a Docker container
We prefer that you run CNV-Sim from the [Docker](http://www.docker.com) container as it has all the dependencies installed inside yet; no need to go through any of the setup scripts included here.

*New to Docker?* Read this [blog post](https://www.toptal.com/devops/getting-started-with-docker-simplifying-devops) to understand how it works
(yet you don't need to be a Docker expert to get CNV-Sim Docker container running!).

*Install Docker:* [Linux](https://docs.docker.com/engine/installation/#/on-linux), [Mac OS](https://docs.docker.com/docker-for-mac/) and [Windows](https://docs.docker.com/docker-for-windows/)

*How to use CNV-Sim from Docker container:*


```shell
docker run -v <absolute_local_path_to_input_directory>:/my_data nabavilab/cnv-sim ./cnv-sim.py -o /my_data/<simulation_name> [OPTIONS] {genome, exome} /my_data/<genome_file> [/my_data/<target_file>]
```

where:

- `<absolute_local_path_to_input_directory>` is the directory where the required input files exist (genome.fa and target.bed if using exome simulation); *MUST* be absolute path 
- `<simulation_name>` is used as the name of the output folder in your data directory, where the output of CNV Sim will be saved
- `<genome_file>` is the FASTA file name for the genome reference
- `<target_file>` is the BED file name for the targets (only if using exome as the simulation type)

Refer to the below [CNV-Sim options](#cnv-sim-options) section for more details on how to use all available command-line options


### Use from Galaxy

CNV-Sim is available as a [Galaxy](https://galaxyproject.org/) tool. 

Download it from the Galaxy Tool Shed: [https://toolshed.g2.bx.psu.edu/view/ahosny/cnvsim/](https://toolshed.g2.bx.psu.edu/view/ahosny/cnvsim/)


## Required Input
### whole genome
- Reference genome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format.

### whole exome
- Reference genome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format.
- Target exons file in [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format. 
The exome consists of introns and exons. The target file here should indicate the start and end positions of exons (not exomes).
 
## CNV-Sim options
```
usage: cnv-sim.py [-h] [-o OUTPUT_DIR_NAME] [-n N_READS] [-l READ_LENGTH]
                  [--cnv_list CNV_LIST] [-g REGIONS_COUNT] [-a AMPLIFICATIONS]
                  [-d DELETIONS] [-min MINIMUM] [-max MAXIMUM]
                  {genome,exome} genome [target]

positional arguments:
  {genome,exome}        simulate copy number variations in whole genome or
                        exome regions
  genome                path to the referece genome file in FASTA format
  target                path to the target regions file in BED format (if
                        using exome) (default: None)

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_DIR_NAME, --output_dir_name OUTPUT_DIR_NAME
                        a name to be used to create the output directory
                        (overrides existing directory with the same name).
                        (default: test)
  -n N_READS, --n_reads N_READS
                        total number of reads without variations (default:
                        10000)
  -l READ_LENGTH, --read_length READ_LENGTH
                        read length (bp) (default: 100)
  --cnv_list CNV_LIST   path to a CNV list file in BED format chr | start |
                        end | variation. If not passed, it is randomly
                        generated using CNV list parameters below (default:
                        None)

CNV list parameters:
  parameters to be used if CNV list is not passed

  -g REGIONS_COUNT, --regions_count REGIONS_COUNT
                        number of CNV regions to be randomly generated
                        (default: 30)
  -a AMPLIFICATIONS, --amplifications AMPLIFICATIONS
                        percentage of amplifications in range [0.0: 1.0].
                        (default: 0.3)
  -d DELETIONS, --deletions DELETIONS
                        percentage of deletions in range [0.0: 1.0]. (default:
                        0.2)
  -min MINIMUM, --minimum MINIMUM
                        minimum number of amplifications/deletions introduced
                        (default: 3)
  -max MAXIMUM, --maximum MAXIMUM
                        maximum number of amplifications/deletions introduced
                        (default: 10)
```

## License
The MIT License (MIT)

Copyright (c) 2016 NabaviLab
