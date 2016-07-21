# Copy Number Variation Simulator (CNV-Sim)
In genomics, Copy Number Variations (CNVs) is a type of structureal variation in a genome where sections of the genome are repeated. 
The number if repetitions (duplications) varies between individuals in the human population.

The Copy Number Variation Simulator (CNV Sim) is a tool used to generate a set of artificial DNA fragments for Next Generation Sequencing (NGS) read simulation.
When aligned back to the reference genome, the artificial generated reads show variations in the CNV regions. Variations can be either amplifications of deletions.

CNV-Sim offers two types of simulation:

- CNV simulation in whole genome. CNV-Sim wraps the functionality of [ART](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/) to introduce variations in the genome.
- CNV simulation in whole exome. CNV-Sim wraps the functionality of [Wessim](https://github.com/sak042/Wessim) to introduce variations in the targets.  

## Use CNV-Sim
### CNV-Sim website
Coming soon!

### Download as a standalone application
Run the setup script appropriate for your operating system to install required dependencies.

### Download as a Docker container
Coming soon!

### Use from Galaxy
Coming soon!

### Use from Bioconductor
Coming soon!

## Required Input
### whole genome
- Reference genome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format.

### whole exome
- Reference genome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format.
- Target exons file in [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format. 
The exome consists of introns and exons. The target file here should indicate the start and end positions of exons (not exomes).
 
## CNV-Sim options
```
usage: cnv-sim.py [-h] [-m NAME] [--cnv_list CNV_LIST] [-n N_READS]
                  [-g REGIONS_COUNT] [-a AMPLIFICATIONS] [-d DELETIONS]
                  [-min MINIMUM] [-max MAXIMUM]
                  {genome,exome} genome [target]

positional arguments:
  {genome,exome}        simulate copy number variations in whole genome or
                        exome regions
  genome                path to the referece genome file in FASTA format
  target                path to the target regions file in BED format (if
                        using exome) (default: None)

optional arguments:
  -h, --help            show this help message and exit
  -m NAME, --name NAME  a name to be used for simulated results. (default:
                        test)
  --cnv_list CNV_LIST   path to a CNV list file in BED format chr | start |
                        end | variation. If not passed, it is randomly
                        generated using CNV list parameters below (default:
                        None)
  -n N_READS, --n_reads N_READS
                        total number of reads without variations (default:
                        10000)

CNV list parameters:
  parameters to be used if CNV list is not passed

  -g REGIONS_COUNT, --regions_count REGIONS_COUNT
                        number of CNV regions to be randomly generated
                        (default: 30)
  -a AMPLIFICATIONS, --amplifications AMPLIFICATIONS
                        percentage of amplifications in range [0.0: 1.0].
                        (default: 0.5)
  -d DELETIONS, --deletions DELETIONS
                        percentage of deletions in range [0.0: 1.0]. (default:
                        0.2)
  -min MINIMUM, --minimum MINIMUM
                        minimum number of amplifications/deletions introduced
                        (default: 10)
  -max MAXIMUM, --maximum MAXIMUM
                        maximum number of amplifications/deletions introduced
                        (default: 15)
```