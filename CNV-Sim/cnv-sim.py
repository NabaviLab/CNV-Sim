#!/usr/bin/python

__author__ = 'Abdelrahman Hosny'

import os.path
import random
import datetime
import argparse
import subprocess

def log(message):
    print '[CNV SIM {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + "] " + message

def readGenome(filename):
    '''
    Read a FASTA genome
    :param filename: genome fasta file name
    :return: a header string and a genome string
    '''
    genome = ''
    header = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
            else:
                header = line.strip()
    return header, genome

def readTargets(filename):
    '''
    Read target file in BED format
    :param filename: target bed file name
    :return: a list of targets (exons). each target is a dictionary of:
             chromosome, start, end, description, score, strand
    '''
    exons = []
    with open(filename, 'r') as tf:
        for line in tf:
            chromosome, start, end = line.strip().split('\t')
            exon = (chromosome, int(start), int(end))
            exons.append(exon)
    return exons

def getCNVMatrix(cnv_list, regions_cout=30):
    '''
    A function to randomly divide sequential exonic targets into whole CNV regions
    :param cnv_list: a target list of exons with Copy Number Variation as the fourth column
    :return: A matrix where rows represent the region index and the first column as a list of targets in this region
    '''
    regions_cout -= 1
    number_of_targets = len(cnv_list)
    comine_size = number_of_targets / regions_cout
    cnv_matrix = []

    i = 0
    while True:
        if i == len(cnv_list):
            break

        first_target = i
        i += comine_size
        if i > len(cnv_list):
            i = len(cnv_list)
        last_target = i
        cnv_matrix.append(cnv_list[first_target:last_target])

    return cnv_matrix

def generateCNVMask(number_of_exons, p_amplify, p_delete, min_variation, max_variation):
    '''
    This function generates random Copy Number Variations mask list
    :param exons: list of exons.
    :param p_amplify: percentage of amplifications introduced
    :param p_delete: percentage of deletions introduced
    :return: a list of the same length as the exons list. each list item
             indicates the variation to be added to the exon in the same position.
             Positive number: number of amplification
             Zero: no change
             -1: delete this exon
    '''

    number_of_amplifications = int(p_amplify * number_of_exons)
    number_of_deletions = int(p_delete * number_of_exons)
    cnv_mask = [0] * number_of_exons

    # generate CNV mask (a list of amplifications and deletions to be applied to the exons)
    while number_of_amplifications > 0:
        choice = random.randrange(0, number_of_exons)
        while cnv_mask[choice] != 0:
            choice = random.randrange(0, number_of_exons)
        cnv_mask[choice] = random.randrange(min_variation, max_variation)     # random amplifications in the range [min_variation, max_variation)
        number_of_amplifications -= 1
    random.shuffle(cnv_mask)
    while number_of_deletions > 0:
        choice = random.randrange(0, number_of_exons)
        while cnv_mask[choice] != 0:
            choice = random.randrange(0, number_of_exons)
        cnv_mask[choice] = -1*random.randrange(min_variation, max_variation)  # random deletions in the range [min_variation, max_variation)
        number_of_deletions -= 1
    random.shuffle(cnv_mask)

    return cnv_mask

def simulateCNV(genome, cnv_matrix, mask):
    '''
    This function creates two new genome references and target files for control and CNV
    :param genome: text string of the genome
    :param cnv_matrix: A matrix where rows represent the region index and the first column as a list of targets in this region
    :param mask: a list to indicate the variations introduced to each CNV region
    :return: control genome, control target list, modified genome and modified target list (in order)
    '''

    control_genome = []
    control_targets = []                # initialize control target list
    cnv_genome = []
    cnv_targets = []                    # initialize cnv target list

    CONTROL_APPEND_INDEX = 0  # a value to re-base the (start, end) of amplified/deleted targets in a region
    CNV_APPEND_INDEX = 0          # a value to re-base the (start, end) of amplified/deleted targets in a region

    for k, cnv_region in enumerate(cnv_matrix):

        number_of_copies = mask[k]

        if number_of_copies > 0:
            # if amplification, it will enter this loop
            for target in cnv_region:
                for _ in range(number_of_copies):
                    amplification = genome[target[1]-100:target[2]+100]
                    cnv_genome.append(amplification)
                    chromosome = target[0]
                    start = CNV_APPEND_INDEX + 100
                    end = CNV_APPEND_INDEX + len(amplification) - 100
                    cnv_targets.append((chromosome, start, end))
                    CNV_APPEND_INDEX += len(amplification)
        elif number_of_copies < 0:
            # if deletion, it will enter this loop
            for target in cnv_region:
                for _ in range(abs(number_of_copies)):
                    deletion = genome[target[1]-100:target[2]+100]
                    control_genome.append(deletion)
                    chromosome = target[0]
                    start = CONTROL_APPEND_INDEX + 100
                    end = CONTROL_APPEND_INDEX + len(deletion) - 100
                    control_targets.append((chromosome, start, end))
                    CONTROL_APPEND_INDEX += len(deletion)

    log("merging results ..")

    return ''.join(control_genome), control_targets, ''.join(cnv_genome), cnv_targets

def targetsLengths(target_list):
    '''
    A function to calcuate the number of base pairs in the target list. This value is to be used to estimate the appropriate number of reads to be passed to WESSIM
    :param target_list: A target list of the form [(chr, start, end), ...]
    :return: integet, # of bp in all targets
    '''
    count = 0
    for target in target_list:
        count += target[2] - target[1]
    return count

def prepareTargetFile(target_file):
    '''
    sort and merge targets in the target file and writes it to '.sorted.merged'
    :param target_file: target exons in BED format
    :return: sorted and merged file name
    '''
    sorted_file = target_file + ".sorted"
    merged_file = sorted_file + ".merged"

    with open(sorted_file, "w") as f:
        subprocess.call(["sort", "-k1,1", "-k2,2n", target_file], stdout=f)
    with open(merged_file, "w") as f:
        subprocess.call(["bedtools", "merge", "-i", sorted_file], stdout=f)

    return merged_file

def callWessim(genome_file, target_file, output_file, number_of_reads, read_length=100, model_file="models/ill100v4_p.gzip", number_of_threads=4):
    '''
    Calls Wessim to generate artificial reads for the targets on the reference genome
    :param genome_file: reference genome file in FASTA format
    :param target_file: target regions file in BED format
    :param output_file: output file name
    :param number_of_reads: the number of reads to generate for this simulation
    :param read_length: the read length
    :param model_file: GemSim's empirical error models used by Wessim for NGS read generation
    :param number_of_threads: the number of threads to run this simulation
    :return: null
    '''
    os.chdir("Wessim")
    subprocess.call(["python", "Wessim1.py", \
                     "-R" + genome_file, \
                     "-B" + target_file, \
                     "-n" + str(number_of_reads), \
                     "-l" + str(read_length), \
                     "-M" + model_file, \
                     "-z", \
                     "-o" + output_file, \
                     "-t" + str(number_of_threads), \
                     "-p"])
    os.chdir("..")

def mergeReads(tmp_dir, output_dir):
    '''
    merges the base reads with normal and cnv
    :return: null
    '''
    base_file_1 = os.path.join(tmp_dir, "base_1.fastq.gz")
    base_file_2 = os.path.join(tmp_dir, "base_2.fastq.gz")
    normal_file_1 = os.path.join(tmp_dir, "control_1.fastq.gz")
    normal_file_2 = os.path.join(tmp_dir, "control_2.fastq.gz")
    cnv_file_1 = os.path.join(tmp_dir, "cnv_1.fastq.gz")
    cnv_file_2 = os.path.join(tmp_dir, "cnv_2.fastq.gz")

    with open(os.path.join(output_dir, "control_1.fastq.gz"), "w") as f:
        subprocess.call(["cat", base_file_1, normal_file_1], stdout=f)
    with open(os.path.join(output_dir, "control_2.fastq.gz"), "w") as f:
        subprocess.call(["cat", base_file_2, normal_file_2], stdout=f)
    with open(os.path.join(output_dir, "cnv_1.fastq.gz"), "w") as f:
        subprocess.call(["cat", base_file_1, cnv_file_1], stdout=f)
    with open(os.path.join(output_dir, "cnv_2.fastq.gz"), "w") as f:
        subprocess.call(["cat", base_file_2, cnv_file_2], stdout=f)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("genome", \
                        help="path to the referece genome file in FASTA format ")
    parser.add_argument("target", \
                        help="path to the target regions file in BED format")

    parser.add_argument("-m", "--name",type=str, default="test", \
                        help="a name to be used for simulated results.")
    # parser.add_argument("--cnv-list", type=str, \
    #                   help="path to a CNV list file in BED format chr | start | end | variation. If not passed, it is randomly generated using --amplifications and --deletions parameters")
    parser.add_argument("-n", "--n_reads", type=int, default=10000, \
                        help="total number of reads without variations")
    parser.add_argument("-a", "--amplifications", type=float, default=0.50, \
                        help="percentage of amplifications in range [0.0: 1.0].")
    parser.add_argument("-d", "--deletions", type=float, default=0.20, \
                        help="percentage of deletions in range [0.0: 1.0].")
    parser.add_argument("-min", "--minimum", type=float, default=10, \
                        help="minimum number of amplifications/deletions introduced")
    parser.add_argument("-max", "--maximum", type=float, default=15, \
                        help="maximum number of amplifications/deletions introduced")

    args = parser.parse_args()

    genome_file = os.path.join(os.getcwd(), args.genome)
    target_file = os.path.join(os.getcwd(), args.target)
    simulation_name = args.name
    number_of_reads = args.n_reads

    output_dir = os.path.join(os.getcwd(), simulation_name)
    tmp_dir = os.path.join(os.getcwd(), "tmp")

    # create a temporary directory for intermediate files
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # create output directory for final results
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    cnv_list_file = os.path.join(output_dir, "CNVList.bed")
    control_genome_file = os.path.join(tmp_dir, "ControlGenome.fa")
    control_target_file = os.path.join(tmp_dir, "ControlTarget.bed")
    cnv_genome_file = os.path.join(tmp_dir, "CNVGenome.fa")
    cnv_target_file = os.path.join(tmp_dir, "CNVTarget.bed")

    base_reads_file = os.path.join(tmp_dir, "base")
    control_reads_file = os.path.join(tmp_dir, "control")
    cnv_reads_file = os.path.join(tmp_dir, "cnv")

    log("loading genome file ..")
    header, genome = readGenome(genome_file)
    log("successfully loaded a genome of length " + `len(genome)`)

    log("sorting and merging targets ..")
    target_file = prepareTargetFile(target_file)

    log("loading target file ..")
    targets = readTargets(target_file)

    log("generating CNV list ..")

    cnv_matrix = getCNVMatrix(targets)
    mask = generateCNVMask(len(cnv_matrix), args.amplifications, args.deletions, args.minimum, args.maximum)

    with open(cnv_list_file, 'w') as f:
        for i, cnv_region in enumerate(cnv_matrix):
            line = cnv_region[0][0] + '\t' \
                   + `cnv_region[0][1]` + '\t' \
                   + `cnv_region[-1][2]` + '\t' \
                   + `mask[i]` + '\n'
            f.write(line)
    log("CNV list saved to " + cnv_list_file)

    # call Wessim to generate reads from the genome file and the target list
    log("generating reads for the target exons ..")
    log("delegating job to Wessim ...")
    callWessim(genome_file, target_file, base_reads_file, number_of_reads)

    # generate genome and target files for the control and the CNV
    log("simulating copy number variations (amplifications/deletions)")
    control_genome, control_targets, cnv_genome, cnv_targets = simulateCNV(genome, cnv_matrix, mask)


    log("saving to the control genome file ..")
    with open(control_genome_file, 'w') as fw:
        fw.write(header + "\n")
        n = 50
        l = [control_genome[i:i + n] for i in range(0, len(control_genome), n)]
        for line in l:
            fw.write(line + "\n")

    log("saving to the control target file ..")
    with open(control_target_file, 'w') as tw:
        for target in control_targets:
            line = target[0] + "\t" + `target[1]` + "\t" + `target[2]` + "\n"
            tw.write(line)

    log("delegating job to Wessim ...")
    callWessim(control_genome_file, control_target_file, control_reads_file, int(args.deletions * number_of_reads))

    log("saving to the CNV genome file ..")
    with open(cnv_genome_file, 'w') as fw:
        fw.write(header + "\n")
        n = 50
        l = [cnv_genome[i:i + n] for i in range(0, len(cnv_genome), n)]
        for line in l:
            fw.write(line + "\n")

    log("saving to the cnv target file ..")
    with open(cnv_target_file, 'w') as tw:
        for target in cnv_targets:
            line = target[0] + "\t" + `target[1]` + "\t" + `target[2]` + "\n"
            tw.write(line)

    log("delegating job to Wessim ...")
    callWessim(cnv_genome_file, cnv_target_file, cnv_reads_file, int(args.amplifications * number_of_reads))

    log("merging results ..")
    mergeReads(tmp_dir,output_dir)

    log("cleaning temporary files")
    subprocess.call(["rm", "-rf", tmp_dir])


if __name__ == '__main__':
    main()



