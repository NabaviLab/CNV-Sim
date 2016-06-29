#!/usr/bin/python

__author__ = 'Abdelrahman Hosny'

import sys
import os.path
import random
import numpy as np
import datetime
import argparse

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
    This function updates the genome with copy number variations indicated in the mask and applied to the CNV matrix
    :param genome: text string of the genome
    :param cnv_list: a list of tuples containing (chromosome, start, end, variation)
    :param genomic_regions: True if the cnv_list passed represents genomic regions, not exonic targets
    :return: control genome, control target list, modified genome and modified target list (in order)
    '''

    control_genome = []
    control_genome.append(genome)       # initialize control genome
    control_targets = []                # initialize control target list
    cnv_genome = []
    cnv_genome.append(genome)           # initialize cnv genome
    cnv_targets = []                    # initialize cnv target list

    CNV_APPEND_INDEX = len(genome)          # a value to re-base the (start, end) of amplified/deleted targets in a region
    CONTROL_APPEND_INDEX = len(genome)      # a value to re-base the (start, end) of amplified/deleted targets in a region

    for i in range(len(cnv_matrix)):
        for j in range(len(cnv_matrix[i])):
            control_targets.append(cnv_matrix[i][j])
            cnv_targets.append(cnv_matrix[i][j])

    for k, cnv_region in enumerate(cnv_matrix):

        number_of_copies = mask[k]
        region_start = cnv_region[0][1] - 100
        region_end = cnv_region[-1][2] + 100

        region_string = genome[region_start:region_end]

        if number_of_copies > 0:
            # if amplification, it will enter this loop
            for target in cnv_region:
                for _ in range(number_of_copies):
                    amplification = genome[target[1]-100:target[2]+100]
                    cnv_genome.append(amplification)
                    chromosome = target[0]
                    start = CNV_APPEND_INDEX
                    end = CNV_APPEND_INDEX + len(amplification)
                    cnv_targets.append((chromosome, start, end))
                    CNV_APPEND_INDEX += len(amplification)
        elif number_of_copies < 0:
            # if deletion, it will enter this loop
            for target in cnv_region:
                for _ in range(abs(number_of_copies)):
                    deletion = genome[target[1]-100:target[2]+100]
                    control_genome.append(deletion)
                    chromosome = target[0]
                    start = CONTROL_APPEND_INDEX
                    end = CONTROL_APPEND_INDEX + len(deletion)
                    control_targets.append((chromosome, start, end))
                    CONTROL_APPEND_INDEX += len(deletion)

    log("merging results ..")

    return ''.join(control_genome), control_targets, ''.join(cnv_genome), cnv_targets

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("genome", \
                        help="path to the referece genome file in FASTA format ")
    parser.add_argument("target", \
                        help="path to the target regions file in BED format")
    parser.add_argument("-n", "--name",type=str, default="test", \
                        help="a name to be used for simulated results.")
    # parser.add_argument("--cnv-list", type=str, \
    #                   help="path to a CNV list file in BED format chr | start | end | variation. If not passed, it is randomly generated using --amplifications and --deletions parameters")
    parser.add_argument("-a", "--amplifications", type=float, default=0.50, \
                        help="percentage of amplifications in range [0.0: 1.0].")
    parser.add_argument("-d", "--deletions", type=float, default=0.0, \
                        help="percentage of deletions in range [0.0: 1.0].")
    parser.add_argument("-min", "--minimum", type=float, default=10, \
                        help="minimum number of amplifications/deletions introduced")
    parser.add_argument("-max", "--maximum", type=float, default=15, \
                        help="maximum number of amplifications/deletions introduced")
    parent_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
    default_output_dir = parent_dir + "/output/cnvsim_output/"
    parser.add_argument("-o", "--output_dir", type=str, default=default_output_dir , \
                        help="output directory to write results to.")
    args = parser.parse_args()

    genome_file = args.genome
    target_file = args.target
    simulation_name = args.name
    cnv_list_file = args.output_dir + simulation_name + '-CNVList.bed'
    control_genome_file = args.output_dir + simulation_name + '-ControlGenome.fa'
    control_target_file = args.output_dir + simulation_name + '-ControlTarget.bed'
    cnv_genome_file = args.output_dir + simulation_name + '-CNVGenome.fa'
    cnv_target_file = args.output_dir + simulation_name + '-CNVTarget.bed'

    log("loading genome file ..")
    header, genome = readGenome(genome_file)
    log("successfully read a genome of length " + `len(genome)`)

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

    log("simulating copy number variations using the generated control files")
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
    log("Control files saved to " + args.output_dir)

    log("saving to the cnv genome file ..")
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
    log("CNV files saved to " + args.output_dir)


if __name__ == '__main__':
    main()


