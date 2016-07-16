#!/usr/bin/python

__author__ = 'Abdelrahman Hosny'

import random
import subprocess
import os

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

    return ''.join(control_genome), control_targets, ''.join(cnv_genome), cnv_targets

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
    os.chdir("lib/Wessim")
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
    os.chdir("../..")