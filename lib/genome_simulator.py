__author__ = 'Abdelrahman Hosny'

import random
import os
import subprocess

def generateCNVMask(mask_length, p_amplify, p_delete, min_variation, max_variation):
    '''
    This function generates random Copy Number Variations mask list
    :param exons: list of regions.
    :param p_amplify: percentage of amplifications introduced
    :param p_delete: percentage of deletions introduced
    :return: a list of the same length as the exons list. each list item
             indicates the variation to be added to the exon in the same position.
             Positive number: number of amplification
             Zero: no change
             -1: delete this exon
    '''

    number_of_amplifications = int(p_amplify * mask_length)
    number_of_deletions = int(p_delete * mask_length)
    cnv_mask = [0] * mask_length

    # generate CNV mask (a list of amplifications and deletions to be applied to the exons)
    while number_of_amplifications > 0:
        choice = random.randrange(0, mask_length)
        while cnv_mask[choice] != 0:
            choice = random.randrange(0, mask_length)
        cnv_mask[choice] = random.randrange(min_variation, max_variation)     # random amplifications in the range [min_variation, max_variation)
        number_of_amplifications -= 1
    random.shuffle(cnv_mask)
    while number_of_deletions > 0:
        choice = random.randrange(0, mask_length)
        while cnv_mask[choice] != 0:
            choice = random.randrange(0, mask_length)
        cnv_mask[choice] = -1*random.randrange(min_variation, max_variation)  # random deletions in the range [min_variation, max_variation)
        number_of_deletions -= 1
    random.shuffle(cnv_mask)

    return cnv_mask

def generateCNVList(chromosome_length, number_of_regions, p_amplify, p_delete, min_variation, max_variation):
    '''

    :param chromosome_length:
    :param number_of_regions:
    :param p_amplify:
    :param p_delete:
    :param min_variation:
    :param max_variation:
    :return:
    '''
    region_length = chromosome_length / number_of_regions
    start = 0
    end = region_length
    cnv_list = []
    mask = generateCNVMask(number_of_regions, p_amplify, p_delete, min_variation, max_variation)
    for i in range(number_of_regions):
        # jump forward start
        jump_start = start + int(random.randrange(int(0.40*region_length), int(0.45*region_length)))

        # jump backward end
        jump_end = end - int(random.randrange(int(0.40*region_length), int(0.45*region_length)))

        cnv_list.append([jump_start, jump_end, mask[i]])
        start = end
        end += region_length

    return cnv_list

def simulateCNV(genome, cnv_list):
    '''
    Simulate the control genome and the CNV genome
    :param genome: original genome sequence
    :param cnv_list: a list of region variations (chromosome, start, end, variation)
    :return: control_genome, cnv_genome
    '''
    control_genome = []
    cnv_genome = []

    for cnv in cnv_list:
        start = cnv[1]
        end = cnv[2]
        variation = cnv[3]
        sequence = genome[start: end]

        if variation > 0:
            # amplification
            amplification = genome[start-100: start] + sequence * variation + genome[end:end+100]
            cnv_genome.append(amplification)
        elif variation < 0:
            # deletion
            deletion = genome[start-100: start] + sequence * variation + genome[end:end+100]
            control_genome.append(deletion)

    return ''.join(control_genome), ''.join(cnv_genome)

def callART(genome_file, output_file, read_length=100, fold_coverage=1):
    '''

    :param genome_file:
    :param output_file:
    :param number_of_reads:
    :param read_length:
    :param number_of_threads:
    :return:
    '''
    os.chdir("lib/ART")
    subprocess.call(["./art_illumina", \
                     "-sam", \
                     "-i", genome_file, \
                     "-l", str(read_length), \
                     "-f", str(fold_coverage), \
                     "-o", output_file])
    os.chdir("../..")