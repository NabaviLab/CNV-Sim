#!/usr/bin/python

__author__ = 'Abdelrahman Hosny'

import sys
import os.path
import random
from operator import itemgetter

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
            chromosome, start, end, description, score, strand = line.strip().split('\t')
            exon = {"chromosome": chromosome, \
                    "start": int(start), \
                    "end": int(end), \
                    "description": description, \
                    "score": score, \
                    "strand": strand}
            exons.append(exon)
    return exons

def generateCNVMask(number_of_exons, p_amplify=0.3, p_delete=0.2):
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
        cnv_mask[choice] = random.randrange(3,11)     # random amplifications in the range [3, 10]
        number_of_amplifications -= 1
    random.shuffle(cnv_mask)
    while number_of_deletions > 0:
        choice = random.randrange(0, number_of_exons)
        while cnv_mask[choice] != 0:
            choice = random.randrange(0, number_of_exons)
        cnv_mask[choice] = -1  # random amplifications in the range [3, 10]
        number_of_deletions -= 1
    random.shuffle(cnv_mask)

    return cnv_mask

def saveCNVList(exons, cnv_mask, filename):
    '''
    Saves the updated exon regions annotated with the amplifications or deletions to be introduced.
    :param exons: a list of exon objects
    :param cnv_mask: a list of changes from the function generateCNVMask
    :return: True when the the write to file process is successful
    '''
    with open(filename, 'w') as f:
        for i in range(0, len(exons)):
            line = exons[i]["chromosome"] + '\t' \
                    + str(exons[i]["start"]) + '\t' \
                    + str(exons[i]["end"]) + '\t' \
                    + exons[i]["description"] + '\t' \
                    + exons[i]["score"] + '\t' \
                    + exons[i]["strand"] + '\t' \
                    + str(cnv_mask[i]) + '\n'
            f.write(line)
    return True

genome_file = sys.argv[1]
target_file = sys.argv[2]
simulation_name = sys.argv[3]

genome = readGenome(genome_file)
exons = readTargets(target_file)
exons = sorted(exons, key=itemgetter("start"))


modified_genome_file = os.path.dirname(os.path.realpath(__file__)) + '/../output/cnvsim_output/'
modified_target_file = os.path.dirname(os.path.realpath(__file__)) + '/../output/cnvsim_output/'
cnv_list_file = os.path.dirname(os.path.realpath(__file__)) + '/../output/cnvsim_output/' + simulation_name + '-CNVList.bed'

cnv_mask = generateCNVMask(len(exons))
saveCNVList(exons, cnv_mask, cnv_list_file)


'''
# now write results to files
with open(modified_genome_file, 'w') as f:
    f.write(header)
    line_width = 50
    chunks = [modified_genome[i:i+line_width]+'\n' for i in range(0, len(genome), line_width)]
    f.writelines(chunks)

with open(modified_target_file, 'w') as f:
    targets = []
    for ex in exons:
        target = ex["chromosome"] + '\t' \
                 + str(ex["start"]) + \
                 '\t' + str(ex["end"]) + \
                 '\t' + ex["description"] + \
                 '\t' + ex["score"] + \
                 '\t' + ex["strand"] + '\n'
        targets.append(target)
    f.writelines(targets)
    '''