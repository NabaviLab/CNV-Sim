#!/usr/bin/python

__author__ = 'Abdelrahman Hosny'

import sys
import os.path
import random
from shutil import copyfile

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
            exon = {"chromosome": chromosome, \
                    "start": int(start), \
                    "end": int(end)}
            exons.append(exon)
    return exons

def generateCNVMask(number_of_exons, p_amplify=0.5, p_delete=0.0):
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
    cnv_list = []
    with open(filename, 'w') as f:
        for i in range(0, len(exons)):
            line = exons[i]["chromosome"] + '\t' \
                    + str(exons[i]["start"]) + '\t' \
                    + str(exons[i]["end"]) + '\t' \
                    + str(cnv_mask[i]) + '\n'
            f.write(line)
            cnv_instance = (exons[i]["chromosome"], exons[i]["start"], exons[i]["end"], cnv_mask[i])
            cnv_list.append(cnv_instance)
    return cnv_list

def simulateCNV(genome, cnv_list):
    '''
    This function updates the genome with copy number variations indicated in the CNV list
    :param genome: text string of the genome
    :param cnv_list: a list of tuples containing (chromosome, start, end, variation)
    :return: modified genome and modified target list
    '''
    ADJUST = 0          # a value used to adjust the start and end positions of all targets
    cnv_genome = genome[:]
    cnv_targets = []
    for cnv in cnv_list:
        start, end = ADJUST + cnv[1], ADJUST + cnv[2]
        number_of_copies = cnv[3]
        exon_string = cnv_genome[start:end]

        if number_of_copies == -1:        # 1 deletion (later on, we will implement multiple deletion)
            # modify the genome
            cnv_genome = cnv_genome[:start] + cnv_genome[end:]

            # modify the global ADJUST value so that the next iteration we capture the correct exon
            ADJUST -= len(exon_string)

        elif number_of_copies >= 0:        # amplifications
            # modify the genome
            for _ in range(number_of_copies):
                cnv_genome = cnv_genome[:end] + exon_string + cnv_genome[end:]

            # modify the target list
            new_exon_start = start
            new_exon_end = end + (len(exon_string) * number_of_copies)
            exon = (cnv[0], new_exon_start, new_exon_end)
            cnv_targets.append(exon)

            # modify the global ADJUST value so that the next iteration we capture the correct exon
            ADJUST += len(exon_string) * number_of_copies

    return cnv_genome, cnv_targets

# Required Input
genome_file = sys.argv[1]
target_file = sys.argv[2]
simulation_name = sys.argv[3]

# Output by CNV Simulator
cnv_list_file = os.path.dirname(os.path.realpath(__file__)) + '/../output/cnvsim_output/' + simulation_name + '-CNVList.bed'
control_genome_file = os.path.dirname(os.path.realpath(__file__)) + '/../output/cnvsim_output/' + simulation_name + '-ControlGenome.fa'
control_target_file = os.path.dirname(os.path.realpath(__file__)) + '/../output/cnvsim_output/' + simulation_name + '-ControlTarget.bed'
cnv_genome_file = os.path.dirname(os.path.realpath(__file__)) + '/../output/cnvsim_output/' + simulation_name + '-CNVGenome.fa'
cnv_target_file = os.path.dirname(os.path.realpath(__file__)) + '/../output/cnvsim_output/' + simulation_name + '-CNVTarget.bed'

print "reading genome file .."
header, genome = readGenome(genome_file)
print "successfully read a genome of length " + str(len(genome))

print "reading target file .."
exons = readTargets(target_file)
print "successfully loaded " + str(len(exons)) + " exonic regions"

print "generating CNV list .."
cnv_mask = generateCNVMask(len(exons))
cnv_list = saveCNVList(exons, cnv_mask, cnv_list_file)
print "CNV list saved to " + cnv_list_file

print "generating the control genome file .."
copyfile(genome_file, control_genome_file)
print "generating the control target file .."
copyfile(target_file, control_target_file)
print "Control files saved to " + os.path.dirname(os.path.realpath(__file__)) + '/../output/cnvsim_output/'

print "simulating copy number variations using the generated control files"
cnv_genome, cnv_target_list = simulateCNV(genome, cnv_list)

print "generating the cnv genome file .."
with open(cnv_genome_file, 'w') as fw:
    fw.write(header + "\n")
    fw.write(cnv_genome)
print "generating the cnv target file .."
with open(cnv_target_file, 'w') as tw:
    for target in cnv_target_list:
        line = target[0] + "\t" + str(target[1]) + "\t" + str(target[2]) + "\n"
        tw.write(line)
print "Control files saved to " + os.path.dirname(os.path.realpath(__file__)) + '/../output/cnvsim_output/'