#!/usr/bin/python

__author__ = 'Abdelrahman Hosny'

import sys
import os.path
import random

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

def generateCNVMask(number_of_exons, p_amplify=0.5, p_delete=0.01):
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
        cnv_mask[choice] = -1*random.randrange(3,11)  # random deletions in the range [3, 10]
        number_of_deletions -= 1
    random.shuffle(cnv_mask)

    return cnv_mask

def saveCNVList(exons, cnv_mask, filename):
    '''
    Saves and returns the updated exon regions annotated with the amplifications or deletions to be introduced.
    :param exons: a list of exon objects
    :param cnv_mask: a list of changes from the function generateCNVMask
    :return: a list of tuples containing (chromosome, start, end, variation)
    '''
    cnv_list = []
    with open(filename, 'w') as f:
        for i in range(0, len(exons)):
            line = exons[i]["chromosome"] + '\t' \
                    + `exons[i]["start"]` + '\t' \
                    + `exons[i]["end"]` + '\t' \
                    + `cnv_mask[i]` + '\n'
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
    print "total number of targets: " + `len(cnv_list)`
    cnv_genome_list = []
    cnv_targets = []
    fragment_start = 0      # a value used to determine the start of the unchanged fragments
    ADJUST = 0              # a value used to adjust the start and end positions of all targets
    for i, cnv in enumerate(cnv_list):
        start, end, number_of_copies = cnv[1], cnv[2], cnv[3]
        exon_string = genome[start:end]

        # add unchanged fragment
        cnv_genome_list.append(genome[fragment_start: start])

        if number_of_copies < 0:    # deletion
            # don't add any exons
            # modify the global ADJUST value so that the next iteration we capture the correct exon
            ADJUST -= len(exon_string)

        elif number_of_copies >= 0:        # amplifications
            # modify the genome
            amplification = exon_string * (number_of_copies + 1)        # to account for the original copy (if number_of_copies = 0)
            cnv_genome_list.append(amplification)

            # modify the target list
            new_exon_start = ADJUST + start
            new_exon_end = ADJUST + end + (len(exon_string) * number_of_copies)

            exon = (cnv[0], new_exon_start, new_exon_end, number_of_copies)
            cnv_targets.append(exon)

            # modify the global ADJUST value so that the next iteration we capture the correct exon
            ADJUST += len(exon_string) * number_of_copies

        fragment_start = end  # the next fragment start will be from the end of that cut

        percentage = i / float(len(cnv_list)) * 100.0
        if  int(i) % 5 == 0:
            print "simulated copy number variations: " + `percentage` + "% .."
    print "merging results .."
    cnv_genome = ''.join(cnv_genome_list)
    return cnv_genome, cnv_targets

def simulateControl(genome, cnv_list):
    '''
    This function amplifies one-copy targets that have one or more deletions marked in the CNV list
    :param genome: text string of the original genome
    :param cnv_list: a list of tuples containing (chromosome, start, end, variation)
    :return: modified genome and modified cnv_list
    '''
    print "marking targets that will be deleted for amplification first .."
    marked_list = []
    for i in range(len(cnv_list)):
        if cnv_list[i][3] < 0:      # deletions
            new_target = (cnv_list[i][0], \
                          cnv_list[i][1], \
                          cnv_list[i][2], \
                          (abs(cnv_list[i][3]) - 1))  # provide many amplifications that will count for deletions after that
        else:                       # no change/amplifications
            new_target = (cnv_list[i][0], \
                          cnv_list[i][1], \
                          cnv_list[i][2], \
                          0)
        marked_list.append(new_target)

    print "simulating the added deletions"
    cnv_genome, cnv_targets = simulateCNV(genome, marked_list)
    return cnv_genome, cnv_targets

def controlListToCNVList(original_cnv_list, control_cnv_list):
    '''
    This function updates the original CNV list with the new coordinates after to-be-deleted targets were amplified
    :param original_cnv_list: a list of tuples containing (chromosome, start, end, variation) generated by saveCNVList
    :param control_cnv_list: a list of tuples containing (chromosome, start, end, variation) generated by simulateControl
    :return: a modified CNV List that reflects the new target regions
    '''
    if len(original_cnv_list) != len(control_cnv_list):
        print 'the original cnv list and the control cnv list have different lengths!'
        print 'aborting now ..'
        exit(0)
    print 'mapping new coordinates from the control cnv list to the original cnv list'
    new_list = []
    for i in range(len(control_cnv_list)):
        if control_cnv_list[i][3] > 0:
            new_target = (original_cnv_list[i][0], \
                          control_cnv_list[i][1], \
                          control_cnv_list[i][2], \
                          -1)       # mark for deletion if it was amplified in the control cnv list
        else:
            new_target = (original_cnv_list[i][0], \
                          control_cnv_list[i][1], \
                          control_cnv_list[i][2], \
                          original_cnv_list[i][3])
        new_list.append(new_target)

    return new_list



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

control_genome, control_cnv_list = simulateControl(genome, cnv_list)

print "generating the control genome file .."
with open(control_genome_file, 'w') as f:
    f.write(header + "\n")
    f.write(control_genome)

print "generating the control target file .."
with open(control_target_file, 'w') as tw:
    for target in control_cnv_list:
        line = target[0] + "\t" + `target[1]` + "\t" + `target[2]` + "\n"
        tw.write(line)

print "Control files saved to " + os.path.dirname(os.path.realpath(__file__)) + '/../output/cnvsim_output/'

cnv_list = controlListToCNVList(cnv_list, control_cnv_list)

print "simulating copy number variations using the generated control files"
cnv_genome, cnv_target_list = simulateCNV(control_genome, cnv_list)

print "generating the cnv genome file .."
with open(cnv_genome_file, 'w') as fw:
    fw.write(header + "\n")
    fw.write(cnv_genome)
print "generating the cnv target file .."
with open(cnv_target_file, 'w') as tw:
    for target in cnv_target_list:
        line = target[0] + "\t" + `target[1]` + "\t" + `target[2]` + "\n"
        tw.write(line)
print "Control files saved to " + os.path.dirname(os.path.realpath(__file__)) + '/../output/cnvsim_output/'

print len(control_genome), len(control_cnv_list)
print len(cnv_genome), len(cnv_target_list)