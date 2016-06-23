#!/usr/bin/python

__author__ = 'Abdelrahman Hosny'

import sys
import os.path
import random
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
            exon = {"chromosome": chromosome, \
                    "start": int(start), \
                    "end": int(end)}
            exons.append(exon)
    return exons

def generateCNVMask(number_of_exons, p_amplify, p_delete):
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
    log("total number of targets: " + `len(cnv_list)`)
    cnv_genome_list = []
    cnv_targets = []
    fragment_start = 0      # a value used to determine the start of the unchanged fragments
    ADJUST = 0              # a value used to adjust the start and end positions of all targets
    for i, cnv in enumerate(cnv_list):
        start, end, number_of_copies = cnv[1], cnv[2], cnv[3]
        exon_string = genome[start:end]

        # add unchanged fragment
        cnv_genome_list.append(genome[fragment_start: start])

        if number_of_copies < 0:    # deletions
            # keep only one copy
            # modify genome
            deletion_cut_at = len(exon_string) / abs(number_of_copies)
            deletion_fragment = exon_string[:deletion_cut_at]           # to keep only one copy of the deletion
            cnv_genome_list.append(deletion_fragment)

            # modify the target list
            new_exon_start = ADJUST + start
            new_exon_end = ADJUST + start + deletion_cut_at
            exon = (cnv[0], new_exon_start, new_exon_end, number_of_copies)
            cnv_targets.append(exon)


            # modify the global ADJUST value so that the next iteration we capture the correct exon
            ADJUST += deletion_cut_at

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
        if  int(i) % 20 == 0:
            log("simulation progress: " + `percentage` + "% ..")
    log("merging results ..")
    cnv_genome = ''.join(cnv_genome_list)
    return cnv_genome, cnv_targets

def simulateControl(genome, cnv_list):
    '''
    This function amplifies one-copy targets that have one or more deletions marked in the CNV list
    :param genome: text string of the original genome
    :param cnv_list: a list of tuples containing (chromosome, start, end, variation)
    :return: modified genome and modified cnv_list
    '''
    log("marking targets that will be deleted for amplification first ..")
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

    log("simulating the added deletions")
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
        log('the original cnv list and the control cnv list have different lengths!')
        log('aborting now ..')
        exit(0)
    log('mapping new coordinates from the control cnv list to the original cnv list')
    new_list = []
    for i in range(len(control_cnv_list)):
        new_target = (original_cnv_list[i][0], \
                      control_cnv_list[i][1], \
                      control_cnv_list[i][2], \
                      original_cnv_list[i][3])
        new_list.append(new_target)

    return new_list

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
                        help="percentage of amplifications in rage [0.0: 1.0].")
    parser.add_argument("-d", "--deletions", type=float, default=0.20, \
                        help="percentage of deletions in range [0.0: 1.0].")
    parent_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
    default_output_dir=parent_dir + "/output/cnvsim_output/"
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

    log("reading genome file ..")
    header, genome = readGenome(genome_file)
    log("successfully read a genome of length " + `len(genome)`)

    log(" reading target file ..")
    exons = readTargets(target_file)
    log("successfully loaded " + str(len(exons)) + " exonic regions")

    log("generating CNV list ..")
    cnv_mask = generateCNVMask(len(exons), args.amplifications, args.deletions)
    cnv_list = saveCNVList(exons, cnv_mask, cnv_list_file)
    log("CNV list saved to " + cnv_list_file)

    control_genome, control_cnv_list = simulateControl(genome, cnv_list)

    log("generating the control genome file ..")
    with open(control_genome_file, 'w') as f:
        f.write(header + "\n")
        f.write(control_genome)

    log("generating the control target file ..")
    with open(control_target_file, 'w') as tw:
        for target in control_cnv_list:
            line = target[0] + "\t" + `target[1]` + "\t" + `target[2]` + "\n"
            tw.write(line)
    log("Control files saved to " + args.output_dir)

    cnv_list = controlListToCNVList(cnv_list, control_cnv_list)

    log("simulating copy number variations using the generated control files")
    cnv_genome, cnv_target_list = simulateCNV(control_genome, cnv_list)

    log("generating the cnv genome file ..")
    with open(cnv_genome_file, 'w') as fw:
        fw.write(header + "\n")
        fw.write(cnv_genome)

    log("generating the cnv target file ..")
    with open(cnv_target_file, 'w') as tw:
        for target in cnv_target_list:
            line = target[0] + "\t" + `target[1]` + "\t" + `target[2]` + "\n"
            tw.write(line)
    log("CNV files saved to " + args.output_dir)

if __name__ == '__main__':
    main()