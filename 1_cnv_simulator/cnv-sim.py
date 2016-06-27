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

def simulateCNV(genome, cnv_list, genomic_regions=False):
    '''
    This function updates the genome with copy number variations indicated in the CNV list
    :param genome: text string of the genome
    :param cnv_list: a list of tuples containing (chromosome, start, end, variation)
    :param genomic_regions: True if the cnv_list passed represents genomic regions, not exonic targets
    :return: control genome, control target list, modified genome and modified target list (in order)
    '''
    log("total number of exonic targets: " + `len(cnv_list)`)
    control_genome = []
    control_targets = []
    cnv_genome = []
    cnv_targets = []

    control_genome.append(genome)       # initialize control genome
    control_targets = cnv_list[:]       # initialize control target list
    cnv_genome.append(genome)           # initialize cnv genome
    cnv_targets = cnv_list[:]           # initialize cnv target list

    control_start_pointer = len(genome) # a pointer used to track correct location of variations in the control target list
    cnv_start_pointer = len(genome)     # a pointer used to track correct location of variations in the cnv target list

    for i, cnv in enumerate(cnv_list):
        chromosome, start, end, number_of_copies = cnv[0], cnv[1], cnv[2], cnv[3]
        if genomic_regions is True:
            start -= 100                    # add 100 bp to the left to avoid alignment problems
            end += 100                      # add 100 bp to the right to avoid alignment problems
        exon_string = genome[start:end]

        if number_of_copies > 0:
            # amplification
            # add it to the CNV genome and target list
            amplification = exon_string * number_of_copies
            amplification_start = cnv_start_pointer
            amplification_end = cnv_start_pointer + len(amplification)

            cnv_genome.append(amplification)
            cnv_targets.append((chromosome, amplification_start, amplification_end, number_of_copies))

            cnv_start_pointer = amplification_end
        else:
            # deletions
            # amplify control and not the cnv
            deletion = exon_string * abs(number_of_copies)
            deletion_start = control_start_pointer
            deletion_end = control_start_pointer + len(deletion)

            control_genome.append(deletion)
            control_targets.append((chromosome, deletion_start, deletion_end, number_of_copies))

            control_start_pointer = deletion_end

        percentage = i / float(len(cnv_list)) * 100.0
        if int(i) % 20 == 0:
            log("simulation progress: " + `percentage` + "% ..")
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
                        help="percentage of amplifications in rage [0.0: 1.0].")
    parser.add_argument("-d", "--deletions", type=float, default=0.20, \
                        help="percentage of deletions in range [0.0: 1.0].")
    parser.add_argument("-min", "--minimum", type=float, default=3, \
                        help="minimum number of amplifications/deletions introduced")
    parser.add_argument("-max", "--maximum", type=float, default=10, \
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

    log("reading genome file ..")
    header, genome = readGenome(genome_file)
    log("successfully read a genome of length " + `len(genome)`)

    log(" reading target file ..")
    exons = readTargets(target_file)
    log("successfully loaded " + str(len(exons)) + " exonic regions")

    log("generating CNV list ..")
    cnv_mask = generateCNVMask(len(exons), args.amplifications, args.deletions, args.minimum, args.maximum)
    cnv_list = saveCNVList(exons, cnv_mask, cnv_list_file)
    log("CNV list saved to " + cnv_list_file)

    log("simulating copy number variations using the generated control files")
    control_genome, control_targets, cnv_genome, cnv_targets = simulateCNV(genome, cnv_list)

    log("saving to the control genome file ..")
    with open(control_genome_file, 'w') as f:
        f.write(header + "\n")
        f.write(control_genome)

    log("saving to the control target file ..")
    with open(control_target_file, 'w') as tw:
        for target in control_targets:
            line = target[0] + "\t" + `target[1]` + "\t" + `target[2]` + "\n"
            tw.write(line)
    log("Control files saved to " + args.output_dir)

    log("saving to the cnv genome file ..")
    with open(cnv_genome_file, 'w') as fw:
        fw.write(header + "\n")
        fw.write(cnv_genome)

    log("saving to the cnv target file ..")
    with open(cnv_target_file, 'w') as tw:
        for target in cnv_targets:
            line = target[0] + "\t" + `target[1]` + "\t" + `target[2]` + "\n"
            tw.write(line)
    log("CNV files saved to " + args.output_dir)

if __name__ == '__main__':
    main()