#!/usr/bin/python

__author__ = 'Abdelrahman Hosny'

import random
import subprocess
import os
import sys
import datetime
import shutil

from . import fileio


def _log(message):
    print '[CNV SIM {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + "] " + message


def getScriptPath():
    return os.path.dirname(os.path.realpath(__file__))


def _generateCNVMatrix(targets_list, regions_count, region_min_length, region_max_length):
    '''
    A function to randomly divide sequential exonic targets into whole CNV regions
    :param targets_list: a list of target exons loaded from the file provided by the user
    :return: A matrix where rows represent the region index and the first column as a list of targets in this region
    '''
    regions_count -= 1
    number_of_targets = len(targets_list)
    combine_size = number_of_targets / regions_count
    cnv_regions = []

    i = 0
    while True:
        if i == len(targets_list):
            break

        region_start = targets_list[i][1]
        region_end = region_start + random.randint(region_min_length, region_max_length) - 1
        cnv_regions.append((targets_list[i][0], region_start, region_end))

        i += combine_size
        if i > len(targets_list):
            i = len(targets_list)

    cnv_matrix = _loadCNVMatrix(targets_list, cnv_regions)

    return cnv_matrix, cnv_regions


def _loadCNVMatrix(targets_list, cnv_regions):
    '''
    A function to map targets to their corresponding CNV regions
    :param targets_list: a list of target exons loaded from the file provided by the user
    :param cnv_regions: a list of CNV regions loaded from the CNV file provided by the user
    :return: A matrix where rows represent the region index and the first column as a list of targets in this region
    '''
    cnv_matrix = []
    for cnv in cnv_regions:
        region_chromosome = cnv[0]
        region_start = cnv[1]
        region_end = cnv[2]

        region_targets = []
        for target in targets_list:
            target_chromosome = target[0]
            target_start = target[1]
            target_end = target[2]
            if target_chromosome == region_chromosome and target_start >= region_start:
                if target_end <= region_end:
                    region_targets.append(target)
                elif target_start <= region_end:
                    truncated_target = (target_chromosome, target_start, region_end)
                    region_targets.append(truncated_target)

        cnv_matrix.append(region_targets)

    return cnv_matrix


def _generateCNVMask(number_of_exons, p_amplify, p_delete, min_variation, max_variation):
    '''
    This function generates random Copy Number Variations mask list
    :param number_of_exons: the length of the target regions.
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


def _simulateCNV(genome, cnv_matrix, mask, read_length):
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

    CONTROL_APPEND_INDEX = 0      # a value to re-base the (start, end) of amplified/deleted targets in a region
    CNV_APPEND_INDEX = 0          # a value to re-base the (start, end) of amplified/deleted targets in a region

    for k, cnv_region in enumerate(cnv_matrix):

        number_of_copies = mask[k]

        if number_of_copies > 0:
            # if amplification
            for target in cnv_region:
                for _ in range(number_of_copies):
                    amplification = genome[target[1]-read_length:target[2]+read_length]
                    cnv_genome.append(amplification)
                    chromosome = target[0]
                    start = CNV_APPEND_INDEX + read_length
                    end = CNV_APPEND_INDEX + len(amplification) - read_length
                    cnv_targets.append((chromosome, start, end))
                    CNV_APPEND_INDEX += len(amplification)
        elif number_of_copies < 0:
            # if deletion
            for target in cnv_region:
                for _ in range(abs(number_of_copies)):
                    deletion = genome[target[1]-read_length:target[2]+read_length]
                    control_genome.append(deletion)
                    chromosome = target[0]
                    start = CONTROL_APPEND_INDEX + read_length
                    end = CONTROL_APPEND_INDEX + len(deletion) - read_length
                    control_targets.append((chromosome, start, end))
                    CONTROL_APPEND_INDEX += len(deletion)

    return ''.join(control_genome), control_targets, ''.join(cnv_genome), cnv_targets


def _callWessim(genome_file, target_file, output_file, number_of_reads, read_length, model_file="models/ill100v4_p.gzip", number_of_threads=4):
    '''
    Calls Wessim to generate artificial reads for the targets on the reference genome
    :param genome_file: reference genome file in FASTA format
    :param target_file: target regions file in BED format
    :param output_file: output file name
    :param number_of_reads: the number of reads to generate for this simulation
    :param read_length: the read length
    :param model_file: GemSim's empirical error models used by Wessim for NGS read generation
    :param number_of_threads: the number of threads to run this simulation
    :return: None
    '''
    os.chdir(os.path.join(getScriptPath(), "Wessim" ))
    subprocess.call(["python", "Wessim1.py", \
                     "-R" + genome_file, \
                     "-B" + target_file, \
                     "-n" + str(number_of_reads), \
                     "-l" + str(read_length), \
                     "-M" + model_file, \
                     "-o" + output_file, \
                     "-t" + str(number_of_threads), \
                     "-p"], stderr=None)
    os.chdir('..')


def simulate_exome_cnv(simulation_parameters, cnv_list_parameters=None):
    '''
    Simulate copy number variations on the passed reference genome based on the given simulation parameters
    :param genome: path to the FASTA file containing the reference genome
    :param simulation_parameters: a dictionary containing parameters for simulation
    :return: simulation status {'simulation_status': bool, 'message': str}
    '''
    _log('simulation type: targeted exome sequencing')

    # create a temporary directory for intermediate files
    if not os.path.exists(simulation_parameters['tmp_dir']):
        os.makedirs(simulation_parameters['tmp_dir'])

    # create output directory for final results
    if not os.path.exists(simulation_parameters['output_dir']):
        os.makedirs(simulation_parameters['output_dir'])

    if simulation_parameters['target_file'] is None:
        _log("ERROR: target file cannot be None for targeted exome simulation!")
        fileio.clean(simulation_parameters['output_dir'])
        fileio.clean(simulation_parameters['tmp_dir'])
        exit()
    target_file = simulation_parameters['target_file']

    # copy genome and target to the tmp folder
    shutil.copyfile(simulation_parameters['genome_file'], os.path.join(simulation_parameters['tmp_dir'], "reference.fa"))
    shutil.copyfile(target_file, os.path.join(simulation_parameters['tmp_dir'], "target.bed"))
    genome_file = os.path.join(simulation_parameters['tmp_dir'], "reference.fa")
    target_file = os.path.join(simulation_parameters['tmp_dir'], "target.bed")

    # initialize variables for temporary files
    control_genome_file = os.path.join(simulation_parameters['tmp_dir'], "ControlGenome.fa")
    control_target_file = os.path.join(simulation_parameters['tmp_dir'], "ControlTarget.bed")
    cnv_genome_file = os.path.join(simulation_parameters['tmp_dir'], "CNVGenome.fa")
    cnv_target_file = os.path.join(simulation_parameters['tmp_dir'], "CNVTarget.bed")
    base_reads_file = os.path.join(simulation_parameters['tmp_dir'], "base")
    control_reads_file = os.path.join(simulation_parameters['tmp_dir'], "control")
    cnv_reads_file = os.path.join(simulation_parameters['tmp_dir'], "cnv")

    _log("loading genome file ..")
    header, genome = fileio.readGenome(genome_file)
    _log("successfully loaded a genome of length " + `len(genome)`)

    _log("loading target file ..")
    _log("sorting and merging targets ..")
    target_file = fileio.prepareTargetFile(target_file)
    targets = fileio.readTargets(target_file)
    _log("successfully loaded " + `len(targets)` + " targets ..")

    if simulation_parameters['cnv_list_file'] is None:
        # CNV list file
        cnv_list_file = os.path.join(simulation_parameters['output_dir'], "copynumber.bed")

        _log("generating CNV list ..")
        cnv_matrix, cnv_regions = _generateCNVMatrix(targets, cnv_list_parameters['regions_count'], \
                                        cnv_list_parameters['minimum_length'], cnv_list_parameters['maximum_length'])
        mask = _generateCNVMask(len(cnv_matrix), cnv_list_parameters['amplifications'], \
                                cnv_list_parameters['deletions'], cnv_list_parameters['minimum_variations'], \
                                               cnv_list_parameters['maximum_variations'])

        with open(cnv_list_file, 'w') as f:
            line = 'chrom\tchr_start\tchrom_stop\tnum_positions\tcopy_number\n'
            f.write(line)
            for i, cnv_region in enumerate(cnv_regions):
                num_positions = cnv_region[2] - cnv_region[1] + 1
                line = cnv_region[0] + '\t' \
                       + `cnv_region[1]` + '\t' \
                       + `cnv_region[2]` + '\t' \
                       + `num_positions` + '\t' \
                       + `mask[i]` + '\n'
                f.write(line)
        _log("randomly generated CNV list saved to " + cnv_list_file)
    else:
        _log("loading CNV list ..")
        with open(simulation_parameters['cnv_list_file'], "r") as f:
            cnv_list = []
            lines = f.readlines()
            lines.pop(0)
            for line in lines:
                chromosome, region_start, region_end, num_positions, variation = line.strip().split("\t")
                cnv_list.append((chromosome, int(region_start), int(region_end), int(variation)))

            cnv_matrix = _loadCNVMatrix(targets, cnv_list)
            mask = map(lambda x: x[3], cnv_list)
        _log("successfully loaded CNV list that contains " + `len(lines)` + " regions ..")

    # call Wessim to generate reads from the genome file and the target list
    _log("generating reads for the target exons ..")
    _log("delegating job to Wessim ...")
    _callWessim(genome_file, target_file, base_reads_file, simulation_parameters['number_of_reads'], simulation_parameters['read_length'])

    # generate genome and target files for the control and the CNV
    _log("simulating copy number variations (amplifications/deletions)")
    control_genome, control_targets, cnv_genome, cnv_targets = _simulateCNV(genome, cnv_matrix, mask, simulation_parameters['read_length'])

    _log("saving to the control genome file ..")
    with open(control_genome_file, 'w') as fw:
        fw.write(header + "\n")
        n = 50
        l = [control_genome[i:i + n] for i in range(0, len(control_genome), n)]
        for line in l:
            fw.write(line + "\n")

    _log("saving to the control target file ..")
    with open(control_target_file, 'w') as tw:
        for target in control_targets:
            line = target[0] + "\t" + `target[1]` + "\t" + `target[2]` + "\n"
            tw.write(line)

    _log("delegating job to Wessim ...")
    _callWessim(control_genome_file, control_target_file, control_reads_file,
                               int(cnv_list_parameters['deletions'] * simulation_parameters['number_of_reads']), simulation_parameters['read_length'])

    _log("saving to the CNV genome file ..")
    with open(cnv_genome_file, 'w') as fw:
        fw.write(header + "\n")
        n = 50
        l = [cnv_genome[i:i + n] for i in range(0, len(cnv_genome), n)]
        for line in l:
            fw.write(line + "\n")

    _log("saving to the CNV target file ..")
    with open(cnv_target_file, 'w') as tw:
        for target in cnv_targets:
            line = target[0] + "\t" + `target[1]` + "\t" + `target[2]` + "\n"
            tw.write(line)

    _log("delegating job to Wessim ...")
    _callWessim(cnv_genome_file, cnv_target_file, cnv_reads_file,
                               int(cnv_list_parameters['amplifications']* simulation_parameters['number_of_reads']), simulation_parameters['read_length'])

    _log("merging results ..")
    fileio.mergeWessimReads(simulation_parameters['tmp_dir'], simulation_parameters['output_dir'])

    _log("cleaning temporary files ..")
    fileio.clean(simulation_parameters['tmp_dir'])

    _log("simulation completed. find results in " + simulation_parameters['output_dir'])