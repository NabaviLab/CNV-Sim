__author__ = 'Abdelrahman Hosny'

import random
import os
import subprocess
import datetime
import shutil

from . import fileio


def _log(message):
    print '[CNV SIM {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + "] " + message


def getScriptPath():
    return os.path.dirname(os.path.realpath(__file__))


def _generateCNVMask(mask_length, p_amplify, p_delete, min_variation, max_variation):
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


def _generateCNVList(chromosome_length, number_of_regions, p_amplify, p_delete, min_variation, max_variation):
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
    mask = _generateCNVMask(number_of_regions, p_amplify, p_delete, min_variation, max_variation)
    for i in range(number_of_regions):
        # jump forward start
        jump_start = start + int(random.randrange(int(0.40*region_length), int(0.45*region_length)))

        # jump backward end
        jump_end = end - int(random.randrange(int(0.40*region_length), int(0.45*region_length)))

        cnv_list.append([jump_start, jump_end, mask[i]])
        start = end
        end += region_length

    return cnv_list


def _simulateCNV(genome, cnv_list, read_length):
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
            amplification = genome[start-read_length: start] + sequence * variation + genome[end:end+read_length]
            cnv_genome.append(amplification)
        elif variation < 0:
            # deletion
            deletion = genome[start-read_length: start] + sequence * variation + genome[end:end+read_length]
            control_genome.append(deletion)

    return ''.join(control_genome), ''.join(cnv_genome)


def _callART(genome_file, output_file, read_length, fold_coverage=1):
    '''
    Calls Wessim to generate artificial reads for the targets on the reference genome
    :param genome_file: reference genome file in FASTA format
    :param output_file: output file name
    :param read_length: the read length
    :param fold_coverage: fold coverage for the reads
    :return: None
    '''
    os.chdir(os.path.join(getScriptPath(), "ART"))
    subprocess.call(["./art_illumina", \
                     "-sam", \
                     "-i", genome_file, \
                     "-l", str(read_length), \
                     "-f", str(fold_coverage), \
                     "-o", output_file], stderr=None)
    os.chdir("..")


def simulate_genome_cnv(simulation_parameters, cnv_list_parameters=None):
    '''
    Simulate copy number variations on the passed reference genome based on the given simulation parameters
    :param simulation_parameters: a dictionary containing parameters for simulation
    :param cnv_list_parameters: a dictionary containing parameters for CNV List creation
    :return: None
    '''
    _log('simulation type: whole genome')

    # create a temporary directory for intermediate files
    if not os.path.exists(simulation_parameters['tmp_dir']):
        os.makedirs(simulation_parameters['tmp_dir'])

    # create output directory for final results
    if not os.path.exists(simulation_parameters['output_dir']):
        os.makedirs(simulation_parameters['output_dir'])

    # copy genome to the tmp folder
    shutil.copyfile(simulation_parameters['genome_file'], os.path.join(simulation_parameters['tmp_dir'], "reference.fa"))
    genome_file = os.path.join(simulation_parameters['tmp_dir'], "reference.fa")

    # initialize variables for temporary files
    control_genome_file = os.path.join(simulation_parameters['tmp_dir'], "ControlGenome.fa")
    cnv_genome_file = os.path.join(simulation_parameters['tmp_dir'], "CNVGenome.fa")
    base_reads_file = os.path.join(simulation_parameters['tmp_dir'], "base")
    control_reads_file = os.path.join(simulation_parameters['tmp_dir'], "control")
    cnv_reads_file = os.path.join(simulation_parameters['tmp_dir'], "cnv")

    _log("loading genome file ..")
    header, genome = fileio.readGenome(genome_file)
    _log("successfully loaded a genome of length " + `len(genome)`)

    if simulation_parameters['cnv_list_file'] is None:
        # CNV list file
        cnv_list_file = os.path.join(simulation_parameters['output_dir'], "CNVList.bed")

        _log("generating CNV list ..")
        cnv_list = _generateCNVList(len(genome), cnv_list_parameters['regions_count'], \
                                    cnv_list_parameters['amplifications'], cnv_list_parameters['deletions'], \
                                    cnv_list_parameters['minimum_variations'], \
                                    cnv_list_parameters['maximum_variations'])
        cnv_list = map(lambda l: [header.replace('>', '')] + l, cnv_list)

        with open(cnv_list_file, 'w') as f:
            for i, cnv_region in enumerate(cnv_list):
                line = cnv_region[0] + '\t' \
                       + `cnv_region[1]` + '\t' \
                       + `cnv_region[2]` + '\t' \
                       + `cnv_region[3]` + '\n'
                f.write(line)
        _log("randomly generated CNV list saved to " + cnv_list_file)

    else:
        _log("loading CNV list ..")
        with open(simulation_parameters['cnv_list_file'], "r") as f:
            cnv_list = []
            lines = f.readlines()
            for line in lines:
                chromosome, region_start, region_end, variation = line.strip().split("\t")
                cnv_list.append((chromosome, int(region_start), int(region_end), int(variation)))

        _log("successfully loaded CNV list that contains " + `len(lines)` + " regions ..")

    # call ART to generate reads from the genome file
    _log("generating reads for the genome ..")
    _log("delegating job to ART ...")
    _callART(genome_file, base_reads_file, simulation_parameters['read_length'])

    _log("simulating copy number variations (amplifications/deletions)")
    control_genome, cnv_genome = _simulateCNV(genome, cnv_list, simulation_parameters['read_length'])

    _log("saving to the control genome file ..")

    with open(control_genome_file, 'w') as fw:
        fw.write(header + "\n")
        n = 50
        l = [control_genome[i:i + n] for i in range(0, len(control_genome), n)]
        for line in l:
            fw.write(line + "\n")

    _log("delegating job to ART ...")
    _callART(control_genome_file, control_reads_file, simulation_parameters['read_length'])

    _log("saving to the CNV genome file ..")
    with open(cnv_genome_file, 'w') as fw:
        fw.write(header + "\n")
        n = 50
        l = [cnv_genome[i:i + n] for i in range(0, len(cnv_genome), n)]
        for line in l:
            fw.write(line + "\n")

    _log("delegating job to ART ...")
    _callART(cnv_genome_file, cnv_reads_file, simulation_parameters['read_length'])

    _log("merging results ..")
    fileio.mergeARTReads(simulation_parameters['tmp_dir'], simulation_parameters['output_dir'])

    _log("cleaning temporary files ..")
    fileio.clean(simulation_parameters['tmp_dir'])

    _log("simulation completed. find results in " + simulation_parameters['output_dir'])