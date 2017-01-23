#!/usr/bin/python

import random
import datetime
import os
import shutil

from . import fileio

__author__ = 'Abdelrahman Hosny'


def _log(message):
    print '[CNV SIM {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + "] " + message


def _changeNucleotide(nucleotide):
    '''
    Replaces the passed nucleotide with a randomly chosen other letter
    :param nucleotide: character A, C, G or T
    :return: a letter other than the one passed
    '''
    letters = ['A', 'C', 'G', 'T']
    if nucleotide in letters:
        letters.remove(nucleotide)
        random.shuffle(letters)
        return letters[0]
    else:
        return nucleotide


def _simulateSNPs(genome, rate):
    '''
    Simulates the SNPs in a given genome reference
    :param genome: a string of the genome sequence
    :param rate: the number of SNPs to be introduced to the genome
    :return: a modified genome string harboring SNPs
    '''
    number_of_snps = int(rate * len(genome))
    positions = [random.randint(0, len(genome)) for _ in range(number_of_snps)]
    modified_genome = list(genome)
    for i in positions:
        original = modified_genome[i]
        snp = _changeNucleotide(original)
        modified_genome[i] = snp
    return "".join(modified_genome)


def simulate_snp(simulation_parameters):
    '''
    Introduce SNPs to the original genome reference as per the rate specified by the user
    :param simulation_parameters: the dictionary of all simulation parameters
    :return: modified dictionary of the simulation parameters
    '''
    _log('simulating SNPs ..')

    # create a temporary directory for intermediate files
    if not os.path.exists(simulation_parameters['tmp_dir']):
        os.makedirs(simulation_parameters['tmp_dir'])

    header, genome = fileio.readGenome(simulation_parameters['genome_file'])

    snp_genome = _simulateSNPs(genome, simulation_parameters['snp'])

    snp_genome_file = os.path.join(simulation_parameters['tmp_dir'], "SNP.fa")

    _log("finished simulating SNPs")
    _log("saving to the SNP genome file")

    with open(snp_genome_file, 'w') as fw:
        fw.write(header + "\n")
        n = 50
        l = [snp_genome[i:i + n] for i in range(0, len(snp_genome), n)]
        for line in l:
            fw.write(line + "\n")

    simulation_parameters['genome_file'] = snp_genome_file
    return simulation_parameters


def _simulateIndels(genome, targets, rate):
    '''
    Introduce indels to the original genome reference as per the rate specified by the user
    :param genome: a string of the genome sequence
    :param targets: list of targets to modify their positions according to the insertions and deletions
    :param rate: the number of indels to be introduced to the genome
    :return: a modified genome string harboring indels, and a modified target list.
    '''
    number_of_indels = int(rate * len(genome))
    positions = [random.randint(0, len(genome)) for _ in range(number_of_indels)]

    # complete here ya Abdo

    return genome, targets

def simulate_indels(simulation_parameters):
    '''
    Introduce indels to the original genome reference as per the rate specified by the user
    :param simulation_parameters: the dictionary of all simulation parameters
    :return: modified dictionary of the simulation parameters
    '''
    _log('simulating SNPs ..')

    # create a temporary directory for intermediate files
    if not os.path.exists(simulation_parameters['tmp_dir']):
        os.makedirs(simulation_parameters['tmp_dir'])

    header, genome = fileio.readGenome(simulation_parameters['genome_file'])

    target_file = fileio.prepareTargetFile(simulation_parameters['target_file'])
    targets = fileio.readTargets(target_file)

    indel_genome, indel_targets = _simulateIndels(genome, targets, simulation_parameters['indel'])

    indel_genome_file = os.path.join(simulation_parameters['tmp_dir'], "indel.fa")
    indel_target_file = os.path.join(simulation_parameters['tmp_dir'], "indel_targets.bed")

    _log("finished simulating indels")
    _log("saving to the indel genome file")

    with open(indel_genome_file, 'w') as fw:
        fw.write(header + "\n")
        n = 50
        l = [indel_genome[i:i + n] for i in range(0, len(indel_genome), n)]
        for line in l:
            fw.write(line + "\n")

    _log("saving to the indel target file ..")
    with open(indel_target_file, 'w') as tw:
        for target in indel_targets:
            line = target[0] + "\t" + `target[1]` + "\t" + `target[2]` + "\n"
            tw.write(line)

    simulation_parameters['genome_file'] = indel_genome_file
    simulation_parameters['target_file'] = indel_target_file

    return simulation_parameters