#!/usr/bin/python

__author__ = 'Abdelrahman Hosny'

import os.path
import datetime
import argparse
import shutil

from cnvsim.fileio import *
from cnvsim.exome_simulator import *
from cnvsim.genome_simulator import *
from cnvsim.sv_simulator import *

class CapitalisedHelpFormatter(argparse.HelpFormatter):
    def add_usage(self, usage, actions, groups, prefix=None):
        if prefix is None:
            prefix = 'Usage: '
            return super(CapitalisedHelpFormatter, self).add_usage(usage, actions, groups, prefix)

def log(message):
    print '[CNV SIM {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + "] " + message


def main():
    parser = argparse.ArgumentParser(add_help=True, formatter_class=CapitalisedHelpFormatter, \
                                     description='Generates NGS short reads that encompass copy number variations in whole genome and targeted exome sequencing')
    parser._positionals.title = 'Positional arguments'
    parser._optionals.title = 'Optional arguments'
    parser.add_argument('-v', '--version', action='version', version = 'CNV-Sim v1.0.0', help = "Show program's version number and exit.")

    parser.add_argument("simulation_type", type=str, choices=['exome', 'genome'], \
                        help="simulate copy number variations in targeted exomes or whole genome")
    parser.add_argument("genome", type=file, \
                        help="path to the referece genome file in FASTA format ")
    parser.add_argument("target", type=file, nargs='?', default=None, \
                        help="path to the target regions file in BED format (if using exome)")

    parser.add_argument("-o", "--output_dir_name",type=str, default="simulation_output", \
                        help="a name to be used to create the output directory (overrides existing directory with the same name)")
    parser.add_argument("-n", "--n_reads", type=int, default=10000, \
                        help="total number of reads without variations")
    parser.add_argument("-l", "--read_length", type=int, default=100, \
                        help="read length (bp)")
    parser.add_argument("--cnv_list", type=file, default=None, \
                        help="path to a CNV list file in BED format chr | start | end | variation. If not passed, it is randomly generated using CNV list parameters below")

    cnv_sim_group = parser.add_argument_group('CNV list parameters', "parameters to be used if CNV list is not passed")
    cnv_sim_group.add_argument("-g", "--regions_count", type=int, default=20, \
                        help="number of CNV regions to be generated randomly")
    cnv_sim_group.add_argument("-r_min", "--region_minimum_length", type=int, default=1000, \
                               help="minimum length of each CNV region")
    cnv_sim_group.add_argument("-r_max", "--region_maximum_length", type=int, default=100000, \
                               help="maximum length of each CNV region")
    cnv_sim_group.add_argument("-a", "--amplifications", type=float, default=0.30, \
                        help="percentage of amplifications in range [0.0: 1.0].")
    cnv_sim_group.add_argument("-d", "--deletions", type=float, default=0.20, \
                        help="percentage of deletions in range [0.0: 1.0].")
    cnv_sim_group.add_argument("-cn_min", "--copy_number_minimum", type=float, default=3, \
                        help="minimum level of variations (copy number) introduced")
    cnv_sim_group.add_argument("-cn_max", "--copy_number_maximum", type=float, default=10, \
                        help="maximum level of variation (copy number) introduced")

    cnv_sim_group = parser.add_argument_group('Tumor parameters', "parameters to describe tumor complexity")
    cnv_sim_group.add_argument("-snp", "--snp_rate", type=float, default=0.1, \
                        help="rate of structural variations represented as SNPs")
    cnv_sim_group.add_argument("-indel", "--indel_rate", type=float, default=0.1, \
                        help="rate of structural variations represented as insertions/deletions")

    args = parser.parse_args()

    simulation_parameters = {}
    simulation_parameters['type'] = args.simulation_type
    simulation_parameters['genome_file'] = args.genome.name
    if args.target is not None:
        simulation_parameters['target_file'] = args.target.name
    else:
        simulation_parameters['target_file'] = None
    simulation_parameters['output_dir'] = os.path.join(os.getcwd(), args.output_dir_name)
    simulation_parameters['number_of_reads'] = args.n_reads
    simulation_parameters['read_length'] = args.read_length
    if args.cnv_list is not None:
        simulation_parameters['cnv_list_file'] = args.cnv_list.name
    else:
        simulation_parameters['cnv_list_file'] = None
    simulation_parameters['tmp_dir'] = os.path.join(os.getcwd(), args.output_dir_name , "tmp")

    cnv_list_parameters = {}
    cnv_list_parameters['regions_count'] = args.regions_count
    cnv_list_parameters['minimum_length'] = args.region_minimum_length
    cnv_list_parameters['maximum_length'] = args.region_maximum_length
    cnv_list_parameters['amplifications'] = args.amplifications
    cnv_list_parameters['deletions'] = args.deletions
    cnv_list_parameters['minimum_variations'] = args.copy_number_minimum
    cnv_list_parameters['maximum_variations'] = args.copy_number_maximum

    if cnv_list_parameters['amplifications'] + cnv_list_parameters['deletions'] > 1.0:
        log("ERROR: percentage of amplifications + percentage of deletions must be less than or equal to 1.0")
        exit()
    cnv_list_parameters['snp'] = args.snp
    cnv_list_parameters['indel'] = args.indel
    
    simulation_parameters = simulate_snp(simulation_parameters)

    if simulation_parameters['type'] == 'genome':
        simulate_genome_cnv(simulation_parameters, cnv_list_parameters)
    else:
        simulate_exome_cnv(simulation_parameters, cnv_list_parameters)


if __name__ == '__main__':
    main()
