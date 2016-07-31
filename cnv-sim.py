#!/usr/bin/python

__author__ = 'Abdelrahman Hosny'

import os.path
import datetime
import argparse
import shutil

from lib import fileio
from lib import exome_simulator
from lib import genome_simulator

def log(message):
    print '[CNV SIM {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + "] " + message

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("simulation_type", type=str, choices=['genome', 'exome'], \
                        help="simulate copy number variations in whole genome or exome regions")
    parser.add_argument("genome", type=file, \
                        help="path to the referece genome file in FASTA format ")
    parser.add_argument("target", type=file, nargs='?', default=None, \
                        help="path to the target regions file in BED format (if using exome)")

    parser.add_argument("-o", "--output_dir",type=str, default="test", \
                        help="a name to be used to create the output directory (overrides existing directory with the same name).")
    parser.add_argument("--cnv_list", type=file, default=None, \
                       help="path to a CNV list file in BED format chr | start | end | variation. If not passed, it is randomly generated using CNV list parameters below")
    parser.add_argument("-n", "--n_reads", type=int, default=10000, \
                        help="total number of reads without variations")

    cnv_sim_group = parser.add_argument_group('CNV list parameters', "parameters to be used if CNV list is not passed")
    cnv_sim_group.add_argument("-g", "--regions_count", type=int, default=30, \
                        help="number of CNV regions to be randomly generated")
    cnv_sim_group.add_argument("-a", "--amplifications", type=float, default=0.50, \
                        help="percentage of amplifications in range [0.0: 1.0].")
    cnv_sim_group.add_argument("-d", "--deletions", type=float, default=0.20, \
                        help="percentage of deletions in range [0.0: 1.0].")
    cnv_sim_group.add_argument("-min", "--minimum", type=float, default=3, \
                        help="minimum number of amplifications/deletions introduced")
    cnv_sim_group.add_argument("-max", "--maximum", type=float, default=10, \
                        help="maximum number of amplifications/deletions introduced")

    args = parser.parse_args()

    genome_file = args.genome.name
    simulation_name = args.output_dir

    cnv_list_file = args.cnv_list
    number_of_reads = args.n_reads
    regions_count = args.regions_count

    output_dir = os.path.join(os.getcwd(), simulation_name)
    tmp_dir = os.path.join(os.getcwd(), "tmp")

    # create a temporary directory for intermediate files
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # create output directory for final results
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if args.simulation_type == 'genome':
        log('simulation type: whole genome')

        # copy genome to the tmp folder
        shutil.copyfile(genome_file, os.path.join(tmp_dir, "reference.fa"))
        genome_file = os.path.join(tmp_dir, "reference.fa")

        # initialize variables for temporary files
        control_genome_file = os.path.join(tmp_dir, "ControlGenome.fa")
        cnv_genome_file = os.path.join(tmp_dir, "CNVGenome.fa")
        base_reads_file = os.path.join(tmp_dir, "base")
        control_reads_file = os.path.join(tmp_dir, "control")
        cnv_reads_file = os.path.join(tmp_dir, "cnv")


        log("loading genome file ..")
        header, genome = fileio.readGenome(genome_file)
        log("successfully loaded a genome of length " + `len(genome)`)

        if cnv_list_file is None:
            # CNV list file
            cnv_list_file = os.path.join(output_dir, "CNVList.bed")

            log("generating CNV list ..")
            cnv_list = genome_simulator.generateCNVList(len(genome), regions_count, args.amplifications, args.deletions, args.minimum, args.maximum)
            cnv_list = map(lambda l:[header.replace('>', '')] + l , cnv_list)

            with open(cnv_list_file, 'w') as f:
                for i, cnv_region in enumerate(cnv_list):
                    line = cnv_region[0] + '\t' \
                           + `cnv_region[1]` + '\t' \
                           + `cnv_region[2]` + '\t' \
                           + `cnv_region[3]` + '\n'
                    f.write(line)
            log("randomly generated CNV list saved to " + cnv_list_file)

        else:
            log("loading CNV list ..")
            with open(cnv_list_file.name, "r") as f:
                cnv_list = []
                lines = f.readlines()
                for line in lines:
                    chromosome, region_start, region_end, variation = line.strip().split("\t")
                    cnv_list.append((chromosome, int(region_start), int(region_end), int(variation)))

            log("successfully loaded CNV list that contains " + `len(lines)` + " regions ..")

        # call ART to generate reads from the genome file
        log("generating reads for the genome ..")
        log("delegating job to ART ...")
        genome_simulator.callART(genome_file, base_reads_file)

        log("simulating copy number variations (amplifications/deletions)")
        control_genome, cnv_genome = genome_simulator.simulateCNV(genome, cnv_list)

        log("saving to the control genome file ..")

        with open(control_genome_file, 'w') as fw:
            fw.write(header + "\n")
            n = 50
            l = [control_genome[i:i + n] for i in range(0, len(control_genome), n)]
            for line in l:
                fw.write(line + "\n")

        log("delegating job to ART ...")
        genome_simulator.callART(control_genome_file, control_reads_file)

        log("saving to the CNV genome file ..")
        with open(cnv_genome_file, 'w') as fw:
            fw.write(header + "\n")
            n = 50
            l = [cnv_genome[i:i + n] for i in range(0, len(cnv_genome), n)]
            for line in l:
                fw.write(line + "\n")


        log("delegating job to ART ...")
        genome_simulator.callART(cnv_genome_file, cnv_reads_file)

        log("merging results ..")
        fileio.mergeARTReads(tmp_dir, output_dir)

        log("cleaning temporary files ..")
        fileio.clean(tmp_dir)

        log("simulation completed. find results in " + output_dir)

    else:
        log('simulation type: whole exome')
        if args.target == None:
            log("ERROR: target file cannot be None for whole exome simulation!")
            fileio.clean(output_dir)
            fileio.clean(tmp_dir)
            exit()
        target_file = args.target.name

        # copy genome and target to the tmp folder
        shutil.copyfile(genome_file, os.path.join(tmp_dir, "reference.fa"))
        shutil.copyfile(target_file, os.path.join(tmp_dir, "target.bed"))
        genome_file = os.path.join(tmp_dir, "reference.fa")
        target_file = os.path.join(tmp_dir, "target.bed")

        # initialize variables for temporary files
        control_genome_file = os.path.join(tmp_dir, "ControlGenome.fa")
        control_target_file = os.path.join(tmp_dir, "ControlTarget.bed")
        cnv_genome_file = os.path.join(tmp_dir, "CNVGenome.fa")
        cnv_target_file = os.path.join(tmp_dir, "CNVTarget.bed")
        base_reads_file = os.path.join(tmp_dir, "base")
        control_reads_file = os.path.join(tmp_dir, "control")
        cnv_reads_file = os.path.join(tmp_dir, "cnv")

        log("loading genome file ..")
        header, genome = fileio.readGenome(genome_file)
        log("successfully loaded a genome of length " + `len(genome)`)

        log("loading target file ..")
        log("sorting and merging targets ..")
        target_file = fileio.prepareTargetFile(target_file)
        targets = fileio.readTargets(target_file)
        log("successfully loaded " + `len(targets)` + " targets ..")

        if cnv_list_file is None:
            # CNV list file
            cnv_list_file = os.path.join(output_dir, "CNVList.bed")

            log("generating CNV list ..")
            cnv_matrix = exome_simulator.generateCNVMatrix(targets, regions_count)
            mask = exome_simulator.generateCNVMask(len(cnv_matrix), args.amplifications, args.deletions, args.minimum, args.maximum)

            with open(cnv_list_file, 'w') as f:
                for i, cnv_region in enumerate(cnv_matrix):
                    line = cnv_region[0][0] + '\t' \
                           + `cnv_region[0][1]` + '\t' \
                           + `cnv_region[-1][2]` + '\t' \
                           + `mask[i]` + '\n'
                    f.write(line)
            log("randomly generated CNV list saved to " + cnv_list_file)
        else:
            log("loading CNV list ..")
            with open(cnv_list_file.name, "r") as f:
                cnv_list = []
                lines = f.readlines()
                for line in lines:
                    chromosome, region_start, region_end, variation = line.strip().split("\t")
                    cnv_list.append((chromosome, int(region_start), int(region_end), int(variation)))

                cnv_matrix = exome_simulator.loadCNVMatrix(targets, cnv_list)
                mask = map(lambda x: x[3], cnv_list)
            log("successfully loaded CNV list that contains " + `len(lines)` + " regions ..")

        # call Wessim to generate reads from the genome file and the target list
        log("generating reads for the target exons ..")
        log("delegating job to Wessim ...")
        exome_simulator.callWessim(genome_file, target_file, base_reads_file, number_of_reads)

        # generate genome and target files for the control and the CNV
        log("simulating copy number variations (amplifications/deletions)")
        control_genome, control_targets, cnv_genome, cnv_targets = exome_simulator.simulateCNV(genome, cnv_matrix, mask)

        log("saving to the control genome file ..")
        with open(control_genome_file, 'w') as fw:
            fw.write(header + "\n")
            n = 50
            l = [control_genome[i:i + n] for i in range(0, len(control_genome), n)]
            for line in l:
                fw.write(line + "\n")

        log("saving to the control target file ..")
        with open(control_target_file, 'w') as tw:
            for target in control_targets:
                line = target[0] + "\t" + `target[1]` + "\t" + `target[2]` + "\n"
                tw.write(line)

        log("delegating job to Wessim ...")
        exome_simulator.callWessim(control_genome_file, control_target_file, control_reads_file, int(args.deletions * number_of_reads))

        log("saving to the CNV genome file ..")
        with open(cnv_genome_file, 'w') as fw:
            fw.write(header + "\n")
            n = 50
            l = [cnv_genome[i:i + n] for i in range(0, len(cnv_genome), n)]
            for line in l:
                fw.write(line + "\n")

        log("saving to the CNV target file ..")
        with open(cnv_target_file, 'w') as tw:
            for target in cnv_targets:
                line = target[0] + "\t" + `target[1]` + "\t" + `target[2]` + "\n"
                tw.write(line)

        log("delegating job to Wessim ...")
        exome_simulator.callWessim(cnv_genome_file, cnv_target_file, cnv_reads_file, int(args.amplifications * number_of_reads))

        log("merging results ..")
        fileio.mergeWessimReads(tmp_dir,output_dir)

        log("cleaning temporary files ..")
        fileio.clean(tmp_dir)

        log("simulation completed. find results in " + output_dir)

if __name__ == '__main__':
    main()