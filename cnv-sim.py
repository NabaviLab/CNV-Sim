#!/usr/bin/python

__author__ = 'Abdelrahman Hosny'

import os.path
import datetime
import argparse
import shutil

from lib import fileio
from lib import exome_simulator

def log(message):
    print '[CNV SIM {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + "] " + message

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("genome", type=file, \
                        help="path to the referece genome file in FASTA format ")
    parser.add_argument("target", type=file, \
                        help="path to the target regions file in BED format")

    parser.add_argument("-m", "--name",type=str, default="test", \
                        help="a name to be used for simulated results.")
    # parser.add_argument("--cnv-list", type=str, \
    #                   help="path to a CNV list file in BED format chr | start | end | variation. If not passed, it is randomly generated using --amplifications and --deletions parameters")
    parser.add_argument("-n", "--n_reads", type=int, default=10000, \
                        help="total number of reads without variations")
    parser.add_argument("-a", "--amplifications", type=float, default=0.50, \
                        help="percentage of amplifications in range [0.0: 1.0].")
    parser.add_argument("-d", "--deletions", type=float, default=0.20, \
                        help="percentage of deletions in range [0.0: 1.0].")
    parser.add_argument("-min", "--minimum", type=float, default=10, \
                        help="minimum number of amplifications/deletions introduced")
    parser.add_argument("-max", "--maximum", type=float, default=15, \
                        help="maximum number of amplifications/deletions introduced")

    args = parser.parse_args()

    genome_file = args.genome.name
    target_file = args.target.name
    simulation_name = args.name
    number_of_reads = args.n_reads

    output_dir = os.path.join(os.getcwd(), simulation_name)
    tmp_dir = os.path.join(os.getcwd(), "tmp")

    # create a temporary directory for intermediate files
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # create output directory for final results
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # copy genome and target to the tmp folder
    shutil.copyfile(genome_file, os.path.join(tmp_dir, "reference.fa"))
    shutil.copyfile(target_file, os.path.join(tmp_dir, "target.bed"))
    genome_file = os.path.join(tmp_dir, "reference.fa")
    target_file = os.path.join(tmp_dir, "target.bed")

    # CNV list file
    cnv_list_file = os.path.join(output_dir, "CNVList.bed")

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

    log("sorting and merging targets ..")
    target_file = fileio.prepareTargetFile(target_file)

    log("loading target file ..")
    targets = fileio.readTargets(target_file)

    log("generating CNV list ..")

    cnv_matrix = exome_simulator.getCNVMatrix(targets)
    mask = exome_simulator.generateCNVMask(len(cnv_matrix), args.amplifications, args.deletions, args.minimum, args.maximum)

    with open(cnv_list_file, 'w') as f:
        for i, cnv_region in enumerate(cnv_matrix):
            line = cnv_region[0][0] + '\t' \
                   + `cnv_region[0][1]` + '\t' \
                   + `cnv_region[-1][2]` + '\t' \
                   + `mask[i]` + '\n'
            f.write(line)
    log("CNV list saved to " + cnv_list_file)

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
    fileio.mergeReads(tmp_dir,output_dir)

    log("cleaning temporary files ..")
    fileio.clean(tmp_dir)

    log("simulation completed. find results in " + output_dir)

if __name__ == '__main__':
    main()