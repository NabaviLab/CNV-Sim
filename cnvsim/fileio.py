#!/usr/bin/python

__author__ = 'Abdelrahman Hosny'

import subprocess
import os

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
            exon = (chromosome, int(start), int(end))
            exons.append(exon)
    return exons

def prepareTargetFile(target_file):
    '''
    sort and merge targets in the target file and writes it to '.sorted.merged'
    :param target_file: target exons in BED format
    :return: sorted and merged file name
    '''
    sorted_file = target_file + ".sorted"
    merged_file = sorted_file + ".merged"

    with open(sorted_file, "w") as f:
        subprocess.call(["sort", "-k1,1", "-k2,2n", target_file], stdout=f)
    with open(merged_file, "w") as f:
        subprocess.call(["bedtools", "merge", "-i", sorted_file], stdout=f)

    return merged_file

def mergeWessimReads(tmp_dir, output_dir):
    '''
    merges the base reads with normal and cnv
    :return: null
    '''
    base_file_1 = os.path.join(tmp_dir, "base_1.fastq.gz")
    base_file_2 = os.path.join(tmp_dir, "base_2.fastq.gz")
    normal_file_1 = os.path.join(tmp_dir, "control_1.fastq.gz")
    normal_file_2 = os.path.join(tmp_dir, "control_2.fastq.gz")
    cnv_file_1 = os.path.join(tmp_dir, "cnv_1.fastq.gz")
    cnv_file_2 = os.path.join(tmp_dir, "cnv_2.fastq.gz")

    with open(os.path.join(output_dir, "control_1.fastq.gz"), "w") as f:
        subprocess.call(["cat", base_file_1, normal_file_1], stdout=f)
    with open(os.path.join(output_dir, "control_2.fastq.gz"), "w") as f:
        subprocess.call(["cat", base_file_2, normal_file_2], stdout=f)
    with open(os.path.join(output_dir, "cnv_1.fastq.gz"), "w") as f:
        subprocess.call(["cat", base_file_1, cnv_file_1], stdout=f)
    with open(os.path.join(output_dir, "cnv_2.fastq.gz"), "w") as f:
        subprocess.call(["cat", base_file_2, cnv_file_2], stdout=f)

def mergeARTReads(tmp_dir, output_dir):
    '''
    merges the base reads with normal and cnv
    :return: null
    '''
    base_file = os.path.join(tmp_dir, "base.fq")
    normal_file = os.path.join(tmp_dir, "control.fq")
    cnv_file = os.path.join(tmp_dir, "cnv.fq")

    with open(os.path.join(output_dir, "control.fq"), "w") as f:
        subprocess.call(["cat", base_file, normal_file], stdout=f)
    with open(os.path.join(output_dir, "cnv.fq"), "w") as f:
        subprocess.call(["cat", base_file, cnv_file], stdout=f)

def clean(tmp_dir):
    subprocess.call(["rm", "-rf", tmp_dir])
