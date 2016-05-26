__author__ = 'Abdelrahman Hosny'

def readGenome(filename):
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

genome_file = 'data/chr1.fa'
target_file = 'data/chr1-sample-target.bed'

modified_genome_file = 'data/chr1.sim.fa'
modified_target_file = 'data/chr1-sample-target.sim.bed'

exons = []
with open(target_file, 'r') as tf:
    for line in tf:
        chromosome, start, end, description, score, strand = line.strip().split('\t')
        exon = {"chromosome": chromosome, \
                "start": int(start), \
                "end": int(end), \
                "description": description, \
                "score": score, \
                "strand": strand}
        exons.append(exon)

# just for testing <adding one amplification only>
# this portion should be changed with our smart algorithm
header, genome = readGenome(genome_file)
last_exon = exons[-1]
modified_genome = genome[:last_exon["end"]] + genome[last_exon["start"]:last_exon["end"]] + genome[last_exon["end"]:]
exons[-1]["end"] = last_exon["end"] + (last_exon["end"] - last_exon["start"])

print len(genome)
print len(modified_genome)
print len(modified_genome) - len(genome)

# now write results to files
with open(modified_genome_file, 'w') as f:
    f.write(header)
    line_width = 50
    chunks = [modified_genome[i:i+line_width]+'\n' for i in range(0, len(genome), line_width)]
    f.writelines(chunks)

with open(modified_target_file, 'w') as f:
    targets = []
    for ex in exons:
        target = ex["chromosome"] + '\t' \
                 + str(ex["start"]) + \
                 '\t' + str(ex["end"]) + \
                 '\t' + ex["description"] + \
                 '\t' + ex["score"] + \
                 '\t' + ex["strand"] + '\n'
        targets.append(target)
    f.writelines(targets)