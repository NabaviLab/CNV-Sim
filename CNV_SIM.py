#! usr/bin/python

import sys
import random
import re
import numpy as np
from subprocess import call
import itertools
import csv


# parser = argparse.ArgumentParser()
# parser.add_argument("Reference", help="You must provide a reference genome")
# args = parser.parse_args()
# print args.Reference
#
# chromosome1 = open(sys.argv[1], "r")

chromosome1 = open("chr1.fa", "r")
regions = open("S04380219_Regions.bed", "r")
# out_cnv_chr1_A1 = open("Cnv_Chr1_A1.fa", "w")
# out_cnv_chr1_A2 = open("Cnv_Chr1_A2.fa", "w")
# cnv_list = open("Cnv_list.csv", "wt")
new_regions = open("target_updates.bed", "w")
modified_regions = open("Modified_Regions.bed", "w")
#Ori_chr1 = open("Chr1_Regions.bed", "w")

head_chr = chromosome1.readline()
head_1_region = regions.readline()
head_2_region = regions.readline()

chr1 = ""
for line in chromosome1:
    lines = line.strip()
    chr1  += lines
    
chr1_A1 = chr1
chr1_A2 = chr1

print len(chr1_A1), "length of A1"
print len(chr1_A2), "length of A2"
####### Extract exons on Chr1 ########
#chr1_region = open("chr1_target.txt", "w")

exons = []
for line in regions.readlines():
    if re.findall('\\bchr1\\b', line, flags=re.IGNORECASE):
        #exon = re.split(',|\t', line)
        #exon = line.rstrip('\n')
        exon = line.split("\t")
        exons.append(exon)
        
print len(exons), "Total number of chr1 exons"

# for thing in exons:
#     if "-\n" in thing[3]:
#         exons.remove(thing) 



# print>>Ori_chr1, head_1_region.strip("\n")
# print>>Ori_chr1, head_2_region.strip("\n")
#
# for item in exons:
#     items = [line.rstrip('\n') for line in item]
#     print>>Ori_chr1, "\t".join(items)


# print>>cnv_list, head



for i in range(1,11):
    head = "chr" + "," + "start" + "," + "end" + "," + "copy" + "," + "type"
    cnv_out = []
    out_cnv_chr1_A1 = open("%s_Cnv_Chr1_A1.fa"% (i), "w")
    out_cnv_chr1_A2 = open("%s_Cnv_Chr1_A2.fa"% (i), "w")
    cnv_list = open("%s_Cnv_list.csv"% (i), "wt")
    print>>cnv_list, head

    ####### Selecting Exons, generating N non-overlapping cnvs ###################

    def rand_exon_pool(exons):
        #final_cnv = []
        cnv_start = np.arange(10, len(exons)-50, 10)
        cnv_long = np.arange(2, 15)
        list_length = 0
        cnv_length = []
        #while len(cnv_length) < 3:
        final_cnv = []
        cnv_S = random.sample(cnv_start, 500)
        cnv_S = sorted(cnv_S)
        #cnv_length = []
        for item in cnv_S:
            # print item, "start"
            # print
            cnv_leng = random.choice(cnv_long)
            # print cnv_leng, "length"
            # print
            end = item+cnv_leng
            # print end, "stop"
            if len(cnv_length) ==0:
                cnv_length.append((item,end))
            if len(cnv_length) > 0:
                if item > cnv_length[-1][1]:
                    cnv_length.append((item,end))
        exon_lines = random.sample(cnv_length, 100) #### Change the number here to get the total number of Simulated CNVs
        exon_lines = sorted(exon_lines)
        return exon_lines

    exon_extract = rand_exon_pool(exons)


    ##### Randomize Amplifications and Deletions #######  Define your amps vs dels ratio here ######

    # del_cnv = random.sample(exon_extract, 50)  ### % Deletions
    # del_cnv.sort()

    del_num = random.choice([86,80,74,63,58])

    del_cnv = exon_extract[del_num:]
    del_cnv.sort()
    print del_cnv, "del items", len(del_cnv)
    print

    for item in del_cnv:
        exon_extract.remove(item)  ####  % Amplifications

    print len(exon_extract), "Amplification"
    #print exon_extract, "left for amp", len(exon_extract)

    ###################################################################

    ########## Exons Amplifications ####################

    seq_chunk = []
    exom = []
    exon_pos = []
    exom_update = []

    for pos in exon_extract:
        exon_chunk = exons[pos[0]:pos[1]]
        #if len(exon_chunk) < 2:
        #exon_chu = [i.split('\n')[0] for i in exon_chunk]
        #print len(exon_chunk), "ori"
        exom.append(exon_chunk)
        #print len(exom), "individual"
        #for line in exon_chunk:
        cnv_start = exon_chunk[0]
        cnv_end = exon_chunk[-1]
        # print cnv_start, "start line"
        # print cnv_end, "end line"
        exon_pos.append((cnv_start[1],cnv_end[2]))
        ex_pos = [cnv_start[1],cnv_end[2]]
        #print ex_pos
        start=int(ex_pos[0])
        end=int(ex_pos[1])
        #print start, end
        copy = random.choice([1,2,3,4,5])
        #print copy, "copy"
        chunk = chr1[start:end]
        if copy == 1:  ## 1.5
            chr1_A1 += chunk
            new_line = "1", str(start),str(end),str(copy),"Amplification"
            cnv_out.append(new_line)
            # print len(exon_chunk), "ori chunk length"
            # print exon_chunk
            exoms = exon_chunk
            # print len(exoms), "length of amp exom", copy
            # print exoms
            exom_update.append(exoms)
            #print>>cnv_list, new_line
        elif copy == 2:
            chr1_A1 += chunk
            chr1_A2 += chunk
            new_line = "1", str(start),str(end),str(copy),"Amplification"
            cnv_out.append(new_line)
            # print len(exon_chunk), "ori"
            # print exon_chunk
            exoms = exon_chunk * copy
            # print len(exoms), "amp"
            # print exoms
            exom_update.append(exoms)
            #print>>cnv_list, new_line
        elif copy == 3:
            chr1_A1 += chunk
            chr1_A2 += chunk*2
            new_line = "1", str(start),str(end),str(copy),"Amplification"
            cnv_out.append(new_line)
            exoms = exon_chunk * copy
            exom_update.append(exoms)
            #print>>cnv_list, new_line
        elif copy == 4:
            chr1_A1 += chunk*2
            chr1_A2 += chunk*2
            new_line = "1", str(start),str(end),str(copy),"Amplification"
            cnv_out.append(new_line)
            exoms = exon_chunk * copy
            exom_update.append(exoms)
            #print>>cnv_list, new_line
        elif copy == 5:
            chr1_A1 += chunk*2
            chr1_A2 += chunk*3
            new_line = "1", str(start),str(end),str(copy),"Amplification"
            cnv_out.append(new_line)
            exoms = exon_chunk * copy
            exom_update.append(exoms)
            #print>>cnv_list, new_line

    #     print len(chr1_A1),"A1 Amp", pos, ex_pos
    #     print len(chr1_A2),"A2 Amp", pos, ex_pos
    #
    # print
    # print len(chr1_A1), "Allel 1 after Amplifications"
    # print len(chr1_A2), "Allel 2 after Amplifications"
    # print len(exom_update), "Amp chunk length"

    #################################################################


    #######  Exon Deletions #########################


    del_seq_chunk = []
    del_exom = []
    del_exon_pos = []
    L = []

    for pos in del_cnv:
        exon_chunk = exons[pos[0]:pos[1]]
        #print len(exons), "before", len(exon_chunk)
        for lines in exon_chunk:
            if lines in exons:
                #print lines
                exons.remove(lines)
        #print len(exons), "after"
        #print exon_chunk, "original chunk"
        #exon_chu = [i.split('\n')[0] for i in exon_chunk]
        #print len(exon_chunk), "ori"
        #print exon_chunk, "starting chunk"
        try:
            cnv_start = exon_chunk[0]
            cnv_end = exon_chunk[-1]
        except:
            print exon_chunk
        # print cnv_start, "start line"
        # print cnv_end, "end line"
        del_exon_pos.append((cnv_start[1],cnv_end[2]))
        ex_pos = [cnv_start[1],cnv_end[2]]
        #for posi in del_exon_pos:
        start=int(ex_pos[0])
        end=int(ex_pos[1])
        #print start, end, "original"
        length = end-start
        L.append(length)
        #print L, "list of deletion size"
        if len(L) >1:
            length = sum(L[0:len(L)-1])
            #print length, "total deletion size"
            start = start - length
            end = end - length
            for lines in exon_chunk:
                #print lines, "original exon position"
                #print
                new_start = int(lines[1]) - length
                new_end = int(lines[2]) - length
                lines[1] = str(new_start)
                lines[2] = str(new_end)
                #print lines, "updated position"
        #print exon_chunk, "modified chunk", (start,end), "new corr"
        del_exom.append(exon_chunk)
        #print start, end, "new coordinates"
        #dele = random.choice([1,2])
        chunk = chr1[start:end]
        # print len(chunk)
        # print
        #if dele == 2:
            # chr1_A1 = re.sub(chunk, '', chr1_A1, flags=re.IGNORECASE)
            # chr1_A2 = re.sub(chunk, '', chr1_A2, flags=re.IGNORECASE)
            # # print len(chr1_A1), "length of A1"
            # # print len(chr1_A2), "length of A2"
            # new_line = "1", str(start), str(end), str(dele), "Homozygous Deletion"
            # cnv_out.append(new_line)
            #print>>cnv_list, new_line
        #elif dele ==1:
            #chr1_A1 = re.sub(chunk, '', chr1_A1, flags=re.IGNORECASE)
        chr1_A2 = re.sub(chunk, '', chr1_A2, flags=re.IGNORECASE)
            # print len(chr1_A1), "length of A1 hetero"
        new_line = "1", str(start), str(end), "1", "Deletion"
        cnv_out.append(new_line)
            #print>>cnv_list, new_line

        # print len(chr1_A1), "A1 dele", pos, ex_pos
        # print len(chr1_A2), "A2 dele", pos, ex_pos

    # print len(del_exom)
    # print
    # print del_exon_pos
    #print L
    print len(chr1_A1), "Allel 1 after deletion"
    print len(chr1_A2), "Allel 2 after deletion"
    # print len(del_exom), "del chunk length"

    ##################################################


    ########### Update Output Target regions #############

    print>>modified_regions, head_1_region.strip("\n")
    print>>modified_regions, head_2_region.strip("\n")

    for item in exons:
        items = [line.rstrip('\n') for line in item]
        print>>modified_regions, "\t".join(items)


    New_exom = del_exom + exom_update
    for item in New_exom:
        for ele in item:
            #print ele
            eles = [line.rstrip('\n') for line in ele]
            print>>new_regions, "\t".join(eles)

    ############ Output Cnv list ##############################

    # new_rand_list = random.sample(cnv_out, 100)
    # new_rand_list = sorted(new_rand_list,key=lambda x: x[1])
    #
    # print new_rand_list

    writer = csv.writer(cnv_list)
    #writer.writerow( ("Chr","Start","End","Copy","Type") )
    for row in cnv_out:
        writer.writerow(row)


    cnv_list.close()

    # for item in cnv_out:
    #     print>>cnv_list, "".join(item)


    new_chr1_A1 = "".join(chr1_A1)
    new_chr1_A2 = "".join(chr1_A2)



    ##########  Outout New Sequence   ###############################################
    head = ">chr1"
    print >>out_cnv_chr1_A1, head
    print >>out_cnv_chr1_A2, head

    n=50
    split_list = [new_chr1_A1[i:i+n] for i in range(0, len(new_chr1_A1), n)]
    for item in split_list:
        print >>out_cnv_chr1_A1, item
    print "Finished Processing Allel 1"

    n=50
    split_list = [new_chr1_A2[i:i+n] for i in range(0, len(new_chr1_A2), n)]
    for item in split_list:
        print >>out_cnv_chr1_A2, item
    print  "Finished Processing Allel 2"
    ###########################################################


    # ###### Generate HTML table for CNV list  #############
    #
    # reader = csv.reader(open("%s_Cnv_list.csv"% (i), "rt"))
    # f_html = open("%s_CNV_list.html"% (i),"w")
    #
    #
    #
    # print>>f_html, """<html>
    # <head>
    #  <title>Simulated CNVs</title>
    # </head>
    # <body>
    # <table border="1">"""
    #
    # print>>f_html, "<tr><th>Chr</th><th>Start</th><th>End</th><th>Copy</th><th>Type</th></tr>"
    #
    # for row in reader:
    #     print>>f_html, "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>" % (row[0], row[1], row[2], row[3], row[4])
    #
    # print>>f_html, """</table></body></html>"""

# ##### Alternative way #################

#rownum=0
# f_html.write('<table>')
# f_html.write('<title><List of Simulated CNVs></title>')

# for row in reader:
#     if rownum == 0:
#         f_html.write('<tr>');
#         for column in row:
#             f_html.write('<th>' + column + '</th>');
#         f_html.write('</tr>')
#     else:
#         f_html.write('<tr>')
#         for column in row:
#             f_html.write('<td>' + column + '</td>')
#         f_html.write('</tr>')
#     rownum += 1
#
#
# f_html.write('</table>')

#######################################











chromosome1.close()
regions.close()
out_cnv_chr1_A1.close()
out_cnv_chr1_A2.close()
#f_html.close()
modified_regions.close()
#Ori_chr1.close()