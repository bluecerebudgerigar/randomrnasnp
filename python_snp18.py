#!/usr/bin/env python

from Bio import SeqIO
import random
import sys
import subprocess
import re
import os
import signal

### Declaration of values based on Paramter values ###
record_dict = SeqIO.index(sys.argv[1], "fasta")
iterations = int(sys.argv[2])
p_value_cutoff = float(sys.argv[3])
transcripts_cutoff = float(sys.argv[4])
quick_mode = sys.argv[5]

### Declare Global Value ###

snp_input_list=[]
snp_input_value=""
score = 0
command=["RNAsnp", "-f","test_fasta","-s","test_snp","-m","2"]
quick_mode_score =  float(iterations - transcripts_cutoff)

    

sequence_length=0

### random_snps is a function that look at choose a
def random_snps(cycle_number):
    snp_input_list=[]
    snp_input=""
    if cycle_number > 1:
        snp_pos=random.sample(range(1,len(sequence)), cycle_number)
        snp_pos.sort()
        for i in snp_pos:
            base=sequence[i-1]
            snp_base=random.choice([x for x in "ATGC" if x != base])
            snp_input = "%s%s%s" % (base, str(i), snp_base)
            snp_input_list.append(snp_input)
            snp_input ="-".join(snp_input_list)
        return snp_input
    else:
        snp_pos=random.randint(1,len(sequence))
        base=sequence[snp_pos-1]
        snp_base=random.choice([x for x in "ATGC" if x != base])
        snp_input = base + str(snp_pos) + snp_base
        return snp_input
        
def test_p(value):
    global hits
    global pos_hits
    handle = open("mode2_output", "r")
    for lines in handle:
        print lines
        entries = re.match('[SNP]', lines)
        if entries:
            pass
        else:
            columns = lines.split('\t')
            columns = [col.strip() for col in columns]
            if float(columns[-1]) > value:
                hits = hits + 1
            elif float(columns[-1]) < value:
                pos_hits = pos_hits + 1 
    return pos_hits
    return hits

### Write Fasta and Write Results File.         

def write_file(object_name, file_name):
    handle = open(file_name, "w+")    
    handle.write(object_name)
    handle.close()
    
def results_file(type_results, delim):
    handle = open("results.txt", "a")
    handle.write(str(type_results) + delim)
    handle.close()

    
header_results="sequence_id\tsequence_length\tcycle"
results_file(header_results,"\n")


    


for keys in record_dict:
    results_file(record_dict[keys].id, "\t")
    i = 0 #keep track of number of snps
    i2 = 0 # master control for snp addition
    #exhaust = 0 might be required to when sequence are too short.
    sequence = str(record_dict[keys].seq)
    sequence=str.upper(sequence)
    sequence = str(sequence)
    write_file(sequence, "test_fasta")
    sequence_length=len(sequence)
#### Generate SNPs ####
    while i2 < 1:
        n = 0
        hits = 0
        pos_hits = 0
        snp_list=[]
        while n < iterations:
            snp_input=random_snps(i)
            if snp_input not in snp_list:
                snp_list.append(snp_input)
                n = n + 1
            else:
                continue
        for snps in snp_list:
            write_file(snps, "test_snp")
            subprocess.call(command, stdout=open("mode2_output", "w+"), stderr=open("/dev/null"))
            test_p(p_value_cutoff)
            print pos_hits 
            if quick_mode == "yes":           
                if pos_hits >= transcripts_cutoff:
                    i2 += 1
                    break
                if hits > quick_mode_score:
                    break
        if quick_mode == "no":
            if pos_hits > transcripts_cutoff:
                i2 += 1
        i += 1
    results_file(sequence_length, "\t")        
    results_file(i, "\n")
    print("run completed")