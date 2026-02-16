#!/usr/bin/env python

import argparse

def get_args():
    parser = argparse.ArgumentParser(description="A program to remove PCR duplicates from sequencing data")
    parser.add_argument("-f", "--file", help="Input sorted sam file", required=True)
    parser.add_argument("-o", "--outfile", help="Output sam file", required=True)
    parser.add_argument("-u", "--umi", help="File that contains a list of the unique UMIs, one per line", required=True)
    return parser.parse_args()

args = get_args()
infile=args.file
outfile=args.outfile
umi=args.umi

outlist = set()
umi_list = set()
current_chrom = ''
duplicate_count = 0
head_lines = 0
uniq_reads = 0
wrong_umis = 0

## Functions:

def extract_umi(info):
    long_umi = info[0]
    split_umi = long_umi.split(":")
    umi = split_umi[-1]
    return umi

def parse_cigar(cigar):
    cigar_parts = [] # list to hold tuples of (number, character) for each part of the CIGAR string
    current_num = ''
    for char in cigar:
        if char.isdigit():
            current_num += char
        else:
            if current_num != '':
                tuple = (current_num, char)
                cigar_parts.append(tuple)
                current_num = '' # reset current_num for the next part of the CIGAR string
    return cigar_parts

def get_pos(cigar_parts, pos, strand):
    if strand == "+": # if strand is positive, check for soft clipping at the start of the read
        if cigar_parts[0][1] =='S':
            start = int(pos) - int(cigar_parts[0][0]) # if there is soft clipping at the start of the read, do position - the number of bases that are soft clipped
        else:
            start = pos
    else:
        cigar_sum = 0
        for tuple in cigar_parts:
            if tuple[1] == "I": # if there is an insertion, skip it because it does not affect the position of the read on the reference genome
                continue
            cigar_sum += int(tuple[0])
        #print(info, "\n",cigar)
        if cigar_parts[0][1] == 'S':
            start = int(pos) + int(cigar_sum) - int(cigar_parts[0][0]) # if there is soft clipping at the start of the read, do position + the number of bases that are mapped to the reference genome - the number of bases that are soft clipped
            #print("soft clip\t",start)
        else:
            start = int(pos) + int(cigar_sum) # otherwise, do position + the number of bases that are mapped to the reference genome
            #print("no soft clip\t",start)
    return start

## Main Code:

with open(umi,"r") as umi_file:
    for line in umi_file:
        line = line.strip()
        umi_list.add(line)

with open(infile, "r") as infile, open(outfile, "w") as out:
    for line in infile:
        line = line.strip()
        if line.startswith("@"): # Automatically write any header line to output, count number of them
            out.write(f'{line}\n')
            head_lines +=1
        else:
            info = line.split("\t")
 			# Extract UMI
            umi = extract_umi(info)
            if umi not in umi_list: # Check if UMI is in the list of unique UMIs, if not add to wrong_umis count
                wrong_umis +=1
                continue
 			# Extract strand:
            flag = int(info[1])
            if ((flag & 16) ==16):
                strand = "-"
            else:
                strand = "+"
 			# CIGAR
            cigar = info[5]
            cigar_parts = parse_cigar(cigar)
 			# Position
            pos = int(info[3])
            start = get_pos(cigar_parts, pos, strand)
			# Chromosome
            chrom = info[2]
            key_tup = (umi, strand, start, chrom) # information to be used to determine if read is a duplicate or not
            # Check if still on the same chrom, if not, wipe the outlist and save memory
            if chrom != current_chrom:
                outlist = set()
                current_chrom = chrom
         	# Fill out output file
            if key_tup not in outlist:
                outlist.add(key_tup)
                out.write(f'{line}\n')
                uniq_reads+=1
            else:
                duplicate_count +=1

#print(outlist)

print(f'Number of header lines: {head_lines}')
print(f'Number of unique reads: {uniq_reads}')
print(f'Number of duplicated reads: {duplicate_count}')
print(f'Number of wrong UMIs: {wrong_umis}')
    
