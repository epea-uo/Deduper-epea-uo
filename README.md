# PCR Deduper

## Overview:

Given a sorted SAM file of uniquely mapped reads, and a text file containing the known UMIs, this program will remove all PCR duplicates (retain only a single copy of each read).
- It accounts for: 
    - all possible CIGAR strings (including adjusting for soft clipping, etc.)
    - Strand
    - Single-end reads
    - Known UMIs

## How to run:

Command line options: 
    
    - -f,--file: designates absolute file path to sorted sam file
    
    - -o, --outfile: designates absolute file path to deduplicated sam file
    
    - -u, --umi: designates file containing the list of UMIs
    
    - -h, --help: prints a help message
    
Example command:
    
    ./pearce_deduper.py -u <UMI.file> -f <in.sam> -o <out.sam>


