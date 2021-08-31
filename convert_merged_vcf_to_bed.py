#!/usr/bin/env python3
import re
import argparse
from os.path import splitext, basename

parser = argparse.ArgumentParser()
parser.add_argument("vcffile")
args = parser.parse_args()

prefix = splitext(basename(args.vcffile))[0]
with open(args.vcffile) as infile, open(f"{prefix}.bed", "w") as outfile:
    outfile.write("#CHROM\tSTART\tEND\tID\tSVTYPE\tSVLEN\tREF\tALT\tNCALLERS\tCALLERS\tNCALLERTYPES\tHAS_ILL\tHAS_LR\tHAS_ASS\tHAS_PAT\tHAS_MAT\tHAS_DIP\tHITS_STR\tHITS_SEGDUP\tHITS_SAT\n")
    for line in infile:
        if not line.startswith("#"):
            cols = line.strip().split("\t")
            chrom = cols[0]
            start = int(cols[1]) - 1
            id = cols[2]
            ref = cols[3]
            alt = cols[4]
            qual = cols[5]
            filter = cols[6]
            info = cols[7]
            format = cols[8]
            sample = cols[9]
            end = int(re.search("END=(\d+)", info)[1])
            svtype = re.search("SVTYPE=(\w+)", info)[1]
            svlen = int(re.search("SVLEN=([-\d]+)", info)[1])
            ncallers = int(re.search("NCALLERS=(\d+)", info)[1])
            callers = re.search(";CALLERS=([,.\w]+)", info)[1]
            has_ill = int(re.search("HAS_ILL=(\d+)", info)[1])
            has_lr = int(re.search("HAS_LR=(\d+)", info)[1])
            has_ass = int(re.search("HAS_ASS=(\d+)", info)[1])
            ncallertypes = has_ill + has_lr + has_ass
            has_pat = int(re.search("HAS_PAT=(\d+)", info)[1])
            has_mat = int(re.search("HAS_MAT=(\d+)", info)[1])
            has_dip = int(re.search("HAS_DIP=(\d+)", info)[1])
            hits_str = int(re.search("HITS_STR=(\d+)", info)[1])
            hits_segdup = int(re.search("HITS_SEGDUP=(\d+)", info)[1])
            hits_sat = int(re.search("HITS_SAT=(\d+)", info)[1])
            outfile.write(f"{chrom}\t{start}\t{end}\t{id}\t{svtype}\t{ref}\t{alt}\t{svlen}\t{ncallers}\t{callers}\t{ncallertypes}\t{has_ill}\t{has_lr}\t{has_ass}\t{has_pat}\t{has_mat}\t{has_dip}\t{hits_str}\t{hits_segdup}\t{hits_sat}\n")
