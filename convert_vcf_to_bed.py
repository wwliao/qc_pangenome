#!/usr/bin/env python3
import re
import argparse
from os.path import splitext, basename

parser = argparse.ArgumentParser()
parser.add_argument("vcffile")
args = parser.parse_args()

prefix = splitext(basename(args.vcffile))[0]
with open(args.vcffile) as infile, open(f"{prefix}.bed", "w") as outfile:
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
            outfile.write(f"{chrom}\t{start}\t{end}\t{id}\t{qual}\t{ref}\t{alt}\t{filter}\t{info}\t{format}\t{sample}\n")
