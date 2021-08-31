#!/usr/bin/env python3
import re
import argparse
from os.path import basename, splitext

from pyfaidx import Fasta

def revcomp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sample", required=True)
parser.add_argument("reference")
parser.add_argument("vcffile")
args = parser.parse_args()

ref = Fasta(args.reference)

root, ext = splitext(basename(args.vcffile))
with open(args.vcffile) as infile, open(f"{root}.explicit{ext}", "w") as outfile:
    for line in infile:
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                outfile.write(line.strip())
                outfile.write(f"\tFORMAT\t{args.sample}\n")
            else:
                outfile.write(line)
        else:
            cols = line.strip().split("\t")
            chrom = cols[0]
            start = int(cols[1])
            id = cols[2]
            qual = cols[5]
            info = cols[7]
            svtype = re.search("SVTYPE=(\w+)", info)[1]
            if svtype in ["DEL", "INS", "INV"]:
                end = int(re.search("END=(\d+)", info)[1])
                
                if svtype == "DEL":
                    ref_allele = ref[chrom][start-1:end].seq
                    alt_allele = ref_allele[0]
                elif svtype == "INS":
                    ref_allele = ref[chrom][start-1].seq
                    alt_allele = "<INS>"
                    info = info.replace(f"END={end}", f"END={start}")
                elif svtype == "INV":
                    ref_allele = ref[chrom][start-1].seq
                    alt_allele = "<INV>"
            elif svtype == "BND":
                ref_allele = ref[chrom][start-1].seq
                alt_allele = cols[4].replace("N", ref_allele)
            outfile.write(f"{chrom}\t{start}\t{id}\t{ref_allele}\t{alt_allele}\t{qual}\t.\t{info}\tGT\t0/1\n")
