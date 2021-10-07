#!/usr/bin/env python3
import re
import argparse
from os.path import basename, splitext
from pyfaidx import Fasta

parser = argparse.ArgumentParser()
parser.add_argument("vcfheader")
parser.add_argument("reference")
parser.add_argument("stable_traversals")
args = parser.parse_args()

ref = Fasta(args.reference)

prefix = splitext(basename(args.stable_traversals))[0]
with open(args.stable_traversals) as infile, open(f"{prefix}.vcf", "w") as outfile, open(f"{prefix}.bug.vcf", "w") as bug_outfile:
    with open(args.vcfheader) as vcfheader:
        for line in vcfheader:
            outfile.write(line)

    for i, line in enumerate(infile):
        if i > 0:
            cols = line.strip().split("\t")
            id = cols[2]
            lv = cols[3]
            ps = cols[4]
            #ref_at = cols[5]
            #alt_at = cols[6]
            #at = f"{ref_at},{alt_at}"
            segs = cols[8].split(",")
            variant_indexes = []
            for idx, seg in enumerate(segs):
                if not seg.startswith(">chr"):
                    variant_indexes.append(idx)

            for j, idx in enumerate(variant_indexes):
                j += 1
                if not segs[idx].startswith("<chr"):
                    m1 = re.search("([><])(chr.+):(\d+)-(\d+)", segs[idx-1])
                    m2 = re.search("([><])(chr.+):(\d+)-(\d+)", segs[idx+1])
                    direct1 = m1[1]
                    chrom1 = m1[2]
                    direct2 = m2[1]
                    chrom2 = m2[2]
                    if direct1 == "<" and direct2 == "<":
                        start1 = int(m2[3])
                        end1 = int(m2[4])
                        start2 = int(m1[3])
                        end2 = int(m1[4])
                    else:
                        start1 = int(m1[3])
                        end1 = int(m1[4])
                        start2 = int(m2[3])
                        end2 = int(m2[4])
                    if chrom1 == chrom2:
                        chrom = chrom1
                        start = end1
                        end = start2
                        seq = segs[idx]
                        ref_len = end - start
                        if segs[idx] != "*":
                            alt_len = len(segs[idx])
                            vlen = alt_len - ref_len
                            if vlen == 0:
                                if alt_len == 1:
                                    vtype = "SNP"
                                else:
                                    vtype = "MNP"
                            elif ref_len == 0:
                                vtype = "INS"
                            elif vlen > 0:
                                vtype = "INS"
                            elif vlen < 0:
                                vtype = "DEL"
                            else:
                                vtype = "UNK"
                        else:
                            vlen = -1*ref_len
                            vtype = "DEL"

                        if vtype == "SNP":
                            ref_allele = ref[chrom][start].seq
                            alt_allele = seq
                            start += 1
                        elif vtype == "MNP":
                            ref_allele = ref[chrom][start:end].seq
                            alt_allele = seq
                            start += 1
                        elif vtype == "INS":
                            if end >= start:
                                ref_allele = ref[chrom][start-1:end].seq
                                alt_allele = ref_allele[0] + seq

                            else:
                                ref_allele = "N"
                                alt_allele = seq
                        elif vtype == "DEL":
                            if end > start:
                                ref_allele = ref[chrom][start-1:end].seq
                                alt_allele = ref_allele[0]
                            # Need to think about it more carefully
                            else:
                                ref_allele = "N"
                                alt_allele = seq
                        elif vtype == "UNK":
                            if end > start:
                                ref_allele = ref[chrom][start-1:end].seq
                                alt_allele = ref_allele[0] + seq
                            # Need to think about it more carefully
                            else:
                                ref_allele = "N"
                                alt_allele = seq

                        if len(variant_indexes) > 1:
                            if vtype in ["SNP", "MNP"]:
                                outfile.write(f"{chrom}\t{start}\t{id}-{j}\t{ref_allele}\t{alt_allele}\t60\t.\t.\tGT\t0/1\n")
                            elif vtype == "UNK":
                                if ref_allele == "N":
                                    bug_outfile.write(f"{chrom}\t{start}\t{id}-{j}\t{ref_allele}\t{alt_allele}\t60\t.\tSVTYPE={vtype};END={end};SVLEN={vlen}\tGT\t0/1\n")
                                else:
                                    outfile.write(f"{chrom}\t{start}\t{id}-{j}\t{ref_allele}\t{alt_allele}\t60\t.\tSVTYPE={vtype};END={end};SVLEN={vlen}\tGT\t0/1\n")
                            else:
                                if vtype in ["DEL", "INS"] and ref_allele == "N":
                                    bug_outfile.write(f"{chrom}\t{start}\t{id}-{j}\t{ref_allele}\t{alt_allele}\t60\t.\tSVTYPE={vtype};END={end};SVLEN={vlen}\tGT\t0/1\n")
                                else:
                                    outfile.write(f"{chrom}\t{start}\t{id}-{j}\t{ref_allele}\t{alt_allele}\t60\t.\tSVTYPE={vtype};END={end};SVLEN={vlen}\tGT\t0/1\n")
                        else:
                            if vtype in ["SNP", "MNP"]:
                                outfile.write(f"{chrom}\t{start}\t{id}\t{ref_allele}\t{alt_allele}\t60\t.\t.\tGT\t0/1\n")
                            elif vtype == "UNK":
                                if ref_allele == "N":
                                    bug_outfile.write(f"{chrom}\t{start}\t{id}-{j}\t{ref_allele}\t{alt_allele}\t60\t.\tSVTYPE={vtype};END={end};SVLEN={vlen}\tGT\t0/1\n")
                                else:
                                    outfile.write(f"{chrom}\t{start}\t{id}\t{ref_allele}\t{alt_allele}\t60\t.\tSVTYPE={vtype};END={end};SVLEN={vlen}\tGT\t0/1\n")
                            else:
                                if vtype in ["DEL", "INS"] and ref_allele == "N":
                                    bug_outfile.write(f"{chrom}\t{start}\t{id}\t{ref_allele}\t{alt_allele}\t60\t.\tSVTYPE={vtype};END={end};SVLEN={vlen}\tGT\t0/1\n")
                                else:
                                    outfile.write(f"{chrom}\t{start}\t{id}\t{ref_allele}\t{alt_allele}\t60\t.\tSVTYPE={vtype};END={end};SVLEN={vlen}\tGT\t0/1\n")
                    else:
                        print(f"{chrom1} != {chrom2}")
                else:
                    m = re.search("[><](chr.+):(\d+)-(\d+)", segs[idx])
                    chrom = m[1]
                    start = int(m[2])
                    end = int(m[3])
                    vlen = end - start + 1
                    vtype = "INV"
                    ref_allele = ref[chrom][start-1:start]
                    alt_allele = "<INV>"
                    if len(variant_indexes) > 1:
                        outfile.write(f"{chrom}\t{start}\t{id}-{j}\t{ref_allele}\t{alt_allele}\t60\t.\tSVTYPE={vtype};END={end};SVLEN={vlen}\tGT\t0/1\n")
                    else:
                        outfile.write(f"{chrom}\t{start}\t{id}\t{ref_allele}\t{alt_allele}\t60\t.\tSVTYPE={vtype};END={end};SVLEN={vlen}\tGT\t0/1\n")
