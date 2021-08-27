#!/usr/bin/env python3
import re
import argparse
from cyvcf2 import VCF

parser = argparse.ArgumentParser()
parser.add_argument("vcffile", help="Single-sample VCF")
args = parser.parse_args()

prefix = re.search("(.+)\.vcf(?:\.gz)*", args.vcffile)[1]
with open(f"{prefix}.allele_traversals.txt", "w") as outfile:
    outfile.write("#CHROM\tPOS\tID\tLV\tPS\tREF_AT\tALT_AT\n")
    for variant in VCF(args.vcffile):
        chrom = variant.CHROM
        pos = variant.POS
        id = variant.ID
        lv = variant.INFO.get("LV")
        ps = "."
        if lv > 0:
            ps = variant.INFO.get("PS")
        gt = set(variant.genotypes[0][:2])
        if (1 in gt) or (2 in gt):
            ats = variant.INFO.get("AT").split(",")
            ref_at = ats[0]
            alt_ats = []
            for idx in sorted(gt):
                if idx != 0:
                    alt_ats.append(ats[idx])

            if len(alt_ats) == 1:
                alt_at = alt_ats[0]
                outfile.write(f"{chrom}\t{pos}\t{id}\t{lv}\t{ps}\t{ref_at}\t{alt_at}\n")
            else:
                for i, alt_at in enumerate(alt_ats):
                    i += 1
                    outfile.write(f"{chrom}\t{pos}\t{id}_{i}\t{lv}\t{ps}\t{ref_at}\t{alt_at}\n")
