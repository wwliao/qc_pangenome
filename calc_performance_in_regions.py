#!/usr/bin/env python3
import os
import argparse
from collections import defaultdict

from intervaltree import IntervalTree
from cyvcf2 import VCF

def get_parser():
    parser = argparse.ArgumentParser(description="Calculate structural variation calling performance in stratification regions based on Truvari outputs")
    parser.add_argument("-d", "--dir", default=".", help="path to Truvari output direcotry (default is current directory)")
    parser.add_argument("bedfile")
    return parser

def build_itree(bedfile):
    regions = defaultdict(IntervalTree)
    with open(bedfile) as infile:
        for line in infile:
            cols = line.strip().split("\t")
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            regions[chrom].addi(start, end)
    return regions

def get_overlap_count(regions, vcffile):
    count = 0
    for variant in VCF(vcffile):
        chrom = variant.CHROM
        start = variant.start
        end = variant.end - 1
        if regions[chrom].overlaps(start) and regions[chrom].overlaps(end):
            count += 1
    return count

def main():
    parser = get_parser()
    args = parser.parse_args()
    regions = build_itree(args.bedfile)
    truth_tp = get_overlap_count(regions, os.path.join(args.dir, "tp-base.vcf"))
    query_tp = get_overlap_count(regions, os.path.join(args.dir, "tp-call.vcf"))
    truth_fn = get_overlap_count(regions, os.path.join(args.dir, "fn.vcf"))
    query_fp = get_overlap_count(regions, os.path.join(args.dir, "fp.vcf"))

    truth_total = truth_tp + truth_fn
    query_total = query_tp + query_fp
    recall = truth_tp/truth_total
    precision = query_tp/query_total
    f1 = 2*recall*precision/(recall + precision)
    print("#Truth TP\tQuery TP\tTruth FN\tQuery FP\tRecall\tPrecision\tF1")
    print(f"{truth_tp}\t{query_tp}\t{truth_fn}\t{query_fp}\t{recall:.6f}\t{precision:.6f}\t{f1:.6f}")

if __name__ == "__main__":
    exit(main())
