#!/usr/bin/env python3
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("truthcov")
parser.add_argument("querycov")
args = parser.parse_args()

truth_tp = 0
truth_fn = 0
with open(args.truthcov) as infile:
    for line in infile:
        cols = line.strip().split("\t")
        cov = float(cols[5])
        if cov >= 0.5:
            truth_tp += 1
        else:
            truth_fn += 1

query_tp = 0
query_fp = 0
with open(args.querycov) as infile:
    for line in infile:
        cols = line.strip().split("\t")
        cov = float(cols[5])
        if cov >= 0.5:
            query_tp += 1
        else:
            query_fp += 1

recall = truth_tp / (truth_tp + truth_fn)*100
precision = query_tp / (query_tp + query_fp)*100
f1 = 2*recall*precision / (recall + precision)
print(f"Recall: {recall:.2f}%")
print(f"Precision: {precision:.2f}%")
print(f"F1: {f1:.2f}%")
