#!/usr/bin/env python3
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--coverage", type=float, default=0.5)
parser.add_argument("truth")
parser.add_argument("query")
parser.add_argument("truthcov")
parser.add_argument("querycov")
args = parser.parse_args()

total_truth = 0
with open(args.truth) as infile:
    for line in infile:
        total_truth += 1

total_query = 0
with open(args.query) as infile:
    for line in infile:
        total_query += 1

truth_tp = 0
with open(args.truthcov) as infile:
    for line in infile:
        cols = line.strip().split("\t")
        cov = float(cols[5])
        if cov >= args.coverage:
            truth_tp += 1

truth_fn = total_truth - truth_tp

query_tp = 0
with open(args.querycov) as infile:
    for line in infile:
        cols = line.strip().split("\t")
        cov = float(cols[5])
        if cov >= args.coverage:
            query_tp += 1

query_fp = total_query - query_tp

recall = truth_tp / (truth_tp + truth_fn)*100
precision = query_tp / (query_tp + query_fp)*100
f1 = 2*recall*precision / (recall + precision)
print(f"Recall: {truth_tp}/({truth_tp} + {truth_fn}) = {recall:.2f}%")
print(f"Precision: {query_tp}/({query_tp} + {query_fp}) = {precision:.2f}%")
print(f"F1: {f1:.2f}%")
