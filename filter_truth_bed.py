#!/usr/bin/env python3
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-n", "--numcallers", metavar="INT", type=int, default=1, help="Consider only variants called by at least INT callers")
parser.add_argument("-t", "--numcallertypes", metavar="INT", type=int, default=1, help="Consider only variants called by at least INT caller types")
parser.add_argument("--has-illumina", action="store_true", help="Consider variants called by illumina-based caller")
parser.add_argument("--has-longread", action="store_true", help="Consider variants called by longread-based caller")
parser.add_argument("--has-assembly", action="store_true", help="Consider variants called by assembly-based caller")
parser.add_argument("-o", "--output", required=True)
parser.add_argument("bedfile")
args = parser.parse_args()

with open(args.bedfile) as infile, open(args.output, "w") as outfile:
    for line in infile:
        if not line.startswith("#"):
            cols = line.strip().split("\t")
            ncallers = int(cols[8])
            ncallertypes = int(cols[10])
            has_ill = bool(int(cols[11]))
            has_lr = bool(int(cols[12]))
            has_asm = bool(int(cols[13]))
            if ncallers >= args.numcallers and ncallertypes >= args.numcallertypes:
                flags = 0
                if args.has_illumina:
                    if has_ill:
                        flags += 1
                else:
                    flags += 1
                if args.has_longread:
                    if has_lr:
                        flags += 1
                else:
                    flags += 1
                if args.has_assembly:
                    if has_asm:
                        flags += 1
                else:
                    flags += 1

                if flags == 3:
                    outfile.write(line)
