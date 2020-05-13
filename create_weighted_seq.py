#!/usr/bin/env python3
import argparse
import math
import sys

parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, help='Input ancestralProbs')
args = parser.parse_args()

outFile = args.input + ".weighted.fasta"
outFH = open(outFile, "w")
seqs = dict()

with open(args.input, "r") as f:
    print("Weighting probabilities for {}\n".format(args.input))
    next(f)
    for line in f:
        e = line.split("\t")
        if e[0] not in seqs:
            seqs[e[0]] = ""
        if e[3] != 1 and e[4] != 1 and e[5] != 1 and e[6] != 1 and e[7] != 1:
            w = dict()
            w['A'] = 0.25/0.307 * float(e[3])
            w['C'] = 0.25/0.186 * float(e[4])
            w['G'] = 0.25/0.212 * float(e[5])
            w['T'] = 0.25/0.295 * float(e[6])
            w['-'] = float(e[7])
            final_base = " "
            final_freq = 0.0
            for base in w:
                if w[base] > final_freq:
                    final_freq = w[base]
                    final_base = base
            seqs[e[0]] += final_base
        else:
            seqs[e[0]] += e[2]
outContent = ""
for key in seqs:
    outContent += ">" + key + "\n" + seqs[key] + "\n"
outFH.write(outContent)
