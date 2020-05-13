#!/usr/bin/env python

import math
from Bio import SeqIO
import pandas as pd
import sys

import matplotlib.pyplot as plt

import logomaker as lm

###
dna = {'A': [1, 0, 0, 0, 0], 'C': [0, 1, 0, 0, 0], 'G': [0, 0, 1, 0, 0], 'T': [0, 0, 0, 1, 0], '-': [0, 0, 0, 0, 1], 'Y': [0, 0.5, 0, 0.5, 0], 'K': [0, 0, 0.5, 0.5, 0], 'S': [0, 0.5, 0.5, 0, 0],  'M': [0.5, 0.5, 0, 0, 0], 'R': [0.5, 0, 0.5, 0, 0], 'W': [0.5, 0, 0, 0.5, 0], 'V': [0.3333, 0.3333, 0.3333, 0, 0],  'H': [0.3333, 0.3333, 0, 0.3333, 0], 'D': [0.3333, 0, 0.3333, 0.3333, 0], 'B': [0, 0.3333, 0.3333, 0.3333, 0], 'N': [0.25, 0.25, 0.25, 0.25, 0]}
dna_df = pd.DataFrame(dna, index=['A','C','G','T','-'])
# print(dna_df)
file = sys.argv[1]
count_df = pd.DataFrame({0: [0, 0, 0, 0, 0]}, index=['A','C','G','T','-'])
for record in SeqIO.parse(file, "fasta"):
  seq = list(str(record.seq).upper())
  for i in range(len(seq)):
      if i in count_df.columns:
          count_df[i] += dna_df[seq[i]]
      else:
          count_df[i] = [0, 0, 0, 0, 0]
          count_df[i] += dna_df[seq[i]]
count_norm = count_df.div(count_df.sum(axis=0), axis=1)

ga_trint = pd.DataFrame({0: [0, 0, 0, 0, 0], 1: [0, 0, 0, 0, 0], 2: [0, 0, 0, 0, 0]}, index=['A','C','G','T','-'])
ct_trint = pd.DataFrame({0: [0, 0, 0, 0, 0], 1: [0, 0, 0, 0, 0], 2: [0, 0, 0, 0, 0]}, index=['A','C','G','T','-'])
gc_trint = pd.DataFrame({0: [0, 0, 0, 0, 0], 1: [0, 0, 0, 0, 0], 2: [0, 0, 0, 0, 0]}, index=['A','C','G','T','-'])
ct_count = 0
gc_count = 0
variable_sites = 0
for (columnName, columnData) in count_norm.iteritems() :
    if columnData['-'] == 0 and columnData['A'] != 1 and columnData['C'] != 1 and columnData['G'] != 1 and columnData['T'] != 1 and columnName < len(seq):
        variable_sites += 1
        if columnName+1 in count_df.columns:
            if columnData['A'] > 0 and columnData['G'] > 0 and columnData['C'] == 0 and columnData['T'] == 0:
                ga_trint[0] += count_df[columnName+1]
                ga_trint[2] += count_df[columnName-1]
                ct_count += 1
                print(str(columnName) + "\tCT\t" + file)
            if columnData['T'] > 0 and columnData['C'] > 0 and columnData['A'] == 0 and columnData['G'] == 0:
                ct_trint[0] += count_df[columnName-1]
                ct_trint[2] += count_df[columnName+1]
                ct_count += 1
                print(str(columnName) + "\tCT\t" + file)
            if columnData['G'] > 0 and columnData['C'] > 0 and columnData['A'] == 0 and columnData['T'] == 0:
                gc_trint[0] += count_df[columnName-1]
                gc_trint[2] += count_df[columnName+1]
                gc_count += 1
ct_logo = lm.Logo(ct_trint.transpose(), color_scheme='classic', baseline_width=0, font_name='Arial', show_spines=False, vsep=.005, width=.95)
ga_logo = lm.Logo(ga_trint.transpose(), color_scheme='classic', baseline_width=0, font_name='Arial', show_spines=False, vsep=.005, width=.95)
gc_logo = lm.Logo(gc_trint.transpose(), color_scheme='classic', baseline_width=0, font_name='Arial', show_spines=False, vsep=.005, width=.95)

ct_logo.fig.savefig(file + '.ct.pdf')
ga_logo.fig.savefig(file + '.ga.pdf')
gc_logo.fig.savefig(file + '.gc.pdf')

outFile = open(file + '.counts.txt', 'w')
outFile.write("\t".join([str(ct_count), str(gc_count), str(variable_sites)]))

# base_motif = motifs.create(seqs, IUPAC.ambiguous_dna)
# pwm = pd.DataFrame(base_motif.counts)
# print(pwm)
# outContent = str(pwm) + "\n"
# print(outContent)
# outFile = file + ".pwm.txt"
# outFH = open(outFile, "w")
# outFH.write(outContent)
