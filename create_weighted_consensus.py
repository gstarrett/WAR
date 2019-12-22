#!/usr/bin/env python

import math
from Bio import SeqIO
from Bio import motifs
import sys

###
file = sys.argv[1]
prefix = sys.argv[2]
outFile = prefix + ".consensus.fasta"
outFH = open(outFile, "w")
seqs = []
for record in SeqIO.parse(file, "fasta"):
  seqs.append(record.seq)

base_motif = motifs.create(seqs)
pwm = base_motif.counts.normalize(pseudocounts={'A':0.307, 'C': 0.186, 'G': 0.212, 'T': 0.294})
outContent = ">" + prefix + "\n" + str(pwm.consensus) + "\n"
outFH.write(outContent)
