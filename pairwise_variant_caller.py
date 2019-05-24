#!/usr/bin/env python

#from Bio import Align
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import numpy as np
#aligner = Align.PairwiseAligner()

inFasta = sys.argv[2]
refFasta = sys.argv[1]

refSeq = SeqIO.read(refFasta, "fasta")
refSeq.seq = Seq(str(refSeq.seq).replace("-",""))

for record in SeqIO.parse(inFasta, "fasta"):
    record.seq = Seq(str(record.seq).replace("-",""))
    alignments = pairwise2.align.globalms(refSeq.seq, record.seq, 2, -1, -3, -.1)
    #print(alignments[0][0])
    refBases = list(alignments[0][0])
    recordBases = list(alignments[0][1])
    i = 1
    max = len(refBases)-1
    #print max
    while i < max:
        if refBases[i] != recordBases[i] and refBases[i] != "-" and recordBases[i] != "-" and refBases[i] != "N" and recordBases[i] != "N":
            strand = 1
            dint5 = recordBases[i-1] + refBases[i]
            dint3 = refBases[i] + recordBases[i+1]
            if refBases[i] == "G" or refBases[i] == "A":
                dint5 = str(Seq(dint3).reverse_complement())
                trint = str(Seq(recordBases[i+1]).reverse_complement()) + "x" + str(Seq(recordBases[i-1]).reverse_complement())
                mut = str(Seq(refBases[i]).reverse_complement()) + ">" + str(Seq(recordBases[i]).reverse_complement())
                strand = -1
            else:
                trint = recordBases[i-1] + "x" + recordBases[i+1]
                mut = refBases[i] + ">" + recordBases[i]
            print "\t".join([record.id, refSeq.id, str(i), dint5, trint, mut, str(strand)])
        i += 1
