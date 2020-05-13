#/usr/bin/env python3

import argparse
from Bio import SeqIO
import re
import sys

parser = argparse.ArgumentParser()
parser.add_argument('fasta', metavar='N', type=str, help='Input fasta file')
parser.add_argument('clstr', metavar='N', type=str, help='Input cdhit clstr file')
parser.add_argument('out', metavar='N', type=str, help='Output file name')
args = parser.parse_args()

record_dict = dict()
with open(args.fasta, "r") as input_handle:
    for record in SeqIO.parse(input_handle, "fasta"):
        name = record.name.split("::")[0]
        # print(name)
        record_dict[name] = record

out_dict = dict()
counter = 0

with open(args.clstr) as fp:
    Lines = fp.readlines()
    for line in Lines:
      if re.match(">Cluster \d+", line):
        counter = re.match(">Cluster (\d+)", line).group(1)
        #print(counter)
        out_dict[counter] = dict()
      elif re.match(".+>(.+)(::)?\.\.\..+", line):
        seq_id = re.match(".+>(.+)(::)?\.\.\..+", line).group(1)
        print("{}\t{}".format(counter,seq_id))
        out_dict[counter][seq_id] = record_dict[seq_id]

for count in out_dict:
  out_file = "".join([args.out,"_",count,".fa"])
  with open(out_file, 'w') as handle:
    for seq_id in out_dict[count]:
      SeqIO.write(out_dict[count][seq_id], handle, "fasta")
