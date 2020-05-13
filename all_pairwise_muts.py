#!/usr/bin/env python3
import argparse
from Bio import AlignIO
from Bio import Phylo
from Bio.Seq import Seq
import scipy.stats as stats

### variant caller
def pairwise_variant_caller(query,count,a_aln,q_aln,path):
  subject = path[count].name
  refBases = list(a_aln[subject].upper())
  recordBases = list(q_aln[query].upper())
  depth = 0.0
  for j in range(count, len(path)) :
      depth += path[j].branch_length
  i = 1
  max = len(refBases)-1
  #print max
  TCK = 0
  VCK = 0
  TDN = 0
  VDN = 0
  variants = []
  while i < max:
      if refBases[i] != recordBases[i] and refBases[i] != "-" and recordBases[i] != "-" and refBases[i] != "N" and recordBases[i] != "N":
          strand = 1
          dint5 = recordBases[i-1] + refBases[i]
          dint3 = refBases[i] + recordBases[i+1]
          if refBases[i] == "G" or refBases[i] == "A":
              dint5 = str(Seq(dint3).reverse_complement())
              trint = str(Seq(refBases[i+1]).reverse_complement()) + "x" + str(Seq(refBases[i-1]).reverse_complement())
              mut = str(Seq(refBases[i]).reverse_complement()) + ">" + str(Seq(recordBases[i]).reverse_complement())
              strand = -1
          else:
              trint = refBases[i-1] + "x" + refBases[i+1]
              mut = refBases[i] + ">" + recordBases[i]
          if (mut == "C>T" or mut == "C>G" or mut == "C>Y" or mut == "C>K" or mut == "C>S") and trint.startswith('T') and '-' not in trint:
              TCK +=1
          elif (mut == "C>T" or mut == "C>G" or mut == "C>Y" or mut == "C>K" or mut == "C>S") and not trint.startswith('T') and '-' not in trint:
              VCK +=1
          elif (mut != "C>T" and mut != "C>G" or mut != "C>Y" or mut != "C>K" or mut != "C>S") and trint.startswith('T') and '-' not in trint:
              TDN +=1
          elif (mut != "C>T" and mut != "C>G" or mut != "C>Y" or mut != "C>K" or mut != "C>S") and not trint.startswith('T') and '-' not in trint:
              VDN +=1
          pos_res = "\t".join([query, subject, str(i), dint5, trint, mut, str(strand)])
          variants.append(pos_res)
      i += 1
  oddsratio, pvalue = stats.fisher_exact([[TCK, TDN], [VCK, VDN]])
  summary = [query, subject, str(depth), str(TCK), str(VCK), str(TDN), str(VDN), str(pvalue), str(oddsratio), "\n"]
  return(variants,summary)

### Main

parser = argparse.ArgumentParser()
parser.add_argument('tips', type=str, help='Input tips of interest (seqnames)')
parser.add_argument('aln', type=str, help='Alignment in fasta format')
parser.add_argument('tree', type=str, help='Ancestral tree')
parser.add_argument('states', type=str, help='Ancestral states fasta')
args = parser.parse_args()

tree = Phylo.read(args.tree, 'newick')

alignment = AlignIO.read(args.aln, "fasta")
aln_dict = dict()
for record in alignment :
    aln_dict[record.id] = str(record.seq)

states = AlignIO.read(args.states, "fasta")
states_dict = dict()
for record in states :
    states_dict[record.id] = str(record.seq)

variant_out = open(args.states + ".variants.txt", "w")
summary_out = open(args.states + ".summary.txt", "w")
split_tips = args.tips.split(",")

for tip in split_tips:
  # print(tip)
  node_path = tree.get_path(tip)
  node_path.pop(0)
  path_depth = len(node_path)
  for n in range(path_depth-1):
    var, res = pairwise_variant_caller(tip, n, states_dict, aln_dict, node_path)
    summary_out.write("\t".join(res))
    variant_out.write("\n".join(var))

# get tip distance to each parent node ## remember to not compare root node
