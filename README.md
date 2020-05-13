# WAR
## Weighted Ancestor Reconstruction

In recent years it has been observed that many viral genomes have unbalanced nucleotide composition likely due to effects of error prone polymerases, antiviral deaminases, and other selective pressures. With sufficient sampling of a population and evolutionary distance from extant members, ancestor state reconstruction approaches likely will create a genome with nucleotide composition that reflects the mutational processes acting on viral genomes. However, if sampling is limited and/or if the distance from the extant species to the desired ancestor is short, there is a likely chance that the ancestor state will not accurately reflect the expected genome composition at the nucleotide level, especially when there is a known strong nucleotide bias, such as has been observed in alphapaplillomaviruses and polyomaviruses.

WAR was developed to study intra-type point mutations differences in alphapapillomaviruses. Alphapapillomaviruses induce the mutagenic antiviral enzymes APOBEC3A and APOBEC3B, which are important for somatic mutations in many cancer types. These enzymes also act on the viral genome and the di- and trinucleotide composition of these genomes reflects this evolutionary conflict. Acute APOBEC3 mutagenesis have been observed through deep sequencing of infected tissues and cervical intraepithelial neoplasia. In the absence of these data for more exotic carcinogenic types, it is more challenging

### Instructions
#### Get base frequencies for viral genomes (optional for weighting ancestral reconstruction)
* Input all viral genomes as a multi fasta file to calculate the average single base frequency for a particular virus (or type, species, etc) of interest
  * `get_base_freqs.pl <input.fasta>`

#### Generate ancestral state reconstructions
* Remove identical duplicated sequences
* Cluster sequences by cutoff using cd-hit-est (example HPV types >=90% nt similarity) and split using `split_cdhit_clusters.py <input.fasta> <cdhit.clstr> <out_prefix>`
* Align sequences using MAFFT
* Build a Maximum Likelihood tree of viral isolates using RAxML-NG from MAFFT alignment
* Input tree into RAxML-NG ancestor to make an ancestral tree, ancestral base probabilities, and predicted ancestral states
* Weight ancestral probabilities to counteract unbalanced mutagenic processes (optional)
* Input ancestral sequeces

#### Generate weighted consensus ancestor (optional)
* Input ancestral reconstruction probabilities and previously generated base frequencies into `create_weighted_seq.py` to inversely weight against base imbalances
  * `create_weighted_seq.py <ancestral_probabilities>`

#### Call variants between extant and ancestor genomes
* Requires a file of line separated sequence names, alignment in fasta format, tree with ancestral nodes, and an ancestral states fasta file to call variants versus ancestor sequences. This will also create a summary file and calculate significant APOBEC3 T[C>K] mutagenesis per sequence relative to ancestor
  * `all_pairwise_muts.py <tips.txt> <alignment.fasta> <tree.newick> <ancestors.fasta>`

#### Get pwm of flanking bases for all C/T variable sites and plot logo
* Input alignment in fasta format into `get_pwm.py` to get flanking bases for all sites that are only variable for C:G and T:A to determine nucleotide context and assess mutagenic processes
  * `get_pwm.py <alignment.fa>`

### Dependencies
* Bioperl (get_base_freqs.pl will be rewritten to use biopython eventually to eliminate this dependency)
* Biopython
* pandas
* Logomaker
