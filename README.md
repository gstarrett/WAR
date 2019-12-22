# WAR
## Weighted Ancestor Reconstruction

In recent years it has been observed that many viral genomes have unbalanced nucleotide composition likely due to effects of error prone polymerases, antiviral deaminases, and other selective pressures. With sufficient sampling of a population and evolutionary distance from extant members, ancestor state reconstruction approaches likely will create a genome with nucleotide composition that reflects the mutational processes acting on viral genomes. However, if sampling is limited and/or if the distance from the extant species to the desired ancestor is short, there is a likely chance that the ancestor state will not accurately reflect the expected genome composition at the nucleotide level, especially when there is a known strong nucleotide bias, such as has been observed in alphapaplillomaviruses and polyomaviruses.

WAR was developed to study intra-type point mutations differences in alphapapillomaviruses. Alphapapillomaviruses induce the mutagenic antiviral enzymes APOBEC3A and APOBEC3B, which are important for somatic mutations in many cancer types. These enzymes also act on the viral genome and the di- and trinucleotide composition of these genomes reflects this evolutionary conflict. Acute APOBEC3 mutagenesis have been observed through deep sequencing of infected tissues and cervical intraepithelial neoplasia. In the absence of these data for more exotic carcinogenic types, it is more challenging 

### Instructions
#### Get base frequencies for viral genomes
* Input all viral genomes as a multi fasta file to calculate the average single base frequency for a particular virus (or type, species, etc) of interest
  * `get_base_freqs.pl <input.fasta>`

#### Generate ancestral state reconstructions
* Align viral isolates using MAFFT
* Remove identical duplicated sequences
* Build a Maximum Likelihood tree of viral isolates us MEGA X from MAFFT alignment
* Input tree or selected subtree into FastML
* Select an ancestral node sufficiently far back for ancestral state reconstruction (this will vary depending on the input data and sequence of interest)
* Make 100 ancestral reconstructions of this node

#### Generate weighted consensus ancestor
* Input ancestral reconstructions in fasta format and previously generated base frequencies into `create_weighted_consensus.py`
  * `create_weighted_consensus.py <base_freqs_output> <ancestral_reconstruction.fasta> <output_prefix>`

#### Call variants between extant and ancestor genomes
* Requires a 2 single sequence fasta inputs (ancestor reference and sequence of interest)
  * `pairwise_variant_caller.py <ancestor.fasta> <extant.fasta>`

### Dependencies
* Bioperl (get_base_freqs.pl will be rewritten to use biopython eventually to eliminate this dependency)
* Biopython
