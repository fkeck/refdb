


refdb_filter_seq_length(x, min_len, max_len)
refdb_filter_seq_ambiguous(x, max = 3, char = "N")
refdb_filter_seq_primer()
refdb_filter_seq_homopolymers(x, max_len = 8)
refdb_filter_seq_duplicates(x)
refdb_filter_seq_subsample(x, tax, method, keep = 1)
refdb_filter_seq_stopcodon(x, frame)
refdb_filter_tax_precision(x, min_tax)
refdb_filter_study_precision(x, min_tax)

#Require alignment:
refdb_filter_seq_dist(x, max_dist)

#Require alignment or phylogenetic tree
refdb_filter_seq_phylo(x, phylo, max_dist)

#Require external RDP software
refdb_filter_seq_selfassign(x, max_dist, exec)
