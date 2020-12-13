
# refdb: Reference database manager for metabarcoding

## Import

- NCBI (rentrez)
- BOLD (bold)
- local

`read_refdb(file)` reads native refdb file (compressed archive with a csv and a yaml file)

`refdb_import_NCBI()`
`refdb_import_BOLD()`
`refdb_import_local()`

`refdb_set_fields(id, taxonomy, sequence, gene, coordinates, yaml_file = NULL)`
`refdb_merge(..., taxo_ncbi = FALSE)`


## Clean

- Sequences
    - Clean illegal character (cf. bioseq)
    - Clean gaps
    - Clean multi-N on sides
    - Crop barcode region using primers

`refdb_clean_sequences()`

- Taxonomy
    - Typo detection using pairwise string distance
    - Harmonization (taxize)
    - Delete replicated spaces, tabs and new lines
    - Optionally delete "cf." and "aff." expressions
    - Optionally delete infra specific taxonomic details (varieties, form, morph, sippe, etc.)
    - Delete possible remaining spaces at the beginning/end of the taxonomic name
    - Replace "var" with "var."
    - Delete words containing more than one uppercase and/or numbers, words of one character
    - Force lower case for all characters and upper case for the first character

`refdb_clean_taxonomy()`


## Filter

- Sequence length
- Number of degenerated bases
- Primer detection (forward and reverse)
- Homopolymers max
- Pairwise nucleotide distance to consensus (based on a subsample)
- Adequation phylogeny and taxonomy
- Identical barcodes corresponding to identical species
- Unidentified species (sp. and one word species)
- Reference database self-assignation
- Amino-acid conversion and codon stop detection
-Same species (names) with different taxonomy

`refdb_filter()`

Rating of sequences using a quality score.

Multicriteria general score (Shiny interface)

`refdb_filter_gui()`


## Export
- Export to csv/fasta (compatible mothur, dada2, quime, ...)
- Report and dashboard generator
- Website generator
- Github Actions

`write_refdb(file)`

`refdb_export_csv()`
`refdb_export_dada2()`
`refdb_export_mothur()`

`refdb_dahboard()`
`refdb_website()`

## Plots
- General statistics
- Taxonomy coverage (tree, treemaps)
- Maps

`refdb_plot()`
