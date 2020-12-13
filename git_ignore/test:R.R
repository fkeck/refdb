
# Test import NCBI

query <- "Ephemerellidae COI"
query <- "Fragilaria capucina var. mesolepta"

qn <- refdb_import_NCBI(query)

# Test import bold

qb <- refdb_import_BOLD("Baetidae", geo = "Switzerland")

