#library(tidyverse)
source('scripts/r/functions.R')

# Load GTDB genome metadata table with slight modification (GTDB tax path is splited into separate columns)
gtdb_metadata <- load_gtdb_meta()

df_bac_tax <- gtdb_metadata[[1]]

# Retrieve all CPR genomes
df_cpr <- cpr_rep_genomes(df_bac_tax)

write.table(df_cpr %>% select(accession), "data/gtdb/cpr_metadata.tsv", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)