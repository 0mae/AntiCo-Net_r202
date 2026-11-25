source('scripts/r/functions.R')

# Load GTDB genome metadata table with slight modification (GTDB tax path is splited into separate columns)
gtdb_metadata <- load_gtdb_meta()
df_bac_tax <- gtdb_metadata[[1]]
df_ar_tax <- gtdb_metadata[[2]]

# Retrieve representative high-quality (completeness >= 95%, contamination <= 5%) genomes
hq_genomes <- rep_hq_genomes(df_bac_tax, df_ar_tax)

df_bac_hq <- hq_genomes[[1]]
df_ar_hq <- hq_genomes[[2]]

# Retrieve CPR representative genomes
df_cpr <- cpr_rep_genomes(df_bac_tax)

# Load CheckM result using CPR marker set
df_cpr_cm <- cpr_checkm()

# Retrieve CPR representative high-quality (based on CPR marker set) genomes
df_cpr_hq <- cpr_hq_rep_genomes(df_cpr_cm)

# Merge CPR representative high-quality (based on CPR marker set) genomes and others
hq_merged_genomes <- hq_genomes_merged(df_bac_hq, df_ar_hq, df_cpr_hq)
df_bac_hq_merged <- hq_merged_genomes[[1]]
df_ar_hq_merged <- hq_merged_genomes[[2]]

# Retrieve genomes of only taxa with >= 50 genomes, and generate final genome data set
#genomes_analyzed <- taxa_50_genomes(df_bac_hq_merged, df_ar_hq_merged)
#df_bac_analyzed <- genomes_analyzed[[3]]
#df_ar_analyzed <- genomes_analyzed[[4]]

# Output
write.table(df_bac_hq_merged, "data/tsv/genomes_analyzed_bac.tsv", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(df_ar_hq_merged, "data/tsv/genomes_analyzed_ar.tsv", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
#write.table(df_bac_analyzed, "data/genomes_analyzed_bac.tsv", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
#write.table(df_ar_analyzed, "data/genomes_analyzed_ar.tsv", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)