source('scripts/r/functions.R')

args <- commandArgs(trailingOnly = T)
out_path <- args[1]
n <- as.numeric(args[2])
out_path2 <- sub(".tsv$", paste0("_", as.character(n), ".tsv"), out_path)

df_genomes_analyzed <- read_genomes_analyzed()

df_summarized_Domain <- summarize_taxa("Domain")
df_summarized_Phylum <- summarize_taxa("Phylum")
df_summarized_Class <- summarize_taxa("Class")
df_summarized_Order <- summarize_taxa("Order")
df_summarized_Family <- summarize_taxa("Family")
df_summarized_Genus <- summarize_taxa("Genus")
df_summarized_Species <- summarize_taxa("Species")

df_summarized_taxa <- dplyr::bind_rows(
        df_summarized_Domain,
        df_summarized_Phylum,
        df_summarized_Class,
        df_summarized_Order,
        df_summarized_Family,
        df_summarized_Genus,
        df_summarized_Species
    ) %>%
    mutate(tax_path = sub("$", ";", sub(";NA.*", "", paste(Domain,Phylum,Class,Order,Family,Genus,Species, sep = ";")))) %>%
    arrange(tax_path)

write.table(df_summarized_taxa, out_path, na = "", append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE)

write.table(
    df_summarized_taxa %>% filter(n_genomes >= n),
    out_path2, na = "", append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE
    )
