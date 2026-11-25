source('scripts/r/functions.R')

args <- commandArgs(trailingOnly = T)
filtered_og_joined <- args[1]
OG <- args[2]
out_path <- args[3]
#filtered_og_joined <- "data/test_mmseqs_cog2020/fast.filtered.cog2020.tsv"
#OG <- "cog"
#out_path <- "data/mmseqs_tsv_filtered_cog2020/"

out_file <- sub(".tsv", ".assigned.tsv", sub(".*/", "", filtered_og_joined))
df_filtered_og_joined <- read.csv(filtered_og_joined, header = TRUE, sep = "\t", quote = "") %>%
	as_tibble() %>%
	dplyr::rename(is_multi_hits_h = "is_multi_domain_2")

if (OG == "cog") {
	df_filtered_og_assigned <- assign_cog() %>% left_join(read_cog2020_def_fixed(), by = "cog_id")
} else if (OG == "kegg") {
	df_filtered_og_assigned <- assign_ko() %>% left_join(read_ko_def(), by = "ko_id")
}

write.table(df_filtered_og_assigned, paste0(out_path, out_file), append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE)
