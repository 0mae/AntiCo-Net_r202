source('scripts/r/functions.R')

args <- commandArgs(trailingOnly = T)

OG <- args[1]
each_tax <- args[2]
out_path <- args[3]
#OG <- "cog"
#each_tax <- "p__Proteobacteria"

if (OG == "cog"){

	dt_genomes_ogs_tidy <- read.csv(
		paste0("data/proteins_raw_gtdbtax_cog2020_tidy/proteins_raw_gtdbtax_cog2020_tidy_", each_tax, ".tsv"),
		header = TRUE, sep = "\t", quote = ""
	) %>%
	dplyr::rename(og_id = "cog_id") %>%
	mutate(presence = count > 0) %>%
	data.table() %>%
	lazy_dt()

} else if (OG == "kegg") {

	dt_genomes_ogs_tidy <- read.csv(
		paste0("data/proteins_raw_gtdbtax_ko_tidy/proteins_raw_gtdbtax_ko_tidy_", each_tax, ".tsv"),
		header = TRUE, sep = "\t", quote = ""
	) %>%
	dplyr::rename(og_id = "ko_id") %>%
	mutate(presence = count > 0) %>%
	data.table() %>%
	lazy_dt()

}

dt_pivot <- dt_genomes_ogs_tidy %>% select(-tax,-count) %>%
	pivot_wider(names_from = "og_id", values_from = presence, values_fill = FALSE)

mx_genome_by_og <- dt_pivot %>% select(-accession) %>% as.data.frame() %>% as.matrix()
rownames(mx_genome_by_og) <- as_tibble(dt_pivot)$accession

mx_og_by_genome <- t(mx_genome_by_og)

mx_all_true <- matrix(TRUE, nrow = dim(mx_genome_by_og)[1], ncol = dim(mx_genome_by_og)[2])
rownames(mx_all_true) <- rownames(mx_genome_by_og)
colnames(mx_all_true) <- colnames(mx_genome_by_og)

mx_genome_by_og_rev <- (mx_all_true - mx_genome_by_og) > 0
mx_og_by_genome_rev <- t(mx_genome_by_og_rev)

mx_TT_present <- mx_og_by_genome %*% mx_genome_by_og 
mx_TT <- mx_TT_present
mx_TT[lower.tri(mx_TT, diag = TRUE)] <- NA

og_present <- diag(mx_TT_present)

mx_TF <- mx_og_by_genome %*% mx_genome_by_og_rev
mx_TF[lower.tri(mx_TF, diag = TRUE)] <- NA

mx_FT <- mx_og_by_genome_rev %*% mx_genome_by_og
mx_FT[lower.tri(mx_FT, diag = TRUE)] <- NA

mx_FF <- mx_og_by_genome_rev %*% mx_genome_by_og_rev
mx_FF[lower.tri(mx_FF, diag = TRUE)] <- NA

df_TT <- mx_TT %>%
	as.data.frame() %>%
	rownames_to_column(var = "og_id_A") %>%
	as_tibble() %>%
	pivot_longer(
		cols = -og_id_A,
		names_to = "og_id_B",
		values_to = "count",
		values_drop_na = TRUE
	) %>%
	mutate(type = "TT")

df_TF <- mx_TF %>%
	as.data.frame() %>%
	rownames_to_column(var = "og_id_A") %>%
	as_tibble() %>%
	pivot_longer(
		cols = -og_id_A,
		names_to = "og_id_B",
		values_to = "count",
		values_drop_na = TRUE
	) %>%
	mutate(type = "TF")

df_FT <- mx_FT %>%
	as.data.frame() %>%
	rownames_to_column(var = "og_id_A") %>%
	as_tibble() %>%
	pivot_longer(
		cols = -og_id_A,
		names_to = "og_id_B",
		values_to = "count",
		values_drop_na = TRUE
	) %>%
	mutate(type = "FT")

df_FF <- mx_FF %>%
	as.data.frame() %>%
	rownames_to_column(var = "og_id_A") %>%
	as_tibble() %>%
	pivot_longer(
		cols = -og_id_A,
		names_to = "og_id_B",
		values_to = "count",
		values_drop_na = TRUE
	) %>%
	mutate(type = "FF")

df_contingency <- rbind(df_TT, df_TF, df_FT, df_FF) %>%
	pivot_wider(names_from = "type", values_from = "count", values_fill = 0) %>%
	mutate(n_genomes = TT + TF + FT + FF) %>%
	left_join(
		tibble(og_id_A = names(og_present), n_present_A = og_present),
		by = "og_id_A"
	) %>%
	left_join(
		tibble(og_id_B = names(og_present), n_present_B = og_present),
		by = "og_id_B"
	) %>%
	mutate(
		contrapositive_ratio = (TF + FT) / n_genomes,
		cooccurring_ratio = (TT + FF) / n_genomes,
		ratio_present_A = n_present_A / n_genomes,
		ratio_present_B = n_present_B / n_genomes,
		tax = each_tax
	)
	
write.table(df_contingency, out_path, append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE)
