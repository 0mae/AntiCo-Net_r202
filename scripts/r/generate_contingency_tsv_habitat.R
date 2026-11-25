source('scripts/r/functions.R')

args <- commandArgs(trailingOnly = T)

OG <- args[1]
each_tax <- args[2]
out_path <- args[3]
hp_thre <- as.numeric(args[4])
#OG <- "cog"
#each_tax <- "p__Proteobacteria"
#hp_thre <- 5

if (OG == "cog" & each_tax == "k__Prokaryotes"){

	dt_genomes_ogs_tidy <- fread(
		paste0("data/tsv/proteins_raw_gtdbtax_cog2020_tidy.tsv"),
		header = TRUE, sep = "\t", quote = ""
	) %>%
	dplyr::rename(og_id = "cog_id", tax = "Domain") %>%
	mutate(presence = count > 0) %>%
	select(-Phylum:-Species) %>%
	lazy_dt()

} else if (OG == "cog"){

	dt_genomes_ogs_tidy <- read.csv(
		paste0("data/proteins_raw_gtdbtax_cog2020_tidy/proteins_raw_gtdbtax_cog2020_tidy_", each_tax, ".tsv"),
		header = TRUE, sep = "\t", quote = ""
	) %>%
	dplyr::rename(og_id = "cog_id") %>%
	mutate(presence = count > 0) %>%
	data.table() %>%
	lazy_dt()

} else if (OG == "kegg" & each_tax == "k__Prokaryotes"){

	dt_genomes_ogs_tidy <- fread(
		paste0("data/tsv/proteins_raw_gtdbtax_ko_tidy.tsv"),
		header = TRUE, sep = "\t", quote = ""
	) %>%
	dplyr::rename(og_id = "ko_id", tax = "Domain") %>%
	mutate(presence = count > 0) %>%
	select(-Phylum:-Species) %>%
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

df_habitat_pref <-
	read.csv("data/tsv/gtdb_ssu_prokatlas_habitat_pref.tsv", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
	as_tibble()

intersect_genomes <- intersect(as_tibble(dt_pivot)$accession, df_habitat_pref$accession)

mx_genome_by_og <- dt_pivot %>% select(-accession) %>% as.data.frame() %>% as.matrix()
rownames(mx_genome_by_og) <- as_tibble(dt_pivot)$accession
mx_genome_by_og <- mx_genome_by_og[intersect_genomes,]

mx_og_by_genome <- t(mx_genome_by_og)

mx_all_true <- matrix(TRUE, nrow = dim(mx_genome_by_og)[1], ncol = dim(mx_genome_by_og)[2])
rownames(mx_all_true) <- rownames(mx_genome_by_og)
colnames(mx_all_true) <- colnames(mx_genome_by_og)

mx_genome_by_og_rev <- (mx_all_true - mx_genome_by_og) > 0
mx_og_by_genome_rev <- t(mx_genome_by_og_rev)

mx_habitat_pref <- df_habitat_pref %>%
	select(-accession) %>%
	as.matrix()
rownames(mx_habitat_pref) <- df_habitat_pref$accession

mx_genome_by_habitat <- mx_habitat_pref[intersect_genomes,] > hp_thre
mx_habitat_by_genome <- t(mx_genome_by_habitat)

mx_all_true_habitat <- matrix(TRUE, nrow = dim(mx_genome_by_habitat)[1], ncol = dim(mx_genome_by_habitat)[2])
rownames(mx_all_true_habitat) <- rownames(mx_genome_by_habitat)
colnames(mx_all_true_habitat) <- colnames(mx_genome_by_habitat)

mx_genome_by_habitat_rev <- (mx_all_true_habitat - mx_genome_by_habitat) > 0
mx_habitat_by_genome_rev <- t(mx_genome_by_habitat_rev)

mx_TT_present <- mx_habitat_by_genome %*% mx_genome_by_og
mx_TT <- mx_TT_present

mx_TF <- mx_habitat_by_genome %*% mx_genome_by_og_rev

mx_FT <- mx_habitat_by_genome_rev %*% mx_genome_by_og

mx_FF <- mx_habitat_by_genome_rev %*% mx_genome_by_og_rev

df_TT <- mx_TT %>%
	as.data.frame() %>%
	rownames_to_column(var = "habitat") %>%
	as_tibble() %>%
	pivot_longer(
		cols = -habitat,
		names_to = "og_id",
		values_to = "count",
		values_drop_na = TRUE
	) %>%
	mutate(type = "TT")

df_TF <- mx_TF %>%
	as.data.frame() %>%
	rownames_to_column(var = "habitat") %>%
	as_tibble() %>%
	pivot_longer(
		cols = -habitat,
		names_to = "og_id",
		values_to = "count",
		values_drop_na = TRUE
	) %>%
	mutate(type = "TF")

df_FT <- mx_FT %>%
	as.data.frame() %>%
	rownames_to_column(var = "habitat") %>%
	as_tibble() %>%
	pivot_longer(
		cols = -habitat,
		names_to = "og_id",
		values_to = "count",
		values_drop_na = TRUE
	) %>%
	mutate(type = "FT")

df_FF <- mx_FF %>%
	as.data.frame() %>%
	rownames_to_column(var = "habitat") %>%
	as_tibble() %>%
	pivot_longer(
		cols = -habitat,
		names_to = "og_id",
		values_to = "count",
		values_drop_na = TRUE
	) %>%
	mutate(type = "FF")

df_contingency <- rbind(df_TT, df_TF, df_FT, df_FF) %>%
	pivot_wider(names_from = "type", values_from = "count", values_fill = 0) %>%
	mutate(n_genomes_hp = TT + TF + FT + FF) %>%
	mutate(
		contrapositive_ratio_hp = (TF + FT) / n_genomes_hp,
		cooccurring_ratio_hp = (TT + FF) / n_genomes_hp,
		tax = each_tax,
		tt = TT,
		tf = TF,
		ft = FT,
		ff = FF
	) %>%
	nest(
		contingency = c(tt,tf,ft,ff),
	) %>%
	mutate(
		p_value = map(contingency, ~ fisher.test(matrix(unlist(.x), nrow = 2, byrow=T))$p.value)
	) %>%
	select(-contingency) %>%
	unnest(c(p_value)) %>%
	mutate(
		q_value = p.adjust(p_value, method = "BH")
	)

write.table(df_contingency, out_path, append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE)
