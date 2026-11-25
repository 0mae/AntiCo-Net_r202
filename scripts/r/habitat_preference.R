source('scripts/r/functions.R')

tbl <- fread("data/tsv/gtdb_ssu_prokatlas.tsv", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
	as_tibble() %>%
	rename(accession = "qseqid", biome = "sseqid")

# Calculate habitat preference

tbl_pivot <- tbl %>%
	group_by(accession,biome) %>%
	summarize(count = dplyr::n()) %>%
	ungroup() %>%
	pivot_wider(names_from = "biome", values_from = "count", values_fill = 0)

rowsum <- rowSums(select(tbl_pivot, -accession))

df_genome_nom <- cbind(select(tbl_pivot, accession), tbl_pivot %>% select(-accession) %>%
    mutate_all(~(.x / rowsum))) %>% as_tibble()

prokatlas_weight <- read.table('db/prokatlas/weight_correction.txt', header = TRUE,  sep = '\t', quote = "") %>%
    as_tibble()

prokatlas_weight$environment <- factor(prokatlas_weight$environment, levels = colnames(select(df_genome_nom, -accession)), ordered = TRUE)

prokatlas_weight_sort <- prokatlas_weight %>% arrange(environment) %>% rename(biome = "environment")

prokatlas_weight_wide <- prokatlas_weight_sort[prokatlas_weight_sort$biome %in% colnames(select(df_genome_nom, -accession)), ] %>% select(-total) %>% pivot_wider(names_from = "biome", values_from = "weight")

prokatlas_weight_expand <- tibble(key = rep("a", length = dim(df_genome_nom)[1])) %>%
    left_join(mutate(prokatlas_weight_wide, key = "a"), by = "key") %>% select(-key)

df_genome_weighted <- cbind(select(df_genome_nom, accession), select(df_genome_nom, -accession) * prokatlas_weight_expand) %>% as_tibble()

rowsum_genome_weighted <- rowSums(select(df_genome_weighted, -accession))

df_habitat_preference <- cbind(select(df_genome_weighted, accession), df_genome_weighted %>% select(-accession) %>% mutate_all(~(.x * 100 / rowsum_genome_weighted))) %>% as_tibble()

write.table(df_habitat_preference, file = "data/tsv/gtdb_ssu_prokatlas_habitat_pref.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

# Statistics of habitat preference

df_genomes_analyzed <- read_genomes_analyzed()

df_tax_color <- rbind(
		df_genomes_analyzed %>%
			filter(Domain == "d__Bacteria") %>%
			group_by(Phylum) %>%
			summarize(n_genomes = dplyr::n()) %>%
			arrange(desc(n_genomes)) %>%
			head(10) %>%
			mutate(tax_color = Phylum),
		df_genomes_analyzed %>%
			filter(Domain == "d__Archaea") %>%
			group_by(Phylum) %>%
			summarize(n_genomes = dplyr::n()) %>%
			mutate(tax_color = "d__Archaea")
	)

df_habitat_preference_taxa <- df_genomes_analyzed %>% select(accession,Domain:Species) %>%
	left_join(df_habitat_preference, by = "accession") %>%
	left_join(df_tax_color, by = "Phylum")

df_habitat_preference_taxa$tax_color[is.na(df_habitat_preference_taxa$tax_color)] <- "Others"

df_habitat_preference_taxa$tax_color <- factor(df_habitat_preference_taxa$tax_color, levels = c("d__Archaea", df_tax_color$tax_color[1:10], "Others"))

df_ngenomes_taxa <- df_genomes_analyzed %>% select(accession,Domain:Species) %>%
	mutate(domain = Domain) %>%
	pivot_longer(cols = Domain:Species, names_to = "tax_level", values_to = "tax") %>%
	group_by(tax_level, domain, tax) %>%
	summarize(n_genomes = dplyr::n()) %>%
	ungroup()

df_ngenomes_taxa$tax_level <- factor(df_ngenomes_taxa$tax_level, levels = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

p_n_genomes_habitat <- df_habitat_preference_taxa %>%
	pivot_longer(cols = Domain:Species, names_to = "tax_level", values_to = "tax") %>%
	filter(is.na(aquatic)) %>%
	group_by(tax) %>%
	summarize(n_genomes_habitat = dplyr::n()) %>%
	ungroup() %>%
	right_join(df_ngenomes_taxa, by = "tax") %>%
	mutate(ratio_genome_habitat = 1 - (n_genomes_habitat / n_genomes)) %>%
	ggplot(aes(x = tax_level, y = ratio_genome_habitat, color = n_genomes)) +
	geom_jitter() +
	geom_violin(fill = NA) +
	facet_grid(~ domain) +
	scale_color_viridis(trans = "log") +
	labs(title = "Ratio of habitat assigned-genomes", subplot = "Each plot denotes a taxon", x = "Taxon level", y = "Ratio", color = "Number of genomes\nin each taxon")
ggsave("png/prokatlas_n_genomes.png", p_n_genomes_habitat, height = 4, width = 10)

p_habitat_pref_dist <- df_habitat_preference_taxa %>%
	filter(!is.na(aquatic)) %>%
	pivot_longer(cols = c(-accession:-Species,-n_genomes,-tax_color), names_to = "environments", values_to = "habitat_preference") %>%
	mutate(
		hp_binned = cut(habitat_preference, seq(0, 100, by = 1))
	) %>%
	filter(!is.na(hp_binned)) %>%
	distinct(accession,hp_binned) %>%
	group_by(hp_binned) %>%
	summarize(n_genomes = dplyr::n()) %>%
	ggplot(aes(x = hp_binned, y = n_genomes)) +
	geom_bar(stat = "identity") +
	scale_fill_manual(values = cbp12) +
	theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))
ggsave("png/prokatlas_habitat_pref_dist.png", p_habitat_pref_dist, height = 5, width = 12)

p_habitat_pref_n_genomes <- df_habitat_preference_taxa %>%
	filter(!is.na(aquatic)) %>%
	pivot_longer(cols = c(-accession:-Species,-n_genomes,-tax_color), names_to = "environments", values_to = "habitat_preference") %>%
	mutate(
		hp_00 = habitat_preference > 0,
		hp_05 = habitat_preference > 5,
		hp_10 = habitat_preference > 10,
		hp_20 = habitat_preference > 20,
		hp_40 = habitat_preference > 40,
		hp_80 = habitat_preference > 80,
		hp_90 = habitat_preference > 90
	) %>%
	select(accession,tax_color,hp_00:hp_90) %>%
	pivot_longer(cols = c(-accession,-tax_color), names_to = "hp_binned", values_to = "presence") %>%
	filter(presence == TRUE) %>%
	distinct(accession,tax_color,hp_binned) %>%
	group_by(tax_color,hp_binned) %>%
	summarize(n_genomes = dplyr::n()) %>%
	ggplot(aes(y = n_genomes, x = hp_binned, fill = tax_color)) +
	geom_bar(stat = "identity") +
	scale_fill_manual(values = cbp12)
ggsave("png/prokatlas_habitat_pref_n_genomes.png", p_habitat_pref_n_genomes, height = 4, width = 7)

p_habitat_pref <- df_habitat_preference_taxa %>%
	filter(!is.na(aquatic)) %>%
	pivot_longer(cols = c(-accession:-Species,-n_genomes,-tax_color), names_to = "environments", values_to = "habitat_preference") %>%
	mutate(
		hp_00 = habitat_preference > 0,
		hp_10 = habitat_preference > 10,
		hp_20 = habitat_preference > 20,
		hp_40 = habitat_preference > 40,
		hp_90 = habitat_preference > 90
	) %>%
	pivot_longer(cols = hp_00:hp_90, names_to = "hp_binned", values_to = "presence") %>%
	filter(presence == TRUE) %>%
	group_by(environments, tax_color, hp_binned) %>%
	summarize(n_genomes = dplyr::n()) %>%
	mutate(percent_genomes = 100 * n_genomes / dim(df_genomes_analyzed)[1]) %>%
	ggplot(aes(x = percent_genomes, y = environments, fill = tax_color)) +
	geom_bar(stat = "identity") +
	facet_grid(~ hp_binned) +
	scale_fill_manual(values = cbp12)
ggsave("png/prokatlas_habitat_pref.png", p_habitat_pref, height = 12, width = 10)

p_habitat_pref_binary <- df_habitat_preference_taxa %>%
	filter(!is.na(aquatic)) %>%
	pivot_longer(cols = c(-accession:-Species,-n_genomes,-tax_color), names_to = "environments", values_to = "habitat_preference") %>%
	mutate(
		hp_00 = habitat_preference > 0,
		hp_10 = habitat_preference > 10,
		hp_20 = habitat_preference > 20,
		hp_40 = habitat_preference > 40,
		hp_90 = habitat_preference > 90
	) %>%
	pivot_longer(cols = hp_00:hp_90, names_to = "hp_binned", values_to = "presence") %>%
	filter(presence == TRUE) %>%
	group_by(environments, tax_color, hp_binned) %>%
	summarize(n_genomes = dplyr::n()) %>%
	mutate(percent_genomes = 100 * n_genomes / dim(df_genomes_analyzed)[1]) %>%
	ggplot(aes(size = percent_genomes, x = tax_color, y = environments, color = tax_color)) +
	geom_point() +
	facet_grid(~ hp_binned) +
	scale_color_manual(values = cbp12) +
	theme(axis.text.x = element_blank())
ggsave("png/prokatlas_habitat_pref_binary.png", p_habitat_pref_binary, height = 12, width = 10)





#p_habitat_pref <- df_habitat_preference_taxa %>%
#	filter(!is.na(aquatic)) %>%
#	pivot_longer(cols = c(-accession:-Species,-n_genomes,-tax_color), names_to = "environments", values_to = "habitat_preference") %>%
#	filter(habitat_preference > 0) %>%
#	group_by(environments, tax_color) %>%
#	summarize(n_genomes = dplyr::n()) %>%
#	mutate(percent_genomes = 100 * n_genomes / dim(df_genomes_analyzed)[1]) %>%
#	ggplot(aes(x = percent_genomes, y = environments, fill = tax_color)) +
#	geom_bar(stat = "identity") +
#	scale_fill_manual(values = cbp12)
#ggsave("png/prokatlas_habitat_pref.png", p_habitat_pref, height = 12, width = 7)





	#head() %>% print(width = Inf)