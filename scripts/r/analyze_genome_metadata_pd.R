source('scripts/r/functions.R')

df_analyzed_genomes <- read_genomes_analyzed()

df_analyzed_genomes_stats <- df_analyzed_genomes %>%
	select(
		coding_bases,
		coding_density,
		contig_count,
		gc_count,
		gc_percentage,
		genome_size,
		protein_count,
		Domain,
		Phylum,
		Class,
		Order,
		Family,
		Genus
	) %>%
	mutate(Kingdom = "k__Prokaryotes") %>%
	pivot_longer(cols = c(-coding_bases:-protein_count), names_to = "tax_level", values_to = "tax") %>%
	select(-tax_level) %>%
	group_by(tax) %>%
	summarise(across(everything(), list(mean = mean, min = min, max = max, sd = sd)))

df_taxa_summarized <- read.csv("data/tsv/taxa_summarized.tsv", header = TRUE, sep = "\t", quote = "") %>% as_tibble()

gtdb_tree_ar <- read_gtdb_tree_ar()
gtdb_tree_bac <- read_gtdb_tree_bac()

df_analyzed_genomes_longer <- df_analyzed_genomes %>% select(accession, Domain:Genus) %>%
	pivot_longer(cols = -accession, names_to = "tax_level", values_to = "tax") %>%
	select(-tax_level)

df_trees_ar <- df_taxa_summarized %>% filter(Domain == "d__Archaea", n_genomes >= 50) %>%
	select(tax) %>%
	left_join(df_analyzed_genomes_longer, by = "tax") %>%
	group_by(tax) %>%
	nest() %>%
	mutate(
		set_diff = map(data, ~ setdiff(gtdb_tree_ar$tip.label, .$accession)),
		trees = map(set_diff, ~ drop.tip(gtdb_tree_ar, .)),
		faith_pd = map(trees, ~ sum(.$edge.length)),
		n_leaves = map(trees, ~ length(.$tip.label))
	) %>%
	ungroup() %>%
	unnest(c(faith_pd, n_leaves)) %>%
	select(tax, faith_pd, n_leaves)

df_trees_bac <- df_taxa_summarized %>% filter(Domain == "d__Bacteria", n_genomes >= 50) %>%
	select(tax) %>%
	left_join(df_analyzed_genomes_longer, by = "tax") %>%
	group_by(tax) %>%
	nest() %>%
	mutate(
		set_diff = map(data, ~ setdiff(gtdb_tree_bac$tip.label, .$accession)),
		trees = map(set_diff, ~ drop.tip(gtdb_tree_bac, .)),
		faith_pd = map(trees, ~ sum(.$edge.length)),
		n_leaves = map(trees, ~ length(.$tip.label))
	) %>%
	ungroup() %>%
	unnest(c(faith_pd, n_leaves)) %>%
	select(tax, faith_pd, n_leaves)

df_taxa_summarized_stats_tree <-
	left_join(df_taxa_summarized, rbind(df_trees_ar, df_trees_bac), by = "tax")

prok <- df_taxa_summarized_stats_tree %>% filter(tax_level == "Domain") %>%
	select(n_genomes, faith_pd, n_leaves) %>%
	colSums()

df_taxa_summarized_stats_tree <- dplyr::bind_rows(
		tibble(
			tax = "k__Prokaryotes",
			tax_level = "Kingdom",
			n_genomes = prok["n_genomes"],
			faith_pd = prok["faith_pd"],
			n_leaves = prok["n_leaves"]
		),
		df_taxa_summarized_stats_tree 
	) %>%
	left_join(df_analyzed_genomes_stats, by = "tax")

write.table(df_taxa_summarized_stats_tree, "data/tsv/taxa_summarized_stats_pd.tsv", append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE, na = "")

write.table(
	df_taxa_summarized_stats_tree %>% filter(n_genomes >= 50),
	"data/tsv/taxa_summarized_50_stats_pd.tsv", append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE, na = ""
)
