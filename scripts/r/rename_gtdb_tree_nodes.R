source('scripts/r/functions.R')

df_genomes_analyzed <- read_genomes_analyzed()
df_taxa <- read.csv("data/tsv/taxa_summarized_50.tsv", header = TRUE, sep = "\t", quote = "") %>% as_tibble()

# Archaea

gtdb_tree_ar <- read_gtdb_tree_ar()

gtdb_tree_analyzed_ar <- gtdb_tree_analyzed_ar()
df_gtdb_tree_node_labeled_ar <- makeNodeLabel(gtdb_tree_analyzed_ar) %>% as_tibble()

df_mrca_ar <- df_taxa %>% filter(Domain == "d__Archaea") %>%
	select(tax, tax_level) %>%
	mutate(
		genomes = map(tax,
			~ df_genomes_analyzed %>%
				filter(Domain == .x | Phylum == .x | Class == .x | Order == .x | Family == .x | Genus == .x) %>%
				select(accession) %>% pull()
		),
		mrca = map(genomes, ~ MRCA(gtdb_tree_analyzed_ar, .x))
	) %>%
	select(-genomes) %>%
	unnest(c(mrca))

df_mrca_ar$tax_level <- factor(df_mrca_ar$tax_level, levels = c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Domain"), ordered = TRUE)
df_mrca_reduced_ar <- df_mrca_ar %>% arrange(tax_level) %>%
	distinct(mrca, .keep_all = TRUE) %>%
	arrange(mrca)
df_gtdb_tree_node_labeled_ar$label[which(df_gtdb_tree_node_labeled_ar$node %in% df_mrca_reduced_ar$mrca)] <- df_mrca_reduced_ar$tax

write.tree(as.phylo(df_gtdb_tree_node_labeled_ar), "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree")

# Bacteria

gtdb_tree_bac <- read_gtdb_tree_bac()

gtdb_tree_analyzed_bac <- gtdb_tree_analyzed_bac()
df_gtdb_tree_node_labeled_bac <- makeNodeLabel(gtdb_tree_analyzed_bac) %>% as_tibble()

df_mrca_bac <- df_taxa %>% filter(Domain == "d__Bacteria") %>%
	select(tax, tax_level) %>%
	mutate(
		genomes = map(tax,
			~ df_genomes_analyzed %>%
				filter(Domain == .x | Phylum == .x | Class == .x | Order == .x | Family == .x | Genus == .x) %>%
				select(accession) %>% pull()
		),
		mrca = map(genomes, ~ MRCA(gtdb_tree_analyzed_bac, .x))
	) %>%
	select(-genomes) %>%
	unnest(c(mrca))

df_mrca_bac$tax_level <- factor(df_mrca_bac$tax_level, levels = c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Domain"), ordered = TRUE)
df_mrca_reduced_bac <- df_mrca_bac %>% arrange(tax_level) %>%
	distinct(mrca, .keep_all = TRUE) %>%
	arrange(mrca)
df_gtdb_tree_node_labeled_bac$label[which(df_gtdb_tree_node_labeled_bac$node %in% df_mrca_reduced_bac$mrca)] <- df_mrca_reduced_bac$tax

write.tree(as.phylo(df_gtdb_tree_node_labeled_bac), "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree")
