# Define functions used in the analysis

# Load tidyverse

library(tidyverse)
library(viridis)
library(igraph)
library(parallel)
library(plyranges)
library(ggupset)
library(patchwork)
library(Biostrings)
library(DT)
library(data.table)
library(dtplyr)
library(RSQLite)
library(ggstar)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(tidytree)
library(treeio)
library(ggnetwork)
library(plotly)
library(htmlwidgets)
library(ape)
library(scales)
library(GGally)
library(hexbin)
library(coRdon)
library(ggpointdensity)

# ggplot setting

## Chenge default setting

theme_set(
    theme_bw() +
    theme(
        text = element_text(family = "Arial", size = 12),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.title.y = element_text(face = "bold", size = 13),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        plot.caption = element_text(hjust = 0.5, size = 10),
    )
)

## Color blind Palette
### 8 colors
cbp8 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
### 12 colors
cbp12 <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

# Genome retrieval

## Load GTDB genome metadata tables
## Only one modification from raw files: GTDB tax path is splited into separate columns

load_gtdb_meta <- function() {
	df_bac <- read.csv("raw/gtdb/bac120_metadata_r202.tsv", header = TRUE, sep = "\t", quote = "")
	df_ar <- read.csv("raw/gtdb/ar122_metadata_r202.tsv", header = TRUE, sep = "\t", quote = "")
	df_ar_tax <- df_ar %>% as_tibble() %>% mutate(gtdb_taxonomy_copy = gtdb_taxonomy) %>%
		separate(gtdb_taxonomy_copy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", extra = "drop")
	df_bac_tax <- df_bac %>% as_tibble() %>% mutate(gtdb_taxonomy_copy = gtdb_taxonomy) %>%
		separate(gtdb_taxonomy_copy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", extra = "drop")
	return(list(df_bac_tax, df_ar_tax))
}

## Retrieve representative high-quality genomes

rep_hq_genomes <- function(df_bac_tax, df_ar_tax) {
	df_ar_hq <- df_ar_tax %>% filter(gtdb_representative == "t", checkm_completeness >= 95, checkm_contamination <= 5)
	df_bac_hq <- df_bac_tax %>% filter(gtdb_representative == "t", checkm_completeness >= 95, checkm_contamination <= 5)
	return(list(df_bac_hq, df_ar_hq))
}

## Retrieve CPR representative genomes for precise estimation of completeness and coverage

cpr_rep_genomes <- function(df_bac_tax) {
	# Retrieve all p__Patescibacteria representative genomes
	df_cpr <- df_bac_tax %>% filter(Phylum == "p__Patescibacteria", gtdb_representative == "t")
	return(df_cpr)
}

## Load CheckM result using CPR marker set

cpr_checkm <- function() {
	df_cpr_cm <- read.csv("data/gtdb/cpr_completeness_contamination.tsv", header = TRUE, sep = "\t", quote = "") %>% as_tibble()
	return(df_cpr_cm)
}

## Retrieve CPR representative high-quality (based on CPR marker set) genomes

cpr_hq_rep_genomes <- function(df_cpr_cm) {
	df_cpr_hq <- df_cpr_cm %>% filter(cpr_completeness >= 95, cpr_contamination <= 5)
	return(df_cpr_hq)
}

## Merge CPR representative high-quality (based on CPR marker set) genomes and others

hq_genomes_merged <- function(df_bac_hq, df_ar_hq, df_cpr_hq) {
	# Retrieve representative high-quality CPR genome rows from GTDB metadata table
	df_bac_tax_cpr <- column_to_rownames(df_bac_tax, "accession")[df_cpr_hq$accession,] %>% rownames_to_column("accession") %>% as_tibble()
	df_bac_tax_cpr2 <- cbind(df_bac_tax_cpr, df_cpr_hq[c("cpr_completeness", "cpr_contamination")]) %>% as_tibble()
	# Retrieve representative high=quality genomes other than CPR
	df_bac_tax_non_cpr <- df_bac_hq %>% filter(Phylum != "p__Patescibacteria") %>% mutate(cpr_completeness = NA, cpr_contamination = NA)
	# Merge with CPRs and generate final genome data set to analyze
	df_bac_hq_merged <- rbind(df_bac_tax_non_cpr, df_bac_tax_cpr2) %>% as_tibble()
	# Archaea was not modified
	df_ar_hq_merged <- df_ar_hq %>% mutate(cpr_completeness = NA, cpr_contamination = NA)
	return(list(df_bac_hq_merged, df_ar_hq_merged))
}

## Load data/tsv/genomes_analyzed_ar.tsv and data/tsv/genomes_analyzed_bac.tsv

read_genomes_analyzed <- function() {
	df_genomes_analyzed <- rbind(
		read.csv("data/tsv/genomes_analyzed_ar.tsv", header = TRUE, sep = "\t", quote = "") %>% as_tibble(),
		read.csv("data/tsv/genomes_analyzed_bac.tsv", header = TRUE, sep = "\t", quote = "") %>% as_tibble()
	)
	return(df_genomes_analyzed)
}

## Summarize taxa

summarize_taxa <- function(each_tax) {
	df_each_tax_n_genomes <- df_genomes_analyzed %>% select(Domain:Species) %>%
		mutate(tax_level = each_tax,tax = .data[[each_tax]]) %>%
		group_by(tax, tax_level) %>%
		summarize(n_genomes = dplyr::n())
	df_each_tax_taxpath <- df_genomes_analyzed %>% select(Domain:.data[[each_tax]]) %>%
		mutate(tax = .data[[each_tax]]) %>%
		distinct()
	df_each_tax_summarized <- inner_join(df_each_tax_n_genomes, df_each_tax_taxpath, by = "tax")
	return(df_each_tax_summarized)
}

## Retrieve genomes of only taxa with >= 50 genomes, and generate final genome data set

taxa_50_genomes <- function(df_bac_hq_merged, df_ar_hq_merged) {
	# Bacteria
	bac_p_50 <- df_bac_hq_merged %>% group_by(Phylum) %>% summarize(tax_genomes = dplyr::n()) %>% filter(tax_genomes >= 50) %>% dplyr::rename(taxon = Phylum) %>% mutate(tax_level = "Phylum")
	bac_c_50 <- df_bac_hq_merged %>% group_by(Class) %>% summarize(tax_genomes = dplyr::n()) %>% filter(tax_genomes >= 50) %>% dplyr::rename(taxon = Class) %>% mutate(tax_level = "Class")
	bac_o_50 <- df_bac_hq_merged %>% group_by(Order) %>% summarize(tax_genomes = dplyr::n()) %>% filter(tax_genomes >= 50) %>% dplyr::rename(taxon = Order) %>% mutate(tax_level = "Order")
	bac_f_50 <- df_bac_hq_merged %>% group_by(Family) %>% summarize(tax_genomes = dplyr::n()) %>% filter(tax_genomes >= 50) %>% dplyr::rename(taxon = Family) %>% mutate(tax_level = "Family")
	bac_g_50 <- df_bac_hq_merged %>% group_by(Genus) %>% summarize(tax_genomes = dplyr::n()) %>% filter(tax_genomes >= 50) %>% dplyr::rename(taxon = Genus) %>% mutate(tax_level = "Genus")
	bac_50_taxon <- rbind(bac_p_50, bac_c_50, bac_o_50, bac_f_50, bac_g_50)
	# Archaea
	ar_p_50 <- df_ar_hq_merged %>% group_by(Phylum) %>% summarize(tax_genomes = dplyr::n()) %>% filter(tax_genomes >= 50) %>% dplyr::rename(taxon = Phylum) %>% mutate(tax_level = "Phylum")
	ar_c_50 <- df_ar_hq_merged %>% group_by(Class) %>% summarize(tax_genomes = dplyr::n()) %>% filter(tax_genomes >= 50) %>% dplyr::rename(taxon = Class) %>% mutate(tax_level = "Class")
	ar_o_50 <- df_ar_hq_merged %>% group_by(Order) %>% summarize(tax_genomes = dplyr::n()) %>% filter(tax_genomes >= 50) %>% dplyr::rename(taxon = Order) %>% mutate(tax_level = "Order")
	ar_f_50 <- df_ar_hq_merged %>% group_by(Family) %>% summarize(tax_genomes = dplyr::n()) %>% filter(tax_genomes >= 50) %>% dplyr::rename(taxon = Family) %>% mutate(tax_level = "Family")
	ar_g_50 <- df_ar_hq_merged %>% group_by(Genus) %>% summarize(tax_genomes = dplyr::n()) %>% filter(tax_genomes >= 50) %>% dplyr::rename(taxon = Genus) %>% mutate(tax_level = "Genus")
	ar_50_taxon <- rbind(ar_p_50, ar_c_50, ar_o_50, ar_f_50, ar_g_50)
	# Select genomes of these taxa from representative hq data set, and generate final data set for analysis
	df_bac_p <- bac_p_50 %>% select(taxon) %>% left_join(mutate(df_bac_hq_merged, taxon = Phylum), by = "taxon")
	df_bac_c <- bac_c_50 %>% select(taxon) %>% left_join(mutate(df_bac_hq_merged, taxon = Class), by = "taxon")
	df_bac_o <- bac_o_50 %>% select(taxon) %>% left_join(mutate(df_bac_hq_merged, taxon = Order), by = "taxon")
	df_bac_f <- bac_f_50 %>% select(taxon) %>% left_join(mutate(df_bac_hq_merged, taxon = Family), by = "taxon")
	df_bac_g <- bac_g_50 %>% select(taxon) %>% left_join(mutate(df_bac_hq_merged, taxon = Genus), by = "taxon")
	df_bac_analyzed <- rbind(df_bac_p, df_bac_c, df_bac_o, df_bac_f, df_bac_g) %>% select(-taxon) %>% distinct()
	df_ar_p <- ar_p_50 %>% select(taxon) %>% left_join(mutate(df_ar_hq_merged, taxon = Phylum), by = "taxon")
	df_ar_c <- ar_c_50 %>% select(taxon) %>% left_join(mutate(df_ar_hq_merged, taxon = Class), by = "taxon")
	df_ar_o <- ar_o_50 %>% select(taxon) %>% left_join(mutate(df_ar_hq_merged, taxon = Order), by = "taxon")
	df_ar_f <- ar_f_50 %>% select(taxon) %>% left_join(mutate(df_ar_hq_merged, taxon = Family), by = "taxon")
	df_ar_g <- ar_g_50 %>% select(taxon) %>% left_join(mutate(df_ar_hq_merged, taxon = Genus), by = "taxon")
	df_ar_analyzed <- rbind(df_ar_p, df_ar_c, df_ar_o, df_ar_f, df_ar_g) %>% select(-taxon) %>% distinct()
	return(list(bac_50_taxon, ar_50_taxon, df_bac_analyzed, df_ar_analyzed))
}

# Proteins

## Length (aa)

protein_length <- function() {
	df_protein_len <- read.csv("data/tsv/proteins_raw_len_aa.tsv", header = FALSE, sep = "\t", quote = "") %>% as_tibble() %>% dplyr::rename(prodigal_acc = "V1", protein_aa = "V2")
	return(df_protein_len)
}

protein_len_density <- function(df_protein_len) {
	p_protein_len <- df_protein_len %>% ggplot(aes(x = protein_aa)) +
		ggtitle("Protein length distribution") +
		xlab("Protein length (aa)") +
		ylab("Density") +
		geom_density() +
		scale_x_log10(
			breaks = scales::trans_breaks("log10", function(x) 10^x),
   			labels = scales::trans_format("log10", scales::math_format(10^.x))
		) +
		annotation_logticks(sides = "b") +
		theme_bw() +
		theme(plot.title = element_text(hjust = 0.5), text = element_text(family = "Arial", size = 12))
	return(p_protein_len)
}

sequence_length <- function(in_file, seqtype) {
	if (seqtype == "na") {
		seq <- readDNAStringSet(in_file, "fasta")
	} else if (seqtype == "aa") {
		seq <- readAAStringSet(in_file, "fasta")
	}
	return(seq)
}

sequence_len_distribution <- function(each_seq) {
	p_protein_len <- seq %>% width() %>% as_tibble() %>%
		ggplot(aes(x = value)) +
		geom_density() +
		scale_x_log10(
			breaks = scales::trans_breaks("log10", function(x) 10^x),
   			labels = scales::trans_format("log10", scales::math_format(10^.x))
		) +
		annotation_logticks(sides = "b")
	return(p_protein_len)
}
# Phylogenetic tree

## Load GTDB representative (trimmed) tree

read_gtdb_tree_ar <- function () {
	gtdb_tree_ar <- read.tree("raw/gtdb/ar122_r202_ete3fmt5.tree")
	return(gtdb_tree_ar)
}

read_gtdb_tree_bac <- function () {
	gtdb_tree_bac <- read.tree("raw/gtdb/bac120_r202_ete3fmt5.tree")
	return(gtdb_tree_bac)
}

read_gtdb_tree_ar_fmt1 <- function () {
	gtdb_tree_ar <- read.tree("raw/gtdb/ar122_r202_ete3fmt1.tree")
	return(gtdb_tree_ar)
}

read_gtdb_tree_bac_fmt1 <- function () {
	gtdb_tree_bac <- read.tree("raw/gtdb/bac120_r202_ete3fmt1.tree")
	return(gtdb_tree_bac)
}

## Retrieve leaves of analyzed genomes

gtdb_tree_analyzed_ar <- function() {
	gtdb_tree_analyzed_ar <- drop.tip(
		gtdb_tree_ar,
		setdiff(
			gtdb_tree_ar$tip.label,
			df_genomes_analyzed %>% filter(Domain == "d__Archaea") %>% select(accession) %>% pull()
		)
	)
}

gtdb_tree_analyzed_bac <- function() {
	gtdb_tree_analyzed_bac <- drop.tip(
		gtdb_tree_bac,
		setdiff(
			gtdb_tree_bac$tip.label,
			df_genomes_analyzed %>% filter(Domain == "d__Bacteria") %>% select(accession) %>% pull()
		)
	)
}

# COG2020

## Load db/cog2020/cog-20.cog.csv

read_cog2020_csv <- function(n) {
	df_cog2020 <- read.csv("db/cog2020/cog-20.cog.csv", header = FALSE, sep = ",", quote = "", nrows = n) %>% as_tibble() %>%
		dplyr::rename(gene_id = "V1", ncbi_assembly_id = "V2", protein_id = "V3", protein_length = "V4", footprint_on_protein = "V5",
			len_footprint_on_protein = "V6", cog_id = "V7", reserved = "V8", cog_membership_class = "V9", PSI_BLAST_bitscore = "V10",
			PSI_BLAST_eval = "V11", cog_profile_len = "V12", footprint_on_cog = "V13") %>% as_tibble()
	return(df_cog2020)
}

## Load db/cog2020_fixed/cog-20.cog.tsv

read_cog2020_tsv <- function() {
	df_cog2020tsv <- read.csv("db/cog2020_fixed/cog-20.cog.tsv", header = TRUE, sep = "\t", quote = "") %>% as_tibble()
	return(df_cog2020tsv)
}

read_cog2020_def_fixed_funcat_text2 <- function () {
	df_fun20tab <- read.csv("db/cog2020/fun-20.tab", header = FALSE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>% as_tibble()
	colnames(df_fun20tab) <- c("cog_fun_cat","cog_fun_cat_color","cog_fun_cat_text")
	v_fun20tab <- df_fun20tab$cog_fun_cat_text
	names(v_fun20tab) <- df_fun20tab$cog_fun_cat
	df_cog2020_def <- read.csv("db/cog2020_fixed/cog-20.def.utf8.tsv", header = FALSE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>% as_tibble()
	colnames(df_cog2020_def) <- c("cog_id", "cog_fun_cat", "cog_name", "cog_gene", "cog_fun_path", "cog_pubmed_id", "cog_pdb_id")
	df_cog2020_def <- df_cog2020_def %>%
		mutate(
			data = strsplit(cog_fun_cat, ""),
			cog_fun_cat_text = map(data, ~ paste(paste(.x, v_fun20tab[.x], sep = ", "), collapse = ";"))
		) %>%
		select(-data) %>%
		unnest(c(cog_fun_cat_text)) %>%
		mutate(
			gene_name_custom = map2(cog_id, cog_gene, ~
				if (is.na(.y)) {
					.x
				} else {
					.y
				}
			)
		) %>% unnest(c(gene_name_custom))
	return(df_cog2020_def)
}

## Transform cog table: split footprint columns

cog2020_split_footprint <- function(n) {
	each_fp <- strsplit(df_cog2020[n,]$footprint_on_protein, "-")[[1]] %>% as.numeric()
	each_fc <- strsplit(df_cog2020[n,]$footprint_on_cog, "-")[[1]] %>% as.numeric()
	cbind(df_cog2020[n,], tibble(protein_start = each_fp[1], protein_end = each_fp[2], cog_start = each_fc[1], cog_end = each_fc[2],
		len_footprint_on_cog = each_fc[2] - each_fc[1] + 1))
}

## Load db/cog2020/cog-20.def.tab

read_cog2020_def <- function () {
	df_cog2020_def <- read.csv("db/cog2020/cog-20.def.tab", header = FALSE, sep = "\t", quote = "") %>% as_tibble()
	colnames(df_cog2020_def) <- c("cog_id", "cog_fun_cat", "cog_name", "cog_gene", "cog_fun_path", "cog_pubmed_id", "cog_pdb_id")
	return(df_cog2020_def)
}

read_cog2020_def_na <- function () {
	df_cog2020_def <- read.csv("db/cog2020/cog-20.def.tab", header = FALSE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>% as_tibble()
	colnames(df_cog2020_def) <- c("cog_id", "cog_fun_cat", "cog_name", "cog_gene", "cog_fun_path", "cog_pubmed_id", "cog_pdb_id")
	return(df_cog2020_def)
}

read_cog2020_def_fixed <- function () {
	df_cog2020_def <- read.csv("db/cog2020_fixed/cog-20.def.utf8.tsv", header = FALSE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>% as_tibble()
	colnames(df_cog2020_def) <- c("cog_id", "cog_fun_cat", "cog_name", "cog_gene", "cog_fun_path", "cog_pubmed_id", "cog_pdb_id")
	return(df_cog2020_def)
}

## Load db/cog2020_fixed/cog-20.seqname_cogid.tsv

read_cog2020_seqname_cogid <- function () {
	df_cog2020_seqname_cogid_multi_cogids <- read.csv("db/cog2020_fixed/cog-20.seqname_cogid.tsv", header = TRUE, sep = "\t", quote = "") %>% as_tibble()
	return(df_cog2020_seqname_cogid_multi_cogids)
}

## Join vector of COGs with db/cog2020/cog-20.def.tab

join_cog2020_def <- function (cogs) {
	df_cogs <- tibble(cog_id = cogs) %>% inner_join(df_cog2020_def, by = "cog_id")
	return(df_cogs)
}

# KEGG

## Load db/kegg/genes/ko/ko_genes.list

read_ko_genes_list <- function () {
	df_ko_genes_list <- read.csv("db/kegg/genes/ko/ko_genes.list", header = FALSE, sep = "\t", quote = "") %>% as_tibble()
	colnames(df_ko_genes_list) <- c("ko_id", "gene_id")
	return(df_ko_genes_list)
}

read_ko_lists <- function(ko_list, ko_list_header) {
	df_ko_list <- read.csv(ko_list, header = FALSE, sep = "\t", quote = "") %>% as_tibble()
	colnames(df_ko_list) <- ko_list_header
	return(df_ko_list)
}

read_ko_def <- function () {
	df_ko_def <- read.csv("db/kegg_fixed/kegg_ko_def.tsv", header = TRUE, sep = "\t", quote = "") %>% as_tibble()
	return(df_ko_def)
}

read_ko_def_na <- function () {
	df_ko_def <- read.csv("db/kegg_fixed/kegg_ko_def.tsv", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
		as_tibble()
	return(df_ko_def)
}

read_ko_def_na2 <- function () {
	df_ko_def <- read.csv("db/kegg_fixed/kegg_ko_def.tsv", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
		as_tibble() %>%
		mutate(
			gene_name_custom = map2(ko_id, ko_symbol, ~
				if (is.na(.y)) {
					.x
				} else {
					.y
				}
			)
		) %>% unnest(c(gene_name_custom))
	return(df_ko_def)
}

# MMSeqs2

## Load search result file

mmseqs_search_result <- function(mmseqs_result) {
	each_mmseqs_search_result <- read.csv(mmseqs_result, header = FALSE, sep = "\t", quote = "", col.names = c("query","target","fident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits","qlen","qcov","tlen","tcov")) %>% as_tibble()
	return(each_mmseqs_search_result)
}

## Load COG2020 joind MMSeqs search result file

read_mmseqs_cog2020 <- function(mmseqs_cog2020) {
	df_mmseqs_cog2020 <- read.csv(mmseqs_cog2020, header = TRUE, sep = "\t", quote = "") %>% as_tibble()
	return(df_mmseqs_cog2020)
}

## Estimation of search results

estimate_mmseqs_search_result <- function(mmseqs_search_tsv) {
	df_mmseqs_search_result <- mmseqs_search_result(mmseqs_search_tsv)
	n_hits <- df_mmseqs_search_result %>% nrow()
	n_query <- df_mmseqs_search_result %>% select(query) %>% distinct() %>% nrow()
	n_query_target <- df_mmseqs_search_result %>% select(query,target) %>% distinct() %>% nrow()
	n_hits_filtered <- df_mmseqs_search_result %>% filter(qcov >= 0.8, tcov >= 0.8, bits >= 50) %>% nrow()
	n_query_filtered <- df_mmseqs_search_result %>% filter(qcov >= 0.8, tcov >= 0.8, bits >= 50) %>% select(query) %>% distinct() %>% nrow()
	tibble(tsv = mmseqs_search_tsv, duplicated_hits = nrow(df_mmseqs_search_result) != n_query_target, n_hit_query = n_query, n_hit_query_filtered = n_query_filtered, n_hits = n_hits, n_hits_filtered = n_hits_filtered)
}

estimate_mmseqs_search_result_cdd <- function(mmseqs_search_cdd_tsv) {
	df_mmseqs_search_cdd_result <- mmseqs_search_result(mmseqs_search_cdd_tsv)
	df_mmseqs_search_cdd_result <- df_mmseqs_search_cdd_result[grep("^COG", df_mmseqs_search_cdd_result$target),]
	n_hits <- df_mmseqs_search_cdd_result %>% nrow()
	n_query <- df_mmseqs_search_cdd_result %>% select(query) %>% distinct() %>% nrow()
	n_query_target <- df_mmseqs_search_cdd_result %>% select(query,target) %>% distinct() %>% nrow()
	n_hits_filtered <- df_mmseqs_search_cdd_result %>% filter(tcov >= 0.8, bits >= 50) %>% nrow()
	n_query_filtered <- df_mmseqs_search_cdd_result %>% filter(tcov >= 0.8, bits >= 50) %>% select(query) %>% distinct() %>% nrow()
	tibble(tsv = mmseqs_search_cdd_tsv, duplicated_hits = nrow(df_mmseqs_search_cdd_result) != n_query_target, n_hit_query = n_query, n_hit_query_filtered = n_query_filtered, n_hits = n_hits, n_hits_filtered = n_hits_filtered)
}

estimate_diamond_search_result <- function(diamond_search_tsv) {
	df_diamond_search_result <- mmseqs_search_result(diamond_search_tsv)
	n_hits <- df_diamond_search_result %>% nrow()
	n_query <- df_diamond_search_result %>% select(query) %>% distinct() %>% nrow()
	n_query_target <- df_diamond_search_result %>% select(query,target) %>% distinct() %>% nrow()
	n_hits_filtered <- df_diamond_search_result %>% filter(qcov >= 80, tcov >= 80, bits >= 50) %>% nrow()
	n_query_filtered <- df_diamond_search_result %>% filter(qcov >= 80, tcov >= 80, bits >= 50) %>% select(query) %>% distinct() %>% nrow()
	tibble(tsv = diamond_search_tsv, duplicated_hits = nrow(df_diamond_search_result) != n_query_target, n_hit_query = n_query, n_hit_query_filtered = n_query_filtered, n_hits = n_hits, n_hits_filtered = n_hits_filtered)
}

## Filter MMSeqs search results

filter_mmseqs <-  function() {
	rng_search_result <- df_search_result %>%
		dplyr::rename(start = "qstart", end = "qend") %>%
		nest(data = -query) %>%
		mutate(
			rang = map(data, ~ as_iranges(.)),
			reduced = map(rang, ~ reduce_ranges(.)),
			paired = map2(reduced, rang, ~ as_tibble(pair_overlaps(.x, .y)))
		) %>%
		rowwise() %>%
		mutate(is_multi_domain_2 = length(reduced) > 1)
	df_search_result_filtered <- rng_search_result %>% select(query, paired, is_multi_domain_2) %>%
		unnest(cols = paired) %>%
		nest(data = c(-query, -ranges.x.start, -ranges.x.end)) %>%
		mutate(
			filtered = map(data, ~ filter(., bits >= max(bits) * 0.9))
		) %>%
		select(-data) %>%
		unnest(cols = filtered)
#	df_search_result_filtered <- df_search_result %>%
#		nest(data = -query) %>%
#		rowwise() %>%
#		mutate(
#			max_bit = max(data$bits),
#			filtered = list(data %>% filter(bits >= 0.8 * max_bit))
#		) %>%
#		ungroup() %>%
#		select(-data) %>%
#		unnest(c(filtered))
	return(df_search_result_filtered)
}

# DIAMOND

## Load search result file

diamond_search_result <- function(diamond_result) {
	each_diamond_search_result <- read.csv(diamond_result, header = FALSE, sep = "\t", quote = "", col.names = c("query","target","fident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits","qlen","q_cov","tlen","t_cov")) %>%
	as_tibble() %>%
	mutate(qcov = q_cov / 100, tcov = t_cov / 100) %>%
	select("query","target","fident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits","qlen","qcov","tlen","tcov")
	return(each_diamond_search_result)
}

# Assign OGs

## Select OGs in each range

select_ogs_each_range <- function(rang) {
	df_overlaps <- join_overlap_intersect(each_protein_rgs, rgs[rang]) %>% arrange(dplyr::desc(bits)) %>% filter(bits >= max(bits)*0.8) %>% as.data.frame()
	ovlp_start <- as.data.frame(rgs[rang])$start
	ovlp_end <- as.data.frame(rgs[rang])$end
	ovlp_width <- as.data.frame(rgs[rang])$width
	ovlp_ogs <- df_overlaps$target %>% unique() %>% sort() %>% paste(collapse =";")
	if (length(unique(df_overlaps$target)) > 1) {
		og_ovlp <- TRUE
	} else {
		og_ovlp <- FALSE
	}
	df_overlaps_grouped <- df_overlaps %>% group_by(target) %>%
		summarize(
			fident = dplyr::first(fident),
			alnlen = dplyr::first(alnlen),
			mismatch = dplyr::first(mismatch),
			gapopen = dplyr::first(gapopen),
			protein_start = dplyr::first(start),
			protein_end = dplyr::first(end),
			og_start = dplyr::first(tstart),
			og_end = dplyr::first(tend),
			evalue = dplyr::first(evalue),
			bits = dplyr::first(bits),
			protein_len = dplyr::first(qlen),
			protein_cov = dplyr::first(qcov),
			og_len = dplyr::first(tlen),
			og_cov = dplyr::first(tcov),
		) %>% dplyr::rename(og_assigned = "target") %>% arrange(desc(bits))
	df_range_each <- data.frame(protein_acc = each_protein, og_assigned = df_overlaps_grouped$og_assigned, og_overlaps = og_ovlp, ovlp_ogs = ovlp_ogs, og_multi_domain = og_multi_domain, protein_ovlp_start = ovlp_start, protein_ovlp_end = ovlp_end, protein_ovlp_width = ovlp_width) %>% cbind(select(df_overlaps_grouped, -og_assigned))
}

## Assign OGs from mmseqs results

assign_cog <-  function() {
	df_mmseqs_cog2020_assigned <- df_filtered_og_joined %>%
		nest(data = c(-query, -ranges.x.start, -ranges.x.end, -ranges.x.width, -is_multi_hits_h)) %>%
		mutate(
			max_bit = map(data, ~ max(.$bits)),
			is_multi_domain = map(data, ~ sum(.$is_multi_domain) > 0),
			is_multi_hits_v = map(data, ~ nrow(distinct(., cog_ids)) > 1),
			cog_id = map(data, ~ arrange(distinct(., cog_id))),
			cog_ids = map(data, ~ paste(unique(sort(.$cog_id)), collapse = ";"))
		) %>%
		select(-data) %>%
		unnest(cols = c(max_bit, is_multi_domain, cog_id, cog_ids, is_multi_hits_v))
	return(df_mmseqs_cog2020_assigned)
}

assign_ko <-  function() {
	df_mmseqs_ko_assigned <- df_filtered_og_joined %>%
		nest(data = c(-query, -ranges.x.start, -ranges.x.end, -ranges.x.width, -is_multi_hits_h)) %>%
		mutate(
			max_bit = map(data, ~ max(.$bits)),
			is_multi_domain = map(data, ~ sum(.$is_multi_domain) > 0),
			is_multi_hits_v = map(data, ~ nrow(distinct(., ko_ids)) > 1),
			ko_id = map(data, ~ arrange(distinct(., ko_id))),
			ko_ids = map(data, ~ paste(unique(sort(.$ko_id)), collapse = ";"))
		) %>%
		select(-data) %>%
		unnest(cols = c(max_bit, is_multi_domain, ko_id, ko_ids, is_multi_hits_v))
	return(df_mmseqs_ko_assigned)
}

load_proteins_raw_gtdbtax_og_tidy <- function(path_tidy) {
	dt_og_tidy <- fread(path_tidy, header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL"))
	return(dt_og_tidy)
}

# Retrieve genes

retrieve_gene_by_cog <- function(each_cog) {
	con <- dbConnect(SQLite(), "data/sql/cpnet.db")
	res <- dbSendQuery(
		con,
		paste0(
			"SELECT * FROM proteins_raw_gtdbtax_cog2020 WHERE cog_id = '",
			each_cog,
			"'"
		)
	)
	df_each_cog <- dbFetch(res)
	dbClearResult(res)
dbListFields(con, "proteins_raw_gtdbtax_cog2020")
dbListTables(con)
}

retrieve_gene_by_ko <- function() {
}

# Draw Network

## Load data/tsv/network_edge_stats_*.tsv

read_network_edge_stats <- function(OG, edge_thre, edgetype) {
	if (OG == "cog") {
		df_network <- read.table(in_file, header = TRUE, sep = '\t', quote = "", na.strings=c("", "NULL")) %>%
			select(-cog_fun_cat_A:-n_cog_fun_cats_B) %>%
			filter(n_tax >= edge_thre, edge_type == edgetype) %>%
			as_tibble()
	} else if (OG == "kegg") {
		df_network <- fread(in_file, header = TRUE, sep = '\t', quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>%
			select(-ko_symbol_A:-n_ko_cogs_B, -n_ko_ecs_A, -n_ko_ecs_B, -n_ko_rns_A, -n_ko_rns_B, -n_ko_paths_num_specific_A, -n_ko_paths_num_specific_B, -n_ko_brs_specific_A, -n_ko_brs_specific_B) %>%
			filter(n_tax >= edge_thre, edge_type == edgetype)
	}
	return(df_network)
}

## Convert graph object to data frame
make_ggnetwork_df <- function(g_network) {
  set.seed(1)
  df_g_network <- g_network %>%
    ggnetwork(layout = igraph::layout_with_fr(g_network)) %>%
    as_tibble() %>%
    mutate(
      x_same = x == xend,
      y_same = y == yend,
      v_e = map2(x_same, y_same, ~ if (.x == TRUE & .y == TRUE) {"vertex"} else {"edge"})
    ) %>%
    unnest(c(v_e))
}

graph2df <- function(OG, col_name) {
	if (OG == "cog" & col_name == "cog_fun_cat") {
		set.seed(1)
		df_g_network <- g_network %>%
			ggnetwork(layout = igraph::layout_with_fr(g_network)) %>%
			as_tibble() %>%
			mutate(
				x_same = x == xend,
				y_same = y == yend,
				v_e = map2(x_same, y_same, ~ if (.x == TRUE & .y == TRUE) {"vertex"} else {"edge"}),
				v_shape = map(.data[[col_name]], ~ if (length(grep("(R|S)", .x)) >= 1) {"23"} else {"21"}),
				v_color = .data[[col_name]],
			) %>%
			unnest(c(v_e, v_shape))
	} else if (OG == "cog") {
		set.seed(1)
		df_g_network <- g_network %>%
			ggnetwork(layout = igraph::layout_with_fr(g_network)) %>%
			as_tibble() %>%
			mutate(
				x_same = x == xend,
				y_same = y == yend,
				v_e = map2(x_same, y_same, ~ if (.x == TRUE & .y == TRUE) {"vertex"} else {"edge"}),
				v_shape = map(.data[[col_name]], ~ if (is.na(.x)) {"23"} else {"21"}),
				v_color = .data[[col_name]]
			) %>%
			unnest(c(v_e, v_shape))
	} else if (OG == "kegg") {
		set.seed(1)
		df_g_network <- g_network %>%
			ggnetwork(layout = igraph::layout_with_fr(g_network)) %>%
			as_tibble() %>%
			mutate(
				x_same = x == xend,
				y_same = y == yend,
				v_e = map2(x_same, y_same, ~ if (.x == TRUE & .y == TRUE) {"vertex"} else {"edge"}),
				v_shape = map(.data[[col_name]], ~ if (is.na(.x)) {"23"} else {"21"}),
				v_color = .data[[col_name]]
			) %>%
			unnest(c(v_e, v_shape))
	}
	return(df_g_network)
}

graph2df_habitat <- function(OG, col_name) {
	if (OG == "cog" & col_name == "cog_fun_cat") {
		set.seed(1)
		df_g_network <- g_network %>%
			ggnetwork(layout = igraph::layout_with_fr(g_network)) %>%
			as_tibble() %>%
			mutate(
				x_same = x == xend,
				y_same = y == yend,
				v_e = map2(x_same, y_same, ~ if (.x == TRUE & .y == TRUE) {"vertex"} else {"edge"}),
				v_shape = map2(
					.data[[col_name]],
					node_type,
					~ if (length(grep("(R|S)", .x)) >= 1 & .y == "OG") {"23"}
						else if (.y == "Habitat") {"22"}
						else {"21"}
				),
				v_color = .data[[col_name]],
			) %>%
			unnest(c(v_e, v_shape))
	} else if (OG == "kegg") {
		set.seed(1)
		df_g_network <- g_network %>%
			ggnetwork(layout = igraph::layout_with_fr(g_network)) %>%
			as_tibble() %>%
			mutate(
				x_same = x == xend,
				y_same = y == yend,
				v_e = map2(x_same, y_same, ~ if (.x == TRUE & .y == TRUE) {"vertex"} else {"edge"}),
				v_shape = map(.data[[col_name]], ~ if (is.na(.x)) {"23"} else {"21"}),
				v_shape = map2(
					.data[[col_name]],
					node_type,
					~ if (is.na(.x) & .y == "OG") {"23"}
						else if (.y == "Habitat") {"22"}
						else {"21"}
				),
				v_color = .data[[col_name]]
			) %>%
			unnest(c(v_e, v_shape))
	}
	return(df_g_network)
}

graph2df_plotly <- function(OG, col_name) {
	if (OG == "cog") {
		df_g_network_plotly <- df_g_network %>%
			rowwise() %>%
			mutate(
				text = if (v_e == "vertex") {
						paste0(
							"<b>OG ID: </b>", name, "<br>",
							"<b>  COG name: </b>", cog_name, "<br>",
							"<b>  COG gene: </b>", cog_gene, "<br>",
							"<b>  COG category: </b>", cog_fun_cat, "<br>",
							"<b>  COG pathway: </b>", cog_fun_path, "<br>",
							"<b>  COG PubMed: </b>", cog_pubmed_id, "<br>",
							"<b>  COG PDB: </b>", cog_pdb_id
						)
					} else {
						""
					}
			)
	} else if (OG == "kegg") {
		df_g_network_plotly <- df_g_network %>%
			rowwise() %>%
			mutate(
				text = if (v_e == "vertex") {
						paste0(
							"<b>OG ID: </b>", name, "<br>",
							"<b>  KO symbol: </b>", ko_symbol, "<br>",
							"<b>  KO name: </b>", ko_name, "<br>",
							"<b>  COG ID: </b>", ko_cog, "<br>",
							"<b>  EC number: </b>", ko_ec, "<br>",
							"<b>  KEGG pathway ID: </b>", ko_path_num_specific, "<br>",
							"<b>  KEGG brite: </b>", ko_br_specific
						)
					} else {
						""
					}
			)
	}
	return(df_g_network_plotly)
}
## Select edge color

ggplot_network <- function(col_name, category) {
	if (category == "discrete") {
		p_network <- ggplot() +
			geom_segment(
				aes(x = x, y = y, xend = xend, yend = yend, color = .data[[col_name]], size = n_tax),
				data = df_g_network %>% filter(v_e == "edge")
			) +
			geom_point(
				aes(x = x, y = y, fill = v_color, shape = v_shape),
				size = 2.7,
				color = "white",
				data = df_g_network %>% filter(v_e == "vertex")
			) +
			scale_fill_viridis_d(na.value = "black") +
			scale_shape_manual(values = c(21, 23)) +
			scale_color_manual(values = cbp8, na.value = "grey50") +
			coord_fixed(ratio = 1) +
			guides(fill = "none", shape = "none")
	} else if (category == "continuous") {
		p_network <- ggplot() +
			geom_segment(
				aes(x = x, y = y, xend = xend, yend = yend, color = .data[[col_name]], size = n_tax),
				data = df_g_network %>% filter(v_e == "edge")
			) +
			geom_point(
				aes(x = x, y = y, fill = v_color, shape = v_shape),
				size = 2.7,
				color = "white",
				data = df_g_network %>% filter(v_e == "vertex")
			) +
			scale_fill_viridis_d() +
			scale_shape_manual(values = c(21, 23)) +
			scale_color_viridis_c(option = "plasma", begin = 0, end = 0.9) +
			coord_fixed(ratio = 1) +
			guides(fill = "none", shape = "none")
	}
	return(p_network)
}

ggplot_network_habitat <- function(col_name, category) {
	if (category == "discrete") {
		p_network <- ggplot() +
			geom_segment(
				aes(x = x, y = y, xend = xend, yend = yend, color = .data[[col_name]], size = n_tax),
				data = df_g_network %>% filter(v_e == "edge")
			) +
			geom_point(
				aes(x = x, y = y, fill = v_color, shape = v_shape),
				size = 2.7,
				color = "white",
				data = df_g_network %>% filter(v_e == "vertex")
			) +
			scale_fill_viridis_d(na.value = "black") +
			scale_shape_manual(values = c(21, 22, 23)) +
			scale_color_manual(values = cbp8, na.value = "grey50") +
			coord_fixed(ratio = 1) +
			guides(fill = "none", shape = "none")
	} else if (category == "continuous") {
		p_network <- ggplot() +
			geom_segment(
				aes(x = x, y = y, xend = xend, yend = yend, color = .data[[col_name]], size = n_tax),
				data = df_g_network %>% filter(v_e == "edge")
			) +
			geom_point(
				aes(x = x, y = y, fill = v_color, shape = v_shape),
				size = 2.7,
				color = "white",
				data = df_g_network %>% filter(v_e == "vertex")
			) +
			scale_fill_viridis_d(na.value = "black") +
			scale_shape_manual(values = c(21, 22, 23)) +
			scale_color_viridis_c(option = "plasma", begin = 0, end = 0.9, na.value = "black") +
			coord_fixed(ratio = 1) +
			guides(fill = "none", shape = "none")
	}
	return(p_network)
}

ggplotly_network <- function(vertex_color, edge_color) {
	vertex_colors <- sort(unique(
			df_g_network %>% select(.data[[vertex_color]]) %>% pull()
		))
	v_e_cols <- c(cbp8[1], cbp8[2], viridis(length(vertex_colors)))
	names(v_e_cols) <- c("FALSE", "TRUE", vertex_colors)
	v_e_color <- v_e_cols[
			levels(factor(c(
				df_g_network_plotly %>% select(.data[[edge_color]]) %>% pull(),
				df_g_network_plotly %>% select(.data[[vertex_color]]) %>% pull()
			)))
		]
	p <- df_g_network_plotly %>%
		ggplot(aes(x = x, y = y, xend = xend, yend = yend, text = text)) +
		geom_edges(
			aes(color = .data[[edge_color]], size = n_tax)
		) +
		geom_nodes(
			aes(color = .data[[vertex_color]], shape = v_shape),
			size = 2.5
		) +
		scale_shape_manual(values = c(16, 18)) +
		scale_color_manual(values = v_e_color, na.value = "black")
	P <- ggplotly(
			p + geom_nodetext(aes(label = name), size = 1) + coord_fixed(),
			tooltip = c("text")
		)
	return(P)
}

## Edge color: Categorical (TRUE or FALSE)

draw_network_edge_category_cog <- function () {
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend, color = e_color_tf, size = n_tax),
			data = df_g_network %>% filter(v_e == "edge")
		) +
		geom_point(
			aes(x = x, y = y, fill = cog_fun_cat, shape = v_shape),
			size = 2.7,
			color = "white",
			data = df_g_network %>% filter(v_e == "vertex")
		) +
		scale_fill_viridis_d() +
		scale_shape_manual(values = c(21, 23)) +
		scale_color_manual(values = cbp8) +
		guides(fill = "none", shape = "none")
	return(p_network)
}

draw_network_edge_category_kegg <- function () {
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend, color = e_color_tf, size = n_tax),
			data = df_g_network %>% filter(v_e == "edge")
		) +
		geom_point(
			aes(x = x, y = y, fill = ko_path_num_specific, shape = v_shape),
			size = 2.7,
			color = "white",
			data = df_g_network %>% filter(v_e == "vertex")
		) +
		scale_fill_viridis_d(na.value = "black") +
		scale_shape_manual(values = c(21, 23)) +
		scale_color_manual(values = cbp8) +
		guides(fill = "none", shape = "none")
	return(p_network)
}

draw_network_edge_category_kegg_edge <- function () {
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend, color = e_color_tf, size = n_tax),
			data = df_g_network %>% filter(v_e == "edge")
		) +
		scale_fill_viridis_d() +
		scale_color_manual(values = cbp8) +
		guides(fill = "none", shape = "none")
	return(p_network)
}

# Network stats

network_structure_tsv2wide <- function(df_net_struct_og, net_indicator) {
	df_net_struct_og_wide <- df_net_struct_og %>% select(tax, tax_level, Domain, edge_type, n_genomes_log, .data[[net_indicator]]) %>%
		pivot_wider(names_from = edge_type, values_from = .data[[net_indicator]]) %>%
		mutate(x = "contrapositive", xend = "cooccurring")
	return(df_net_struct_og_wide)
}

network_structure_2cols <- function(df_net_struct_og, df_net_struct_og_wide, net_indicator) {
	p_net_struct_og <- ggplot() +
		geom_segment(
			data = df_net_struct_og_wide,
			aes(x = x, y = contrapositive, xend = xend, yend = cooccurring, color = Domain),
			alpha = 0.4
		) +
		geom_point(
			data = df_net_struct_og,
			#aes(x = edge_type, y = .data[[net_indicator]], color = Domain, fill = n_genomes_log, shape = Domain),
			aes(x = edge_type, y = .data[[net_indicator]], fill = n_genomes_log, shape = Domain),
			color = "white",
			size = 2.5,
			alpha = 0.8
		) +
		scale_color_manual(values = cbp8, na.value = "black") +
		scale_fill_viridis(end = 0.95) +
		scale_shape_manual(values = c(22, 24), na.value = 21) +
		facet_grid(. ~ tax_level) +
		theme(
			axis.text.x = element_text(angle = 45, hjust = 1),
			panel.grid = element_blank()
		) +
		guides(
			shape = guide_legend(override.aes = list(fill = "black")),
			color = guide_legend(override.aes = list(alpha = 1))
		)
	return(p_net_struct_og)
}

network_structure_2cols_log <- function(df_net_struct_og, df_net_struct_og_wide, net_indicator) {
	p_net_struct_og <- ggplot() +
		geom_segment(
			data = df_net_struct_og_wide,
			aes(x = x, y = contrapositive, xend = xend, yend = cooccurring, color = Domain),
			alpha = 0.5
		) +
		geom_point(
			data = df_net_struct_og,
			#aes(x = edge_type, y = .data[[net_indicator]], color = Domain, shape = Domain),
			aes(x = edge_type, y = .data[[net_indicator]], fill = n_genomes_log, shape = Domain),
			color = "white",
			size = 2.5,
			alpha = 0.8
		) +
		annotation_logticks(sides = "l") +
		scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
		scale_color_manual(values = cbp8, na.value = "black") +
		scale_fill_viridis(end = 0.95) +
		scale_shape_manual(values = c(22, 24), na.value = 21) +
		facet_grid(. ~ tax_level) +
		theme(
			axis.text.x = element_text(angle = 45, hjust = 1),
			panel.grid = element_blank()
		) +
		guides(
			shape = guide_legend(override.aes = list(fill = "black")),
			color = guide_legend(override.aes = list(alpha = 1))
		)
	return(p_net_struct_og)
}

return_edge_stats_xlim_max <- function (OG) {
	if (OG == "cog") {
		edge_stats_xlim_max <- 251
	} else if (OG == "kegg") {
		edge_stats_xlim_max <- 286
	}
	return(edge_stats_xlim_max)
}

vis_edge_stats <- function(dt_net_edge_og, shared_profile, edge_stats_xlim_max) {
	p_edge <- dt_net_edge_og %>%
		mutate(n_edge = 1) %>%
		as_tibble() %>%
		ggplot(aes(x = n_tax, y = n_edge, fill = .data[[shared_profile]], color = .data[[shared_profile]])) +
		geom_bar(stat = "identity", position = "fill") +
		geom_vline(xintercept = 8, size = 0.1, color = "red") +
		facet_grid(edge_type ~ .) +
		scale_fill_manual(values = cbp8) +
		scale_color_manual(values = cbp8)
	if (OG == "cog") {
		return(p_edge + xlim(0, edge_stats_xlim_max))
	} else if (OG == "kegg") {
		return(p_edge + xlim(0, edge_stats_xlim_max))
	}
}

vis_edge_stats_category <- function(dt_net_edge_og, profile_A, profile_B, edge_stats_xlim_max) {
	df_net_edge <- dt_net_edge_og %>% select(n_tax, edge_type, .data[[profile_A]], .data[[profile_B]]) %>%
		as_tibble() %>%
		mutate(
			As = strsplit(.data[[profile_A]], ";"),
			Bs = strsplit(.data[[profile_B]], ";"),
			shared = map2(As, Bs, ~ if (is.na(.x[[1]]) | is.na(.y[[1]])) {NA} else {intersect(.x, .y)}),
			n_edge = 1
		) %>%
		select(-.data[[profile_A]], -.data[[profile_B]], -As, -Bs) %>%
		unnest(c(shared))
	df_net_edge_contra <- df_net_edge %>%
		filter(is.na(shared) == FALSE, edge_type == "contrapositive", n_tax >= 8) %>%
		group_by(shared) %>%
		summarize(count = dplyr::n()) %>%
		arrange(desc(count)) %>%
		head(6)
	df_net_edge_cooccu <- df_net_edge %>%
		filter(is.na(shared) == FALSE, edge_type == "cooccurring", n_tax >= 40) %>%
		group_by(shared) %>%
		summarize(count = dplyr::n()) %>%
		arrange(desc(count)) %>%
		head(6)
	categories <- rbind(df_net_edge_contra, df_net_edge_cooccu) %>% select(shared) %>% arrange() %>% distinct() %>% pull()
	df_net_edge$shared[!(df_net_edge$shared %in% categories)][!is.na(df_net_edge$shared[!(df_net_edge$shared %in% categories)])] <- "Others"
	df_net_edge$shared <- factor(df_net_edge$shared, levels = (c(categories, "Others")))
	p_edge <- df_net_edge %>%
		filter(!is.na(shared)) %>%
		mutate(n_edge = 1) %>%
		as_tibble() %>%
		ggplot(aes(x = n_tax, y = n_edge, fill = shared, color = shared)) +
		geom_bar(stat = "identity", position = "fill") +
		geom_vline(xintercept = 8, size = 0.1, color = "red") +
		facet_grid(edge_type ~ .) +
		scale_fill_manual(values = c(cbp12[1:length(categories)], cbp12[12])) +
		scale_color_manual(values = c(cbp12[1:length(categories)], cbp12[12]))
	if (OG == "cog") {
		return(p_edge + xlim(0, edge_stats_xlim_max))
	} else if (OG == "kegg") {
		return(p_edge + xlim(0, edge_stats_xlim_max))
	}
}

vis_edge_stats_category_kegg <- function(dt_net_edge_og, profile_A, profile_B, edge_stats_xlim_max) {
	if (profile_A == "ko_path_num_specific_A") {
		df_kegg_list <- read.csv("db/kegg/pathway/map_title.tab", header = FALSE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>%
			mutate(shared = sprintf("%05d", V1)) %>%
			rename(shared_mod = "V2")
	}
	df_net_edge <- dt_net_edge_og %>% select(n_tax, edge_type, .data[[profile_A]], .data[[profile_B]]) %>%
		as_tibble() %>%
		mutate(
			As = strsplit(.data[[profile_A]], ";"),
			Bs = strsplit(.data[[profile_B]], ";"),
			shared = map2(As, Bs, ~ if (is.na(.x[[1]]) | is.na(.y[[1]])) {NA} else {intersect(.x, .y)}),
			n_edge = 1
		) %>%
		select(-.data[[profile_A]], -.data[[profile_B]], -As, -Bs) %>%
		unnest(c(shared)) %>%
		left_join(df_kegg_list, by = "shared")
	df_net_edge_contra <- df_net_edge %>%
		filter(is.na(shared_mod) == FALSE, edge_type == "contrapositive", n_tax >= 8) %>%
		group_by(shared_mod) %>%
		summarize(count = dplyr::n()) %>%
		arrange(desc(count)) %>%
		head(6)
	df_net_edge_cooccu <- df_net_edge %>%
		filter(is.na(shared_mod) == FALSE, edge_type == "cooccurring", n_tax >= 40) %>%
		group_by(shared_mod) %>%
		summarize(count = dplyr::n()) %>%
		arrange(desc(count)) %>%
		head(5)
	categories <- rbind(df_net_edge_contra, df_net_edge_cooccu) %>% select(shared_mod) %>% arrange() %>% distinct() %>% pull()
	df_net_edge$shared_mod[!(df_net_edge$shared_mod %in% categories)][!is.na(df_net_edge$shared_mod[!(df_net_edge$shared_mod %in% categories)])] <- "Others"
	df_net_edge$shared_mod <- factor(df_net_edge$shared_mod, levels = (c(categories, "Others")))
	p_edge <- df_net_edge %>%
		filter(!is.na(shared_mod)) %>%
		mutate(n_edge = 1) %>%
		as_tibble() %>%
		ggplot(aes(x = n_tax, y = n_edge, fill = shared_mod, color = shared_mod)) +
		geom_bar(stat = "identity", position = "fill") +
		geom_vline(xintercept = 8, size = 0.1, color = "red") +
		facet_grid(edge_type ~ .) +
		scale_fill_manual(values = c(cbp12[1:length(categories)], cbp12[12])) +
		scale_color_manual(values = c(cbp12[1:length(categories)], cbp12[12]))
	if (OG == "cog") {
		return(p_edge + xlim(0, edge_stats_xlim_max))
	} else if (OG == "kegg") {
		return(p_edge + xlim(0, edge_stats_xlim_max))
	}
}

net_ratio_cooccu_contra <- function (df_net_struct_og, indicator) {
	df_ratio <- df_net_struct_og %>% select(tax, , tax_level, Domain, n_genomes, edge_type, .data[[indicator]]) %>%
		pivot_wider(names_from = edge_type, values_from = .data[[indicator]]) %>%
		mutate(ratio = cooccurring / contrapositive) %>%
		select(tax, tax_level, Domain, n_genomes, ratio)
	return(df_ratio)
}

# Ancestral reconstruction
edge_thre4og_replacement_cooccu <- 6
edge_thre4og_replacement_contra <- 4

edge_stats4og_replacement_count_cooccu <- function (df_net_edge_stats) {
	df_net_edge_stats_cooccu <- df_net_edge_stats %>% filter(edge_type == "cooccurring", n_tax >= edge_thre4og_replacement_cooccu)
	return(df_net_edge_stats_cooccu)
}

edge_stats4og_replacement_count_contra <- function (df_net_edge_stats) {
	df_net_edge_stats_contra <- df_net_edge_stats %>% filter(edge_type == "contrapositive", n_tax >= edge_thre4og_replacement_contra)
	return(df_net_edge_stats_contra)
}

count_og_replacement <- function (x) {
	og_a <- df_net_edge_stats_select[x, "og_id_A"] %>% pull()
	og_b <- df_net_edge_stats_select[x, "og_id_B"] %>% pull()
	og_a <- sub("^ko:", "", og_a)
	og_b <- sub("^ko:", "", og_b)
	if (OG == "cog" & dom == "archaea") {
		og_a_path <- paste0("data/pastml_cog2020_ar/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_cog2020_ar/", og_b, ".txt")
		tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Archaea"
	} else if (OG == "cog" & dom == "bacteria") {
		og_a_path <- paste0("data/pastml_cog2020_bac/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_cog2020_bac/", og_b, ".txt")
		tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Bacteria"
	} else if (OG == "kegg" & dom == "archaea") {
		og_a_path <- paste0("data/pastml_ko_ar/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_ko_ar/", og_b, ".txt")
		tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Archaea"
	} else if (OG == "kegg" & dom == "bacteria") {
		og_a_path <- paste0("data/pastml_ko_bac/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_ko_bac/", og_b, ".txt")
		tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Bacteria"
	}
	try_a <- try(df_og_a <- read.csv(og_a_path, header = TRUE, sep = "\t", quote = "") %>% as_tibble(), silent = TRUE)
	try_b <- try(df_og_b <- read.csv(og_b_path, header = TRUE, sep = "\t", quote = "") %>% as_tibble(), silent = TRUE)
	if (sum(class(try_a) != "try-error") > 0 & sum(class(try_b) != "try-error") > 0) {
		gtdb_node_labeled_tree <- read.tree(tree_path)
		df_genomes <- read_genomes_analyzed() %>% select(accession, Domain:Species) %>% rename(label = "accession")
		df_og_a_pivot <- df_og_a %>% mutate(presence = 1) %>% group_by(node) %>%
			pivot_wider(names_from = all_of(og_a), values_from = presence, values_fill = 0)
		if (dim(df_og_a_pivot)[2] == 3) {
			df_og_a_pivot <- df_og_a_pivot %>% rename(label = "node", og_a_TRUE = "TRUE", og_a_FALSE = "FALSE")
		} else {
			df_og_a_pivot <- df_og_a_pivot %>% rename(label = "node", og_a_TRUE = "TRUE") %>% mutate(og_a_FALSE = 0)
		}
		df_og_b_pivot <- df_og_b %>% mutate(presence = 1) %>% group_by(node) %>%
			pivot_wider(names_from = all_of(og_b), values_from = presence, values_fill = 0)
		if (dim(df_og_b_pivot)[2] == 3) {
			df_og_b_pivot <- df_og_b_pivot %>% rename(label = "node", og_b_TRUE = "TRUE", og_b_FALSE = "FALSE")
		} else {
			df_og_b_pivot <- df_og_b_pivot %>% rename(label = "node", og_b_TRUE = "TRUE") %>% mutate(og_b_FALSE = 0)
		}
		df_tree_mapped <- gtdb_node_labeled_tree %>% as_tibble() %>%
			left_join(df_og_a_pivot, by = "label") %>%
			left_join(df_og_b_pivot, by = "label") %>%
			mutate(
				og_a_presence = map2(og_a_TRUE, og_a_FALSE, ~
					if (.x == 1 & .y == 0){
						1
					} else if (.x == 1 & .y == 1) {
						0.5
					} else if (.x == 0 & .y == 1) {
						0
					}
				),
				og_b_presence = map2(og_b_TRUE, og_b_FALSE, ~
					if (.x == 1 & .y == 0){
						1
					} else if (.x == 1 & .y == 1) {
						0.5
					} else if (.x == 0 & .y == 1) {
						0
					}
				)
			) %>%
			unnest(c(og_a_presence, og_b_presence)) %>%
			mutate(
				state_child = paste0(
					as.character(og_a_presence),
					":",
					as.character(og_b_presence)
				)
			)
		node2state <- df_tree_mapped$state_child
		names(node2state) <- df_tree_mapped$node
		node2label <- df_tree_mapped$label
		names(node2label) <- df_tree_mapped$node
		df_tree_state <- df_tree_mapped %>%
			mutate(
				tree_edge_id = paste(node2label[parent], node2label[node], sep = "--"),
				state_parent = node2state[parent]
			) %>%
			group_by(state_child, state_parent) %>%
			summarize(
				n_tree_edges = dplyr::n(),
				tree_edges = paste(tree_edge_id, collapse = ";")
			) %>%
			ungroup() %>%
			mutate(
				edge_state_trans = map2(
					state_child,
					state_parent,
					~ 
					if (.x == "0:0" & .y == "0:0") {
						"n_N_N"
					} else if (.x == "0:0" & .y == "1:0") {
						"n_N_A"
					} else if (.x == "0:0" & .y == "0:1") {
						"n_N_B"
					} else if (.x == "0:0" & .y == "1:1") {
						"n_N_AB"
					} else if (.x == "0:0") {
						"n_N_U"
					} else if (.x == "1:0" & .y == "0:0") {
						"n_A_N"
					} else if (.x == "1:0" & .y == "1:0") {
						"n_A_A"
					} else if (.x == "1:0" & .y == "0:1") {
						"n_A_B"
					} else if (.x == "1:0" & .y == "1:1") {
						"n_A_AB"
					} else if (.x == "1:0") {
						"n_A_U"
					} else if (.x == "0:1" & .y == "0:0") {
						"n_B_N"
					} else if (.x == "0:1" & .y == "1:0") {
						"n_B_A"
					} else if (.x == "0:1" & .y == "0:1") {
						"n_B_B"
					} else if (.x == "0:1" & .y == "1:1") {
						"n_B_AB"
					} else if (.x == "0:1") {
						"n_B_U"
					} else if (.x == "1:1" & .y == "0:0") {
						"n_AB_N"
					} else if (.x == "1:1" & .y == "1:0") {
						"n_AB_A"
					} else if (.x == "1:1" & .y == "0:1") {
						"n_AB_B"
					} else if (.x == "1:1" & .y == "1:1") {
						"n_AB_AB"
					} else if (.x == "1:1") {
						"n_AB_U"
					}
				)
			) %>%
			unnest(c(edge_state_trans))
		tree_edges_A_B <- df_tree_state %>%
			filter(edge_state_trans == "n_A_B") %>%
			select(tree_edges) %>%
			pull()
		tree_edges_B_A <- df_tree_state %>%
			filter(edge_state_trans == "n_B_A") %>%
			select(tree_edges) %>%
			pull()
		n_A_B <- df_tree_state %>%
			filter(edge_state_trans == "n_A_B") %>%
			select(n_tree_edges) %>%
			pull()
		n_B_A <- df_tree_state %>%
			filter(edge_state_trans == "n_B_A") %>%
			select(n_tree_edges) %>%
			pull()
		if (length(tree_edges_A_B) == 0) {
			tree_edges_A_B <- NA
			n_A_B <- 0
		}
		if (length(tree_edges_B_A) == 0) {
			tree_edges_B_A <- NA
			n_B_A <- 0
		}
		n_replace <- n_A_B + n_B_A
		df_tree_state %>% group_by(edge_state_trans) %>%
			summarize(n_tree_edges = sum(n_tree_edges)) %>%
			ungroup() %>%
			pivot_wider(
				names_from = edge_state_trans,
				values_from = n_tree_edges
			) %>%
			mutate(
				og_id_A = og_a,
				og_id_B = og_b,
				tree_dom = Dom,
				n_replace = n_replace,
				tree_edges_A_B = tree_edges_A_B,
				tree_edges_B_A = tree_edges_B_A
			)
	} else {
		tibble(
			og_id_A = og_a,
			og_id_B = og_b,
			tree_dom = Dom,
			tree_edges_A_B = NA,
			tree_edges_B_A = NA
		)
	}
}

load_og_replacement <- function (OG, dom) {
	df_og_replacement <- read.table(
		paste0("data/tsv/og_displacement_", dom, "_", OG, ".tsv"),
		header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")
	) %>% as_tibble()
	return(df_og_replacement)
}

og_replacement_prokaryotes <- function (df_og_replacement_ar, df_og_replacement_bac) {
	df_og_replacement_prok <- rbind(df_og_replacement_ar, df_og_replacement_bac) %>%
		group_by(og_id_A, og_id_B) %>%
		summarize(
			n_replace_prok = sum(n_replace, na.rm = TRUE),
			n_A_B_prok = sum(n_A_B, na.rm = TRUE),
			n_B_A_prok = sum(n_B_A, na.rm = TRUE),
			n_replace_ar = first(n_replace, na.rm = TRUE),
			n_A_B_ar = first(n_A_B, na.rm = TRUE),
			n_B_A_ar = first(n_B_A, na.rm = TRUE),
			n_replace_bac = last(n_replace, na.rm = TRUE),
			n_A_B_bac = last(n_A_B, na.rm = TRUE),
			n_B_A_bac = last(n_B_A, na.rm = TRUE)
		) %>%
		ungroup() %>%
		arrange(desc(n_replace_prok))
	return(df_og_replacement_prok)
}

combine_with_net_edge_stats <- function (df_for_edge) {
	df_net_edge_og_mod <- right_join(
		as_tibble(dt_net_edge_og) %>% mutate(edge_id = paste(og_id_A, og_id_B, sep = "_")),
		df_for_edge %>% mutate(edge_id = paste(og_id_A, og_id_B, sep = "_")) %>% select(-og_id_A, -og_id_B),
		by = "edge_id"
	)
	return(df_net_edge_og_mod)
}

combine_og_replacement_with_net_edge_stats <- function (df_og_replacement_prok) {
	df_net_edge_og_replace <- right_join(
		as_tibble(dt_net_edge_og) %>% mutate(edge_id = paste(og_id_A, og_id_B, sep = "_")),
		df_og_replacement_prok %>% select(-og_id_A, -og_id_B),
		by = "edge_id"
	) %>%
	filter(
		(edge_type == "cooccurring" & n_tax >= edge_thre4og_replacement_cooccu) |
		(edge_type == "contrapositive" & n_tax >= edge_thre4og_replacement_contra)
	)
	return(df_net_edge_og_replace)
}

og_displacement_freq_tree <- function(sub_network, alpha_vertices, beta_vertices) {
	gtdb_node_labeled_tree <- read.tree(tree_path)
	df_og_displacement <- read.csv(og_displacement_path, header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>% as_tibble()
	g_og_displacement <- df_og_displacement %>% graph_from_data_frame(directed = FALSE)
	df_og_displacement_sub <- g_og_displacement %s% graph_from_edgelist(as_edgelist(sub_network), directed = FALSE) %>%
		as_data_frame() %>% as_tibble()
	df_og_displacement_sub <- df_og_displacement_sub %>% left_join(
			rbind(
				tibble(
					from = alpha_vertices,
					from_group = "alpha"
				),
				tibble(
					from = beta_vertices,
					from_group = "beta"
				)
			),
			by = "from"
		) %>% left_join(
			rbind(
				tibble(
					to = alpha_vertices,
					to_group = "alpha"
				),
				tibble(
					to = beta_vertices,
					to_group = "beta"
				)
			),
			by = "to"
		) %>%
		rowwise() %>%
		mutate(
			tree_edges_alpha_beta =
				if (from_group == "alpha" & to_group == "beta") {
					tree_edges_A_B
				} else if (from_group == "beta" & to_group == "alpha") {
					tree_edges_B_A
				},
			tree_edges_beta_alpha =
				if (from_group == "alpha" & to_group == "beta") {
					tree_edges_B_A
				} else if (from_group == "beta" & to_group == "alpha") {
					tree_edges_A_B
				},
		)
	df_displaced_edges <- full_join(
		df_og_displacement_sub %>%
		select(tree_edges_alpha_beta) %>%
		mutate(tree_edge_id = strsplit(tree_edges_alpha_beta, ";")) %>%
		unnest(tree_edge_id) %>%
		group_by(tree_edge_id) %>%
		summarize(freq_alpha_beta = dplyr::n()),
		df_og_displacement_sub %>%
		select(tree_edges_beta_alpha) %>%
		mutate(tree_edge_id = strsplit(tree_edges_beta_alpha, ";")) %>%
		unnest(tree_edge_id) %>%
		group_by(tree_edge_id) %>%
		summarize(freq_beta_alpha = dplyr::n()),
		by = "tree_edge_id"
	)
	df_displaced_edges$freq_alpha_beta[is.na(df_displaced_edges$freq_alpha_beta)] <- 0
	df_displaced_edges$freq_beta_alpha[is.na(df_displaced_edges$freq_beta_alpha)] <- 0
	df_og_displacement_tree <- gtdb_node_labeled_tree %>%
		as_tibble()
	node2label <- df_og_displacement_tree$label
	names(node2label) <- df_og_displacement_tree$node
	df_og_displacement_tree <- df_og_displacement_tree %>%
		mutate(
			parent_node_label = node2label[parent],
			tree_edge_id = paste(parent_node_label, label, sep = "--")
		) %>%
		left_join(
			df_displaced_edges %>% mutate(total_freq = freq_alpha_beta + freq_beta_alpha),
			by = "tree_edge_id"
		)
	df_og_displacement_tree[is.na(df_og_displacement_tree)] <- 0
	p_og_displacement_tree <- df_og_displacement_tree %>%
		as.treedata() %>%
		ggtree(color = "gray")
	df_genomes <- read_genomes_analyzed() %>% select(accession, Domain:Species) %>% rename(label = "accession")
	df_tax <- df_genomes[which(df_genomes$label %in% p_og_displacement_tree$data[p_og_displacement_tree$data$isTip,]$label), c("label", tax_level)]
	df_tax_top <- df_tax %>% group_by(.data[[tax_level]]) %>% summarize(count = dplyr::n()) %>% arrange(desc(count)) %>%
		head(11) %>% select(-count) %>% mutate(tax_color = .data[[tax_level]])
	df_tax <- df_tax %>% left_join(df_tax_top, by = tax_level)
	df_tax$tax_color[is.na(df_tax$tax_color)] <- "Others"
	df_tax$tax_color <- factor(df_tax$tax_color, levels = c(df_tax_top$tax_color, "Others"))
	df_tax_data <- left_join(p_og_displacement_tree$data, df_tax, by = "label") %>% filter(isTip == TRUE)
	return(list(p_og_displacement_tree, df_tax_data))
}

og_displacement_freq_tree_point <- function(tree_color) {
	p_cog_displacement_tree_out <- p_cog_displacement_tree +
		geom_point(
			data = p_cog_displacement_tree$data %>% filter(.data[[tree_color]] > 0),
			aes(x = branch, y = y, color = .data[[tree_color]]),
			size = 0.2,
		) +
		scale_color_viridis_c(na.value = "gray")
	return(p_cog_displacement_tree_out)
}

og_displacement_freq_tree_point_size <- function() {
	df_cog_displacement_tree_long <- p_cog_displacement_tree$data %>%
		pivot_longer(
			cols = c(freq_alpha_beta, freq_beta_alpha),
			names_to = "displacement_type",
			values_to = "frequency",
		) %>%
		filter(frequency > 0)
	p_cog_displacement_tree_out <- p_cog_displacement_tree +
		geom_point(
			data = df_cog_displacement_tree_long,
			aes(x = branch, y = y, size = frequency, color = displacement_type),
			shape = 16,
			alpha = 0.5,
		) +
		geom_fruit(
			data = df_tax_data %>% rename(ID = "label") %>% select(ID, tax_color) %>% mutate(x_pos = 1),
			geom = geom_star,
			mapping = aes(y = ID, x = x_pos, fill = tax_color),
			starshape = 25,
			color = NA,
			alpha = 0.8,
			pwidth = 0.08
		) +
		scale_fill_manual(values = cbp12, na.value = cbp12[12])
	p_cog_displacement_tree_plot <- df_cog_displacement_tree_long %>%
		ggplot(
			aes(x = x, y = frequency)
		) +
		geom_pointdensity() +
		scale_color_viridis_c() +
		facet_grid(displacement_type ~ .) +
		theme_bw()
	return(list(p_cog_displacement_tree_out, p_cog_displacement_tree_plot))
}

# Tree

ggtree_object <- function(sub_tax) {
	sub_tax
	if (is.na(sub_tax) == TRUE) {
		p_tree <- ggtree(gtdb_node_labeled_tree, branch.length = 'none', layout = 'circular', color = "gray")
	} else if (sub("__.*", "", sub_tax) == "p") {
		p_tree <- ggtree(
			drop.tip(gtdb_node_labeled_tree,
				df_genomes %>%
					filter(.data[["Phylum"]] != sub_tax) %>%
					select(label) %>%
					pull()
			),
			branch.length = 'none', layout = 'circular', color = "gray"
		)
	} else if (sub("__.*", "", sub_tax) == "c") {
		p_tree <- ggtree(
			drop.tip(gtdb_node_labeled_tree,
				df_genomes %>%
					filter(.data[["Class"]] != sub_tax) %>%
					select(label) %>%
					pull()
			),
			branch.length = 'none', layout = 'circular', color = "gray"
		)
	} else if (sub("__.*", "", sub_tax) == "o") {
		p_tree <- ggtree(
			drop.tip(gtdb_node_labeled_tree,
				df_genomes %>%
					filter(.data[["Order"]] != sub_tax) %>%
					select(label) %>%
					pull()
			),
			branch.length = 'none', layout = 'circular', color = "gray"
		)
	} else if (sub("__.*", "", sub_tax) == "f") {
		p_tree <- ggtree(
			drop.tip(gtdb_node_labeled_tree,
				df_genomes %>%
					filter(.data[["Family"]] != sub_tax) %>%
					select(label) %>%
					pull()
			),
			branch.length = 'none', layout = 'circular', color = "gray"
		)
	} else if (sub("__.*", "", sub_tax) == "g") {
		p_tree <- ggtree(
			drop.tip(gtdb_node_labeled_tree,
				df_genomes %>%
					filter(.data[["Genus"]] != sub_tax) %>%
					select(label) %>%
					pull()
			),
			branch.length = 'none', layout = 'circular', color = "gray"
		)
	}
	return(p_tree)
}

# OG distance

edge_thre4og_distance_cooccu <- 50
edge_thre4og_distance_contra <- 7

edge_stats4og_distance_cooccu <- function (df_net_edge_stats) {
	df_net_edge_stats_cooccu <- df_net_edge_stats %>% filter(edge_type == "cooccurring", n_tax >= edge_thre4og_distance_contra)
	return(df_net_edge_stats_cooccu)
}

edge_stats4og_distance_contra <- function (df_net_edge_stats) {
	df_net_edge_stats_contra <- df_net_edge_stats %>% filter(edge_type == "contrapositive", n_tax >= edge_thre4og_distance_cooccu)
	return(df_net_edge_stats_contra)
}

og_distance <- function (n) {
	og_a <- df_net_edge_stats_select[n, ]$og_id_A
	og_b <- df_net_edge_stats_select[n, ]$og_id_B
	if (OG == "cog") {
		path_tmp_og_a <- paste0("data/genes_cog2020/", og_a, ".tsv")
		path_tmp_og_b <- paste0("data/genes_cog2020/", og_b, ".tsv")
		df_og_a <- read.csv(path_tmp_og_a, header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>%
			filter(cog_id == og_a) %>%
			select(gene_id:prodigal_strand) %>%
			rename(start = prodigal_start, end = prodigal_end) %>%
			mutate(contig_id = sub("_[0-9]*$", "", gene_id))
		df_og_b <- read.csv(path_tmp_og_b, header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>%
			filter(cog_id == og_b) %>%
			select(gene_id:prodigal_strand) %>%
			rename(start = prodigal_start, end = prodigal_end) %>%
			mutate(contig_id = sub("_[0-9]*$", "", gene_id))
	} else if (OG == "kegg") {
		path_tmp_og_a <- paste0("data/genes_kegg/", sub("^ko:", "", og_a), ".tsv")
		path_tmp_og_b <- paste0("data/genes_kegg/", sub("^ko:", "", og_b), ".tsv")
		df_og_a <- read.csv(path_tmp_og_a, header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>%
			filter(ko_id == og_a) %>%
			select(gene_id:prodigal_strand) %>%
			rename(start = prodigal_start, end = prodigal_end) %>%
			mutate(contig_id = sub("_[0-9]*$", "", gene_id))
		df_og_b <- read.csv(path_tmp_og_b, header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>%
			filter(ko_id == og_b) %>%
			select(gene_id:prodigal_strand) %>%
			rename(start = prodigal_start, end = prodigal_end) %>%
			mutate(contig_id = sub("_[0-9]*$", "", gene_id))
	}
	rgs_a <- df_og_a %>% group_by(contig_id, accession) %>%
		nest() %>%
		rename(data_A = "data")
	rgs_b <- df_og_b %>% group_by(contig_id) %>%
		select(-accession) %>%
		nest() %>%
		rename(data_B = "data")
	rgs_ab <- inner_join(rgs_a, rgs_b, by = "contig_id")
	if (dim(rgs_ab)[1] == 0) {
		tibble(
			og_id_A = og_a,
			og_id_B = og_b,
			n_AB_genomes = 0,
			min_intergenic_mean = 0,
			min_intergenic_median = 0,
			min_intergenic_max = 0,
			min_intergenic_min = 0,
			min_intergenic_sd = 0
		)
	} else {
		rgs_ab %>%
			mutate(
				min_intergenic = map2(
					data_A,
					data_B,
					~ mcols(distanceToNearest(as_iranges(.x), as_iranges(.y)))$distance %>% min()
				)
			) %>%
			unnest(c(min_intergenic)) %>%
			group_by(accession) %>%
			summarize(min_intergenic = min(min_intergenic)) %>%
			summarize(
				og_id_A = og_a,
				og_id_B = og_b,
				n_AB_genomes = dplyr::n(),
				min_intergenic_mean = mean(min_intergenic),
				min_intergenic_median = median(min_intergenic),
				min_intergenic_max = max(min_intergenic),
				min_intergenic_min = min(min_intergenic),
				min_intergenic_sd = sd(min_intergenic)
			)
	}
}

## Number of edges

total_n_edge <- function (dt_net_edge_og) {
	dt_n_edge <- dt_net_edge_og %>%
		group_by(edge_type, n_tax) %>%
		summarize(n_edge = dplyr::n()) %>%
		ungroup()
	p_n_edge <- dt_n_edge %>% as_tibble() %>%
		ggplot(aes(x = n_tax, y = n_edge, color = edge_type)) +
		geom_point() +
		geom_vline(xintercept = 8, size = 0.1, color = "red") +
		annotation_logticks(sides = "l") +
		scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
		scale_color_manual(values = cbp12) +
		xlim(0, return_edge_stats_xlim_max(OG)) +
		labs(title = paste0("Number of edges in whole ", og_name, " network"), x = "Number of taxa", y = "Number of edges", color = "Appeared in")
	return(p_n_edge)
}

both_n_edge <- function (dt_net_edge_og) {
	dt_net_edge_both <- dt_net_edge_og %>% select(og_id_A, og_id_B, edge_type) %>%
		mutate(presence = 1) %>%
		pivot_wider(names_from = "edge_type", values_from = "presence", values_fill = 0) %>%
		mutate(n_edge_types = contrapositive + cooccurring) %>%
		filter(n_edge_types >= 2)
	p_n_edge_both <- combine_with_net_edge_stats(as_tibble(dt_net_edge_both)) %>%
		select(og_id_A:n_tax, edge_type) %>%
		pivot_wider(names_from = "edge_type", values_from = "n_tax") %>%
		ggplot(aes(x = cooccurring, y = contrapositive)) +
		geom_hex() +
		geom_vline(xintercept = 8, size = 0.1, color = "red") +
		geom_hline(yintercept = 8, size = 0.1, color = "red") +
		scale_fill_viridis(trans = "log") +
		labs(title = paste0("Number of taxa of edges\nshared in both types of networks"), subtitle =  paste0(og_name, " network"), x = "Cooccurring", y = "Contrapositive")
	return(p_n_edge_both)
}

## Codon usage

codon_distance <- function(n) {
	gid <- df_proteins_codon_count[n,]$gene_id
	acc <- df_proteins_codon_count[n,]$accession
	gene_codon <- df_proteins_codon_ratio_vecs[,gid] %>% pull()
	genome_codon <- df_genome_codon_ratio_vecs[,acc] %>% pull()
	tibble(
		gene_id = gid,
		accession = acc,
		codon_distance = sqrt(sum((gene_codon - genome_codon)^2))
	)
}

codon_distance_old <- function (n) {
	og_a <- df_net_edge_stats_select[n, ]$og_id_A
	og_b <- df_net_edge_stats_select[n, ]$og_id_B
	if (OG == "cog") {
		path_tmp_og_a <- paste0("data/genes_cog2020/", og_a, ".tsv")
		path_tmp_og_b <- paste0("data/genes_cog2020/", og_b, ".tsv")
		df_og_a <- read.csv(path_tmp_og_a, header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>%
			filter(cog_id == og_a) %>%
			select(gene_id:accession)
		df_og_b <- read.csv(path_tmp_og_b, header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>%
			filter(cog_id == og_b) %>%
			select(gene_id:accession)
	} else if (OG == "kegg") {
		path_tmp_og_a <- paste0("data/genes_kegg/", sub("^ko:", "", og_a), ".tsv")
		path_tmp_og_b <- paste0("data/genes_kegg/", sub("^ko:", "", og_b), ".tsv")
		df_og_a <- read.csv(path_tmp_og_a, header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>%
			filter(cog_id == og_a) %>%
			select(gene_id:accession)
		df_og_b <- read.csv(path_tmp_og_b, header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>%
			filter(cog_id == og_b) %>%
			select(gene_id:accession)
	}
	df_distance_a <- df_og_a %>%
	#dt_proteins_codon_count %>% select(gene_id,accession) %>%
	#	head(1000) %>%
	#	as_tibble() %>%
		mutate(
			distance = map2(gene_id, accession, ~ sqrt(sum((df_proteins_codon_ratio_vecs[,.x] - df_genome_codon_ratio_vecs[,.y])^2)))
		) %>%
		unnest(c(distance))
	df_distance_b <- df_og_b %>%
	#dt_proteins_codon_count %>% select(gene_id,accession) %>%
	#	head(1000) %>%
	#	as_tibble() %>%
		mutate(
			distance = map2(gene_id, accession, ~ sqrt(sum((df_proteins_codon_ratio_vecs[,.x] - df_genome_codon_ratio_vecs[,.y])^2)))
		) %>%
		unnest(c(distance))
	if (dim(rgs_ab)[1] == 0) {
		tibble(
			og_id_A = og_a,
			og_id_B = og_b,
			n_AB_genomes = 0,
			min_intergenic_mean = 0,
			min_intergenic_max = 0,
			min_intergenic_min = 0,
			min_intergenic_sd = 0
		)
	} else {
		rgs_ab %>%
			mutate(
				min_intergenic = map2(
					data_A,
					data_B,
					~ mcols(distanceToNearest(as_iranges(.x), as_iranges(.y)))$distance %>% min()
				)
			) %>%
			unnest(c(min_intergenic)) %>%
			group_by(accession) %>%
			summarize(min_intergenic = min(min_intergenic)) %>%
			summarize(
				og_id_A = og_a,
				og_id_B = og_b,
				n_AB_genomes = dplyr::n(),
				min_intergenic_mean = mean(min_intergenic),
				min_intergenic_max = max(min_intergenic),
				min_intergenic_min = min(min_intergenic),
				min_intergenic_sd = sd(min_intergenic)
			)
	}
}

# Experiments

## Make growth curve from TECAN plate reader output xlsx files

growth_curves <- function(id, xlsx, max_od, max_time, out) {
	# import files
	id <- sub(".xlsx", ".tsv", xlsx)
	df_id <- read.table(id, header = TRUE, sep = "\t", quote = "")
	csv <- sub("/data/", "/csv/", sub(".xlsx", ".csv", xlsx))
	system(paste0("xlsx2csv -a ", xlsx, " | perl scripts/perl/xlsx2csv.pl > ", csv))
	df_csv <- read.csv(csv, header = TRUE, quote = "")
	# time calculation and pivot_longer
	start <- sort(strptime(paste(df_csv$date, df_csv$time), "%m-%d-%y %H:%M:%S"))[1]
	df_csv_l <- df_csv %>% pivot_longer(c(-alphabet, -date, -time, -mode, -wave, -band, -flash, -temperature), names_to = "number", values_to = "OD") %>%
		mutate(well = sub("X", "", paste0(alphabet, number)), date_time = strptime(paste(date, time), "%m-%d-%y %H:%M:%S"), start_time = start) %>% mutate(hours = as.numeric(date_time - start_time, units = "hours"))
	# make table for plot, colored by category
	df_plot <- inner_join(df_csv_l, df_id, by = "well")
	df_names <- df_id %>% select(name, category) %>% distinct()
	color_n <- df_names %>% select(category) %>% distinct() %>% nrow()
	df_names <- df_names %>% select(category) %>% distinct() %>% mutate(color = viridis(color_n)) %>% right_join(df_names, by = "category")
	df_plot$name <- factor(df_plot$name, levels = df_names$name)
	return(list(df_plot, df_names))
}

growth_curves_v2 <- function(id, xlsx, max_od, max_time, out) {
	# import files
	csv <- sub("/data/", "/csv/", sub(".xlsx", ".csv", xlsx))
	system(paste0("xlsx2csv -a ", xlsx, " | perl scripts/perl/xlsx2csv_v2.pl ", csv, " ", id))
	df_csv <- read.csv(csv, header = TRUE, quote = "")
	df_id <- read.table(id, header = TRUE, sep = "\t", quote = "")
	# time calculation and pivot_longer
	start <- sort(strptime(paste(df_csv$date, df_csv$time), "%m-%d-%y %H:%M:%S"))[1]
	df_csv_l <- df_csv %>% pivot_longer(c(-alphabet, -date, -time), names_to = "number", values_to = "OD") %>%
		mutate(well = sub("X", "", paste0(alphabet, number)), date_time = strptime(paste(date, time), "%m-%d-%y %H:%M:%S"), start_time = start) %>% mutate(hours = as.numeric(date_time - start_time, units = "hours"))
	# make table for plot, colored by category
	df_plot <- inner_join(df_csv_l, df_id, by = "well")
	df_plot$name <- factor(df_plot$name, levels = df_plot$name %>% unique())
	df_plot$strain <- factor(df_plot$strain, levels = df_plot$strain %>% unique())
	df_plot$media <- factor(df_plot$media, levels = df_plot$media %>% unique())
	df_plot$category <- factor(df_plot$category, levels = df_plot$category %>% unique())
	return(df_plot)
}

save_growth_curves_v2 <- function(df_plot, id, xlsx, max_od, max_time, out, height, width, color_by) {
	# Set colorblind-friendly palette
	cbp <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
	if (color_by == "name") {
		# plot with coloring by names
		p <- ggplot(df_plot, aes(x = hours, y = OD, colour = name)) + geom_line() + geom_point() +
			xlim(0, max_time) + ylim(0, max_od) + scale_colour_manual(values = cbp) + theme_bw()
		ggsave(out, p, height = height, width = width)
	} else if (color_by == "strain") {
		# plot with coloring by strains
		p <- ggplot(df_plot, aes(x = hours, y = OD, colour = strain)) + geom_line() + geom_point() +
			xlim(0, max_time) + ylim(0, max_od) + scale_colour_manual(values = cbp) + theme_bw()
		ggsave(out, p, height = height, width = width)
	} else if (color_by == "media") {
		# plot with coloring by media
		p <- ggplot(df_plot, aes(x = hours, y = OD, colour = media)) + geom_line() + geom_point() +
			xlim(0, max_time) + ylim(0, max_od) + scale_colour_manual(values = cbp) + theme_bw()
		ggsave(out, p, height = height, width = width)
	} else if (color_by == "category") {
		# plot with coloring by categories
		p <- ggplot(df_plot, aes(x = hours, y = OD, colour = category)) + geom_line() + geom_point() +
			xlim(0, max_time) + ylim(0, max_od) + scale_colour_manual(values = cbp) + theme_bw()
		ggsave(out, p, height = height, width = width)
	}
}

save_growth_curves_facet_v2 <- function(df_plot, id, xlsx, max_od, max_time, out, height, width, color_by, facet_by) {
	# Set colorblind-friendly palette
	cbp <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
	if (color_by == "strain" && facet_by == "media") {
		# plot with coloring by strains
		p <- ggplot(df_plot, aes(x = hours, y = OD, colour = strain)) + geom_line() + geom_point() +
			xlim(0, max_time) + ylim(0, max_od) + scale_colour_manual(values = cbp) + theme_bw() + facet_wrap(~ media)
		ggsave(out, p, height = height, width = width)
	} else if (color_by == "media" && facet_by == "strain") {
		# plot with coloring by media
		p <- ggplot(df_plot, aes(x = hours, y = OD, colour = media)) + geom_line() + geom_point() +
			xlim(0, max_time) + ylim(0, max_od) + scale_colour_manual(values = cbp) + theme_bw() + facet_wrap(~ strain)
		ggsave(out, p, height = height, width = width)
	}
}

growth_curves_v3 <- function(xlsx) {
	# import files
	id <- sub(".xlsx", ".tsv", xlsx)
	df_id <- read.table(id, header = TRUE, sep = "\t", quote = "")
	csv <- sub("/data/", "/csv/", sub(".xlsx", ".csv", xlsx))
	system(paste0("xlsx2csv -a ", xlsx, " | perl scripts/perl/xlsx2csv_v3.pl ", csv, " ", id))
	df_csv <- data.frame()
	csv_path <- sub("/csv/.*", "/csv/", csv)
	csv_file <- sub(".*/csv/", "", csv)
	for (csv_each in list.files(csv_path, pattern=paste0(csv_file, "."))){
		df_csv_each <- read.csv(paste0(csv_path, csv_each), header = TRUE, quote = "", stringsAsFactors=F)
		if (df_csv_each$alphabet[1] == FALSE) {
			df_csv_each$alphabet <- "F"
		}
		df_csv_each_l <- df_csv_each %>% pivot_longer(c(-alphabet, -date, -time), names_to = "number", values_to = "OD")
		df_csv <- rbind(df_csv, df_csv_each_l)
	}
	# time calculation and pivot_longer
	start <- sort(strptime(paste(df_csv$date, df_csv$time), "%m-%d-%y %H:%M:%S"))[1]
	df_csv_l <- df_csv %>% mutate(well = sub("X", "", paste0(alphabet, number)), date_time = strptime(paste(date, time),
		"%m-%d-%y %H:%M:%S"), start_time = start) %>% mutate(hours = as.numeric(date_time - start_time, units = "hours"))
	# make table for plot, colored by category
	df_plot <- inner_join(df_csv_l, df_id, by = "well")
	df_plot$name <- factor(df_plot$name, levels = df_plot$name %>% unique())
	df_plot$strain <- factor(df_plot$strain, levels = df_plot$strain %>% unique())
	df_plot$media <- factor(df_plot$media, levels = df_plot$media %>% unique())
	df_plot$category <- factor(df_plot$category, levels = df_plot$category %>% unique())
	return(df_plot)
}

draw_growth_curves <- function(id, xlsx, max_od, max_time, out) {
	plot_name_dfs <- growth_curves(id, xlsx, max_od, max_time, out)
	df_plot <- plot_name_dfs[[1]]
	df_names <- plot_name_dfs[[2]]
	# plot
	p <- ggplot(df_plot, aes(x = hours, y = OD, colour = name)) + geom_line() + geom_point() +
		xlim(0, max_time) + ylim(0, max_od) + scale_color_manual(values = df_names$color)
	return(p)
}

save_growth_curves <- function(id, xlsx, max_od, max_time, out, height, width) {
	plot_name_dfs <- growth_curves(id, xlsx, max_od, max_time, out)
	df_plot <- plot_name_dfs[[1]]
	df_names <- plot_name_dfs[[2]]
	# plot
	p <- ggplot(df_plot, aes(x = hours, y = OD, colour = name)) + geom_line() + geom_point() +
		xlim(0, max_time) + ylim(0, max_od) + scale_color_manual(values = df_names$color)
	ggsave(out, p, height = height, width = width)
}

draw_growth_curve_category <- function(id, xlsx, max_od, max_time, out) {
	plot_name_dfs <- growth_curves(id, xlsx, max_od, max_time, out)
	df_plot <- plot_name_dfs[[1]]
	df_names <- plot_name_dfs[[2]]
	# summarize data
	df_summarize <- df_plot %>% group_by(category, hours) %>% summarize(mean_OD = mean(OD), sd_OD = sd(OD)) %>% dplyr::rename(OD = mean_OD, sd = sd_OD)
	# plot
	p <- ggplot(df_summarize, aes(x = hours, y = OD, colour = category)) + geom_line() + geom_point() + geom_errorbar(aes(ymin = OD - sd, ymax = OD + sd, width = 0.1)) +
		xlim(0, max_time) + ylim(0, max_od) + scale_color_viridis(discrete=TRUE)
	return(p)
}

save_growth_curve_category <- function(id, xlsx, max_od, max_time, out, height, width) {
	plot_name_dfs <- growth_curves(id, xlsx, max_od, max_time, out)
	df_plot <- plot_name_dfs[[1]]
	df_names <- plot_name_dfs[[2]]
	# summarize data
	df_summarize <- df_plot %>% group_by(category, hours) %>% summarize(mean_OD = mean(OD), sd_OD = sd(OD)) %>% dplyr::rename(OD = mean_OD, sd = sd_OD)
	# plot
	p <- ggplot(df_summarize, aes(x = hours, y = OD, colour = category)) + geom_line() + geom_point() + geom_errorbar(aes(ymin = OD - sd, ymax = OD + sd, width = 0.1)) +
		xlim(0, max_time) + ylim(0, max_od) + scale_color_viridis(discrete=TRUE)
	ggsave(out, p, height = height, width = width)
}

# Make plots with different colors and types

growth_curves_2param <- function(id, xlsx, max_od, max_time, out) {
	# import files
	df_id <- read.table(id, header = TRUE, sep = "\t", quote = "")
	csv <- sub("/data/", "/csv/", sub(".xlsx", ".csv", xlsx))
	system(paste0("xlsx2csv -a ", xlsx, " | perl scripts/perl/xlsx2csv.pl > ", csv))
	df_csv <- read.csv(csv, header = TRUE, quote = "")
	# time calculation and pivot_longer
	start <- sort(strptime(paste(df_csv$date, df_csv$time), "%m-%d-%y %H:%M:%S"))[1]
	df_csv_l <- df_csv %>% pivot_longer(c(-alphabet, -date, -time, -mode, -wave, -band, -flash, -temperature), names_to = "number", values_to = "OD") %>%
		mutate(well = sub("X", "", paste0(alphabet, number)), date_time = strptime(paste(date, time), "%m-%d-%y %H:%M:%S"), start_time = start) %>% mutate(hours = as.numeric(date_time - start_time, units = "hours"))
	# make table for plot, colored by category
	df_plot <- inner_join(df_csv_l, df_id, by = "well")
	df_names <- df_id %>% select(name, color) %>% distinct()
	color_n <- df_names %>% select(color) %>% distinct() %>% nrow()
	df_names <- df_names %>% select(color) %>% distinct() %>% mutate(color = viridis(color_n)) %>% right_join(df_names, by = "color")
	df_plot$name <- factor(df_plot$name, levels = df_names$name)
	return(list(df_plot, df_names))
}

save_growth_curve_category_2param <- function(id, xlsx, max_od, max_time, out, height, width) {
	plot_name_dfs <- growth_curves_2param(id, xlsx, max_od, max_time, out)
	df_plot <- plot_name_dfs[[1]]
	df_names <- plot_name_dfs[[2]]
	# summarize data
	df_summarize <- df_plot %>% group_by(color, type, hours) %>% summarize(mean_OD = mean(OD), sd_OD = sd(OD)) %>% dplyr::rename(OD = mean_OD, sd = sd_OD)
	# plot
	p <- ggplot(df_summarize, aes(x = hours, y = OD, colour = color, linetype = type)) + geom_line() + geom_point(aes(shape = type)) + geom_errorbar(aes(ymin = OD - sd, ymax = OD + sd, width = 0.1)) +
		xlim(0, max_time) + ylim(0, max_od) + scale_color_viridis(discrete=TRUE) + theme_bw()
	ggsave(out, p, height = height, width = width)
}

# Habitat analysis (GMGC)

gmgc_habitat_stats <- function (ogid, download_date) {
	df_gmgc <- read.csv(paste0("data/gmgc/gmgc_", ogid, "_", download_date, ".csv"), header = TRUE, sep = ",", quote = "\"", na.strings=c("", "NULL")) %>% as_tibble()
	df_gmgc$Habitat[df_gmgc$Habitat == "-"] <- NA
	n_unigene <- df_gmgc %>% distinct(Unigene) %>% nrow()
	n_source <- df_gmgc %>% distinct(Source) %>% nrow()
	n_unigene_progenomes <- df_gmgc %>% filter(Source == "ProGenomes v2") %>% nrow()
	n_unigene_with_habitat <- df_gmgc %>% filter(!is.na(Habitat)) %>% nrow()
	n_source_with_habitat <- df_gmgc %>% filter(!is.na(Habitat)) %>% distinct(Source) %>% nrow()
	df_gmgc_long <- df_gmgc %>%
		filter(!is.na(Habitat)) %>%
		mutate(
			habitat = strsplit(Habitat, ", ")
		) %>%
		unnest(c(habitat))
	n_habitat <- df_gmgc_long %>% distinct(habitat) %>% nrow()
	df_gmgc_stats <- tibble(
		`OG ID` = sub("cog", "COG", ogid),
		`Number of Unigenes` = n_unigene,
		`Number of sources` = n_source,
		`Number of Unigenes from Progenomes` = n_unigene_progenomes,
		`Number of Unigenes with habitats` = n_unigene_with_habitat,
		`Number of sources with habitats` = n_source_with_habitat,
		`Number of habitats` = n_habitat,
	)
	df_gmgc_habitats <- df_gmgc_long %>%
		group_by(habitat) %>%
		summarize(n_unigene = dplyr::n()) %>%
		mutate(
			bar_color = "Others"
		) %>%
		rbind(
			tibble(
				habitat = "Total Unigenes",
				n_unigene = n_unigene_with_habitat,
				bar_color = "Total Unigenes"
			)
		) %>%
		arrange(desc(n_unigene))
	df_gmgc_habitats$habitat <- factor(df_gmgc_habitats$habitat, levels = rev(df_gmgc_habitats$habitat))
	df_gmgc_habitats$bar_color[grep("gut", df_gmgc_habitats$habitat)] <- "Gut"
	df_gmgc_habitats$bar_color[grep("marine", df_gmgc_habitats$habitat)] <- "Marine"
	df_gmgc_habitats$bar_color[grep("soil", df_gmgc_habitats$habitat)] <- "Soil"
	df_gmgc_habitats$bar_color[grep("freshwater", df_gmgc_habitats$habitat)] <- "Freshwater"
	df_gmgc_habitats$bar_color <- factor(df_gmgc_habitats$bar_color, levels = c("Total Unigenes", "Gut", "Marine", "Freshwater", "Soil", "Others"))
	p_gmgc_habitats <- df_gmgc_habitats %>%
		ggplot(aes(x = n_unigene, y = habitat, fill = bar_color)) +
		geom_bar(stat = "identity") +
		scale_fill_manual(values = c("black", cbp8[6], cbp8[5], cbp8[2], cbp8[1], "grey50")) +
		labs(
			title = "GMGC habitat",
			subtitle = sub("cog", "COG", ogid), 
			x = "Number of Unigenes",
			y = "Habitats"
		)
	return(list(df_gmgc_stats, p_gmgc_habitats))
}

# Protein domain

gene_tree_hmmsearch <- function (treefile, hmmfile, seg_size) {
	# Load hmmsearch result
	df_hmm <- read.csv(hmmfile, header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>% as_tibble()

	# Load tree
	tree_og <- read.tree(treefile)
	# Remove double quotation from GS output
	tree_og$tip.label <- gsub("\"", "", tree_og$tip.label)

	# Phyum color
	og <- sub("data/gs_tree/", "", sub(".nwk", "", treefile))
	og_tsv <- paste0("data/tsv/proteins_raw_gtdbtax_", og, "_uniq.tsv.gz")
	df_og_seqs <- read.csv(og_tsv, header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
		as_tibble()
	df_og_seqs$aa_len <- AAStringSet(df_og_seqs$protein_seq) %>% width()
	top11_phyla <- df_og_seqs %>%
		group_by(Phylum) %>%
		summarize(n_seqs = dplyr::n()) %>%
		ungroup() %>%
		arrange(desc(n_seqs)) %>%
		head(11) %>%
		select(Phylum) %>%
		pull()
	df_tax <- tibble(Phylum = c(top11_phyla), tax_color = c(top11_phyla))
	df_tax$tax_color <- factor(df_tax$tax_color, levels = c(top11_phyla, "Others"))
	df_og_seqs_phylum <- df_og_seqs %>%
		left_join(df_tax, by = "Phylum")
	df_og_seqs_phylum$tax_color[is.na(df_og_seqs_phylum$tax_color)] <- "Others"
	df_heatmap <- rbind(
	  df_og_seqs_phylum %>% select(gene_id,tax_color) %>% mutate(gene_pos = 1),
	  df_og_seqs_phylum %>% select(gene_id,tax_color) %>% mutate(gene_pos = 2)
	)

	# Draw tree
	p_tree <- ggtree(tree_og, color = "gray", size = 0.5) +
	  geom_fruit(
	    data = df_heatmap,
	    geom = geom_tile,
	    mapping = aes(y = gene_id, x = gene_pos, fill = tax_color),
	    color = NA,
	    pwidth = 0.15,
	    offset = 0.15
	  ) +
	  scale_fill_manual(values = cbp12)

	# Add domain
	p2 <- facet_plot(p_tree, panel = "Protein domain", data = df_og_seqs_phylum,
	           geom = geom_spoke, mapping = aes(x = 1, angle = 0, radius = aa_len), color = "gray80", size = seg_size)
	p3 <- facet_plot(p2, panel = "Protein domain", data = df_hmm,
	             geom = geom_spoke, mapping = aes(x = from_ali, angle = 0, radius = to_ali - from_ali, color = dom_name), size = seg_size) +
	  scale_color_brewer(palette = "Dark2")
	
	return(p3)
}

gene_tree_hmmsearch_limit <- function (treefile, hmmfile, seg_size, max_aa_len) {
	# Load hmmsearch result
	df_hmm <- read.csv(hmmfile, header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
	  as_tibble() %>% filter(from_ali <= max_aa_len)

	# Load tree
	tree_og <- read.tree(treefile)
	# Remove double quotation from GS output
	tree_og$tip.label <- gsub("\"", "", tree_og$tip.label)

	# Phyum color
	og <- sub("data/gs_tree/", "", sub(".nwk", "", treefile))
	og_tsv <- paste0("data/tsv/proteins_raw_gtdbtax_", og, "_uniq.tsv.gz")
	df_og_seqs <- read.csv(og_tsv, header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
		as_tibble()
	df_og_seqs$aa_len <- AAStringSet(df_og_seqs$protein_seq) %>% width()
	df_og_seqs$aa_len[df_og_seqs$aa_len > max_aa_len] <- max_aa_len
	top11_phyla <- df_og_seqs %>%
		group_by(Phylum) %>%
		summarize(n_seqs = dplyr::n()) %>%
		ungroup() %>%
		arrange(desc(n_seqs)) %>%
		head(11) %>%
		select(Phylum) %>%
		pull()
	df_tax <- tibble(Phylum = c(top11_phyla), tax_color = c(top11_phyla))
	df_tax$tax_color <- factor(df_tax$tax_color, levels = c(top11_phyla, "Others"))
	df_og_seqs_phylum <- df_og_seqs %>%
		left_join(df_tax, by = "Phylum")
	df_og_seqs_phylum$tax_color[is.na(df_og_seqs_phylum$tax_color)] <- "Others"
	df_heatmap <- rbind(
	  df_og_seqs_phylum %>% select(gene_id,tax_color) %>% mutate(gene_pos = 1),
	  df_og_seqs_phylum %>% select(gene_id,tax_color) %>% mutate(gene_pos = 2)
	)

	# Draw tree
	p_tree <- ggtree(tree_og, color = "gray", size = 0.5) +
	  geom_fruit(
	    data = df_heatmap,
	    geom = geom_tile,
	    mapping = aes(y = gene_id, x = gene_pos, fill = tax_color),
	    color = NA,
	    pwidth = 0.15,
	    offset = 0.15
	  ) +
	  scale_fill_manual(values = cbp12)

	# Add domain
	p2 <- facet_plot(p_tree, panel = "Protein domain", data = df_og_seqs_phylum,
	           geom = geom_spoke, mapping = aes(x = 1, angle = 0, radius = aa_len), color = "gray80", size = seg_size)
	p3 <- facet_plot(p2, panel = "Protein domain", data = df_hmm,
	             geom = geom_spoke, mapping = aes(x = from_ali, angle = 0, radius = to_ali - from_ali, color = dom_name), size = seg_size) +
	  scale_color_brewer(palette = "Dark2")
	
	return(p3)
}


gene_tree_hmmsearch_score_bin <- function (treefile, hmmfile, score_bin, metatsv, len_max, len_break) {
  # Load hmmsearch result
  df_hmm <- read.csv(hmmfile, header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
    as_tibble()
  # Load metadata
  df_meta <- read.csv(metatsv, header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
    as_tibble()
  # Load tree
  tree_og <- read.tree(treefile)
  # Remove double quotation from GS output
  tree_og$tip.label <- gsub("\"", "", tree_og$tip.label)
  
  # Phyum color
  phylum_colors <- c(
    "p__Proteobacteria",
    "p__Actinobacteriota",
    "p__Bacteroidota",
    "p__Firmicutes",
    "p__Firmicutes_A",
    "p__Patescibacteria",
    "p__Cyanobacteria",
    "p__Campylobacterota",
    "p__Desulfobacterota",
    "p__Planctomycetota",
    "p__Verrucomicrobiota",
    "Other_Bacteria",
    "d__Archaea"
  )
  df_meta2 <- df_meta %>% select(gene_id,accession,prodigal_start,prodigal_end,Domain,Phylum) %>% 
    mutate(
      seq_len = (prodigal_end - prodigal_start + 1) / 3,
      seq_len_mod = if_else(seq_len >= len_max, len_max, seq_len),
      phylum_color = map2(Domain, Phylum, ~ if_else(.x == "d__Bacteria", .y, "d__Archaea"))
    ) %>% 
    unnest(c(phylum_color))
  df_meta2$phylum_color <- factor(df_meta2$phylum_color, levels = phylum_colors)
  df_meta2$phylum_color[is.na(df_meta2$phylum_color)] <- "Other_Bacteria"
  ## phylum heatmap tsv
  df_heatmap_phylum <- rbind(
    df_meta2 %>% select(gene_id,phylum_color) %>% mutate(gene_pos = 1),
    df_meta2 %>% select(gene_id,phylum_color) %>% mutate(gene_pos = 2)
  )
  
  # Protein domain presence
  df_dom <- df_hmm %>% group_by(target) %>% summarize(score_max = max(score_dom)) %>% mutate(score_bin = cut(score_max, score_bin))
  ## protein domain heatmap tsv
  df_heatmap_dom <- rbind(
    df_dom %>% select(target,score_bin) %>% mutate(gene_pos = 1),
    df_dom %>% select(target,score_bin) %>% mutate(gene_pos = 2)
  )
  
  # Draw tree
  p_tree <- ggtree(tree_og, color = "gray", size = 0.5) +
    geom_fruit(
      data = df_heatmap_phylum,
      geom = geom_tile,
      mapping = aes(y = gene_id, x = gene_pos, fill = phylum_color),
      color = NA,
      pwidth = 0.26,
      offset = 0.18
    ) +
    scale_fill_manual(values = c(cbp12,"black"), drop = FALSE) +
    new_scale_fill() +
    geom_fruit(
      data = df_heatmap_dom,
      geom = geom_tile,
      mapping = aes(y = target, x = gene_pos, fill = score_bin),
      color = NA,
      pwidth = 0.26,
      offset = 0.3
    ) +
    geom_fruit(
      data = df_meta2,
      geom = geom_bar,
      mapping = aes(y = gene_id, x = seq_len_mod),
      orientation = "y", 
      stat = "identity",
      #color = NA,
      pwidth = 0.55,
      offset = 0.15,
      axis.params=list(
        axis       = "x",
        text.size  = 1,
        text.angle = 90,
        hjust      = 1.1,
        vjust      = 0.5,
        nbreak     = len_break,
        line.size = 0.05
      ),
      grid.params=list()
    ) #+
    #ylim(-length(tree_og$tip.label)/40, length(tree_og$tip.label))
  
  return(p_tree)
}

species_tree_heatmap_manual_1 <- function (dom, ogs, og_colors, p_width, text_size) {
  #dom <- "archaea"
  #ogs <- "MgtE_TM_CBS_N;MgtE_TM_CBS_N_gt590;MgtE_TM_CBS_N_N;MgtE_TM_CBS;MgtE_TM_CBS_gt380;MgtE_TM_TM;MgtE_TM;MgtE_TM_gt280;MgtE_TM_lt170;MgtE_TM_N;MgtE_TM_N_N"
  #og_colors <- "MgtE;MgtE;MgtE;MgtE_dN;MgtE_dN;MgtE_dCBS_dN;MgtE_dCBS_dN;MgtE_dCBS_dN;MgtE_dCBS_dN;MgtE_dCBS;MgtE_dCBS"
  if (dom == "archaea") {
    og_path_suffix <- "ar"
    tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Archaea"
  } else if (dom == "bacteria") {
    og_path_suffix <- "bac"
    tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Bacteria"
  }
  # Load tree
  gtdb_node_labeled_tree <- read.tree(tree_path)
  # Draw tree
  p_tree <- ggtree(gtdb_node_labeled_tree, color = "gray", size = 0.5)
  # Load presence/absence matrix
  OGs <- strsplit(ogs, ";") %>% unlist
  OG_colors <- strsplit(og_colors, ";") %>% unlist()
  df_og_matrix <- tibble(accession = gtdb_node_labeled_tree$tip.label)
  for (i in OGs) {
    #i <- "MgtE_TM_CBS_N"
    if (length(grep("^COG", i)) == 1) {
      df_og <- read.csv(paste0("data/pastml_csv_cog2020_", og_path_suffix, "/", i, ".csv"), header = TRUE, sep = ",", quote = "") %>% 
        as_tibble()
      df_og_matrix <- left_join(df_og_matrix, df_og, by = "accession")
    } else if (length(grep("^K", i)) == 1) {
      df_og <- read.csv(paste0("data/pastml_csv_ko_", og_path_suffix, "/", i, ".csv"), header = TRUE, sep = ",", quote = "") %>% 
        as_tibble()
      df_og_matrix <- left_join(df_og_matrix, df_og, by = "accession")
    } else {
      df_og <- read.csv(paste0("data/pastml_csv_manual_", og_path_suffix, "/", i, ".csv"), header = TRUE, sep = ",", quote = "") %>% 
        as_tibble()
      df_og_matrix <- left_join(df_og_matrix, df_og, by = "accession")
    }
  }
  df_presence <- df_og_matrix %>% pivot_longer(cols = c(-accession), names_to = "col_pos", values_to = "presence") %>% 
    left_join(tibble(col_pos = OGs, tl_color = OG_colors), by = "col_pos") %>% 
    mutate(Tile_color = map2(presence, tl_color, ~ if (.x == FALSE) {NA} else {.y})) %>% unnest(Tile_color) %>% 
    select(accession,Tile_color,col_pos) %>% dplyr::rename(label = "accession")
  # Load genome & taxa info
  df_genomes <- read_genomes_analyzed() %>% select(accession, Domain:Species) %>% rename(label = "accession")
  # Phylum color
  tax_level <- "Phylum"
  df_tax <- df_genomes[which(df_genomes$label %in% p_tree$data[p_tree$data$isTip,]$label), c("label", tax_level)]
  df_tax_top <- df_tax %>% group_by(.data[[tax_level]]) %>% summarize(count = dplyr::n()) %>% arrange(desc(count)) %>%
    head(11) %>% select(-count) %>% mutate(Tile_color = .data[[tax_level]])
  df_tax <- df_tax %>% left_join(df_tax_top, by = tax_level) %>% mutate(col_pos = tax_level)
  df_tax$Tile_color[is.na(df_tax$Tile_color)] <- "Others"
  # Make tile color data frame
  df_heatmap <- df_tax %>% select(-all_of(tax_level)) %>% rbind(df_presence)
  df_heatmap$Tile_color <- factor(df_heatmap$Tile_color, levels = c(df_tax_top$Tile_color, "Others", unique(OG_colors)))
  df_heatmap$col_pos <- factor(df_heatmap$col_pos, levels = c(tax_level, OGs))
  # Draw heatmap
  p2 <- p_tree + 
    geom_fruit(
      data = df_heatmap, geom = geom_tile, mapping = aes(y = label, x = col_pos, fill = Tile_color), color = NA,
      axis.params = list(axis = "x", text.size = 3, text.angle = 90, hjust = 1), pwidth = p_width
    ) +
    scale_fill_manual(
      values = c(cbp12[1:length(df_tax_top$Tile_color)], cbp12[12], palette.colors(palette = "Dark 2")), 
      drop = FALSE, na.value = NA
    ) +
    ylim(-length(gtdb_node_labeled_tree$tip.label)/10, length(gtdb_node_labeled_tree$tip.label)) +
    theme(text = element_text(family = "Arial", size = text_size))
  return(p2)
}

species_tree_heatmap_manual_2 <- function (dom, ogs, p_width, text_size) {
  #dom <- "archaea"
  #ogs <- "MgtE_TM_CBS_N;MgtE_TM_CBS_N_gt590;MgtE_TM_CBS_N_N;MgtE_TM_CBS;MgtE_TM_CBS_gt380;MgtE_TM_TM;MgtE_TM;MgtE_TM_gt280;MgtE_TM_lt170;MgtE_TM_N;MgtE_TM_N_N"
  #og_colors <- "MgtE;MgtE;MgtE;MgtE_dN;MgtE_dN;MgtE_dCBS_dN;MgtE_dCBS_dN;MgtE_dCBS_dN;MgtE_dCBS_dN;MgtE_dCBS;MgtE_dCBS"
  if (dom == "archaea") {
    og_path_suffix <- "ar"
    tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Archaea"
  } else if (dom == "bacteria") {
    og_path_suffix <- "bac"
    tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Bacteria"
  }
  # Load tree
  gtdb_node_labeled_tree <- read.tree(tree_path)
  # Draw tree
  p_tree <- ggtree(gtdb_node_labeled_tree, color = "gray", size = 0.5)
  # Load presence/absence matrix
  OGs <- strsplit(ogs, ";") %>% unlist
  #OG_colors <- strsplit(og_colors, ";") %>% unlist()
  df_og_matrix <- tibble(accession = gtdb_node_labeled_tree$tip.label)
  for (i in OGs) {
    #i <- "MgtE_TM_CBS_N"
    if (length(grep("^COG", i)) == 1) {
      df_og <- read.csv(paste0("data/pastml_csv_cog2020_", og_path_suffix, "/", i, ".csv"), header = TRUE, sep = ",", quote = "") %>% 
        as_tibble()
      df_og_matrix <- left_join(df_og_matrix, df_og, by = "accession")
    } else if (length(grep("^K", i)) == 1) {
      df_og <- read.csv(paste0("data/pastml_csv_ko_", og_path_suffix, "/", i, ".csv"), header = TRUE, sep = ",", quote = "") %>% 
        as_tibble()
      df_og_matrix <- left_join(df_og_matrix, df_og, by = "accession")
    } else {
      df_og <- read.csv(paste0("data/pastml_csv_copies_", og_path_suffix, "/", i, ".csv"), header = TRUE, sep = ",", quote = "") %>% 
        as_tibble()
      df_og_matrix <- left_join(df_og_matrix, df_og, by = "accession")
    }
  }
  df_presence <- df_og_matrix %>% pivot_longer(cols = c(-accession), names_to = "col_pos", values_to = "presence") %>% 
    #left_join(tibble(col_pos = OGs, tl_color = OG_colors), by = "col_pos") %>% 
    mutate(Tile_color = map(presence, ~ if (.x == 0) {NA} else if (.x >= 8) {8} else {.x})) %>% unnest(Tile_color) %>% 
    select(accession,Tile_color,col_pos) %>% dplyr::rename(label = "accession")
  # Load genome & taxa info
  df_genomes <- read_genomes_analyzed() %>% select(accession, Domain:Species) %>% rename(label = "accession")
  # Phylum color
  tax_level <- "Phylum"
  df_tax <- df_genomes[which(df_genomes$label %in% p_tree$data[p_tree$data$isTip,]$label), c("label", tax_level)]
  df_tax_top <- df_tax %>% group_by(.data[[tax_level]]) %>% summarize(count = dplyr::n()) %>% arrange(desc(count)) %>%
    head(11) %>% select(-count) %>% mutate(Tile_color = .data[[tax_level]])
  df_tax <- df_tax %>% left_join(df_tax_top, by = tax_level) %>% mutate(col_pos = tax_level)
  df_tax$Tile_color[is.na(df_tax$Tile_color)] <- "Others"
  # Make tile color data frame
  df_heatmap <- df_tax %>% select(-all_of(tax_level)) %>% rbind(df_presence)
  df_heatmap$Tile_color <- factor(df_heatmap$Tile_color, levels = c(df_tax_top$Tile_color, "Others", as.character(1:8)))
  df_heatmap$col_pos <- factor(df_heatmap$col_pos, levels = c(tax_level, OGs))
  # Draw heatmap
  p2 <- p_tree + 
    geom_fruit(
      data = df_heatmap, geom = geom_tile, mapping = aes(y = label, x = col_pos, fill = Tile_color), color = NA,
      axis.params = list(axis = "x", text.size = 3, text.angle = 90, hjust = 1), pwidth = p_width
    ) +
    scale_fill_manual(
      values = c(cbp12[1:length(df_tax_top$Tile_color)], cbp12[12], palette.colors(palette = "Dark 2")), 
      drop = FALSE, na.value = NA
    ) +
    ylim(-length(gtdb_node_labeled_tree$tip.label)/10, length(gtdb_node_labeled_tree$tip.label)) +
    theme(text = element_text(family = "Arial", size = text_size))
  return(p2)
}

species_tree_with_taxa <- function (dom, p_width, off_set) {
  #dom <- "archaea"
  if (dom == "archaea") {
    og_path_suffix <- "ar"
    tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Archaea"
  } else if (dom == "bacteria") {
    og_path_suffix <- "bac"
    tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Bacteria"
  }
  # Load tree
  gtdb_node_labeled_tree <- read.tree(tree_path)
  # Draw tree
  p_tree <- ggtree(gtdb_node_labeled_tree, color = "gray", size = 0.5)
  # Load genome & taxa info
  df_genomes <- read_genomes_analyzed() %>% select(accession, Domain:Species) %>% rename(label = "accession")
  # Phylum color
  tax_level <- "Phylum"
  df_tax <- df_genomes[which(df_genomes$label %in% p_tree$data[p_tree$data$isTip,]$label), c("label", tax_level)]
  df_tax_top <- df_tax %>% group_by(.data[[tax_level]]) %>% summarize(count = dplyr::n()) %>% arrange(desc(count)) %>%
    head(11) %>% select(-count) %>% mutate(Tile_color = .data[[tax_level]])
  df_tax <- df_tax %>% left_join(df_tax_top, by = tax_level) %>% mutate(col_pos = tax_level)
  df_tax$Tile_color[is.na(df_tax$Tile_color)] <- "Others"
  # Make tile color data frame
  df_heatmap <- df_tax %>% select(-all_of(tax_level)) #%>% rbind(df_presence)
  df_heatmap$Tile_color <- factor(df_heatmap$Tile_color, levels = c(df_tax_top$Tile_color, "Others"))
  df_heatmap$col_pos <- factor(df_heatmap$col_pos, levels = c(tax_level))
  # Draw heatmap (taxon)
  p2 <- p_tree + 
    geom_fruit(
      data = df_heatmap, geom = geom_tile, mapping = aes(y = label, x = col_pos, fill = Tile_color), color = NA,
      axis.params = list(axis = "x", text.size = 3, text.angle = 90, hjust = 1), pwidth = p_width, offset = off_set
    ) +
    scale_fill_manual(
      values = c(cbp12[1:length(df_tax_top$Tile_color)], cbp12[12], palette.colors(palette = "Dark 2")), 
      drop = FALSE, na.value = NA
    )
  return(p2)
}

species_tree_with_subtaxa <- function (dom, tax_x, tax_level, p_width, off_set) {
  #dom <- "archaea"
  if (dom == "archaea") {
    og_path_suffix <- "ar"
    tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Archaea"
  } else if (dom == "bacteria") {
    og_path_suffix <- "bac"
    tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Bacteria"
  }
  tax_prefix <- sub("__.*", "", tax_x)
  if (tax_prefix == "d") {
    subtax <- "Domain"
  } else if (tax_prefix == "p") {
    subtax <- "Phylum"
  } else if (tax_prefix == "c") {
    subtax <- "Class"
  } else if (tax_prefix == "o") {
    subtax <- "Order"
  } else if (tax_prefix == "f") {
    subtax <- "Family"
  } else if (tax_prefix == "g") {
    subtax <- "Genus"
  } else if (tax_prefix == "s") {
    subtax <- "Species"
  }
  # Load tree
  gtdb_node_labeled_tree <- read.tree(tree_path)
  # Load genome & taxa info
  df_genomes <- read_genomes_analyzed() %>% select(accession, Domain:Species) %>% rename(label = "accession")
  # Draw tree
  to_drop <- df_genomes %>% filter(.data[[subtax]] != tax_x) %>% select(label) %>% pull()
  p_tree <- ggtree(drop.tip(gtdb_node_labeled_tree, to_drop), color = "gray", size = 0.5)
  # Phylum color
  df_tax <- df_genomes[which(df_genomes$label %in% p_tree$data[p_tree$data$isTip,]$label), c("label", tax_level)]
  df_tax_top <- df_tax %>% group_by(.data[[tax_level]]) %>% summarize(count = dplyr::n()) %>% arrange(desc(count)) %>%
    head(11) %>% select(-count) %>% mutate(Tile_color = .data[[tax_level]])
  df_tax <- df_tax %>% left_join(df_tax_top, by = tax_level) %>% mutate(col_pos = tax_level)
  df_tax$Tile_color[is.na(df_tax$Tile_color)] <- "Others"
  # Make tile color data frame
  df_heatmap <- df_tax %>% select(-all_of(tax_level)) #%>% rbind(df_presence)
  df_heatmap$Tile_color <- factor(df_heatmap$Tile_color, levels = c(df_tax_top$Tile_color, "Others"))
  df_heatmap$col_pos <- factor(df_heatmap$col_pos, levels = c(tax_level))
  # Draw heatmap (taxon)
  p2 <- p_tree + 
    geom_fruit(
      data = df_heatmap, geom = geom_tile, mapping = aes(y = label, x = col_pos, fill = Tile_color), color = NA,
      axis.params = list(axis = "x", text.size = 3, text.angle = 90, hjust = 1), pwidth = p_width, offset = off_set
    ) +
    scale_fill_manual(
      values = c(cbp12[1:length(df_tax_top$Tile_color)], cbp12[12], palette.colors(palette = "Dark 2")), 
      drop = FALSE, na.value = NA
    )
  return(p2)
}

species_tree_heatmap_manual_prokatlas_1 <- function (dom, ogs, og_colors, p_width, text_size, off_set) {
  #dom <- "archaea"
  #ogs <- "MgtE_TM_CBS_N;MgtE_TM_CBS_N_gt590;MgtE_TM_CBS_N_N;MgtE_TM_CBS;MgtE_TM_CBS_gt380;MgtE_TM_TM;MgtE_TM;MgtE_TM_gt280;MgtE_TM_lt170;MgtE_TM_N;MgtE_TM_N_N"
  #og_colors <- "MgtE;MgtE;MgtE;MgtE_dN;MgtE_dN;MgtE_dCBS_dN;MgtE_dCBS_dN;MgtE_dCBS_dN;MgtE_dCBS_dN;MgtE_dCBS;MgtE_dCBS"
  if (dom == "archaea") {
    og_path_suffix <- "ar"
    tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Archaea"
  } else if (dom == "bacteria") {
    og_path_suffix <- "bac"
    tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Bacteria"
  }
  # Load tree
  gtdb_node_labeled_tree <- read.tree(tree_path)
  # Draw tree
  p_tree <- ggtree(gtdb_node_labeled_tree, color = "gray", size = 0.5)
  # Load presence/absence matrix
  OGs <- strsplit(ogs, ";") %>% unlist
  OG_colors <- strsplit(og_colors, ";") %>% unlist()
  df_og_matrix <- tibble(accession = gtdb_node_labeled_tree$tip.label)
  for (i in OGs) {
    #i <- "MgtE_TM_CBS_N"
    if (length(grep("^COG", i)) == 1) {
      df_og <- read.csv(paste0("data/pastml_csv_cog2020_", og_path_suffix, "/", i, ".csv"), check.names = FALSE, header = TRUE, sep = ",", quote = "") %>% 
        as_tibble()
      df_og_matrix <- left_join(df_og_matrix, df_og, by = "accession")
    } else if (length(grep("^K", i)) == 1) {
      df_og <- read.csv(paste0("data/pastml_csv_ko_", og_path_suffix, "/", i, ".csv"), check.names = FALSE, header = TRUE, sep = ",", quote = "") %>% 
        as_tibble()
      df_og_matrix <- left_join(df_og_matrix, df_og, by = "accession")
    } else {
      df_og <- read.csv(paste0("data/pastml_csv_manual_", og_path_suffix, "/", i, ".csv"), check.names = FALSE, header = TRUE, sep = ",", quote = "") %>% 
        as_tibble()
      df_og_matrix <- left_join(df_og_matrix, df_og, by = "accession")
    }
  }
  df_presence <- df_og_matrix %>% pivot_longer(cols = c(-accession), names_to = "col_pos", values_to = "presence") %>% 
    left_join(tibble(col_pos = OGs, tl_color = OG_colors), by = "col_pos") %>% 
    mutate(Tile_color = map2(presence, tl_color, ~ if (.x == FALSE) {NA} else {.y})) %>% unnest(Tile_color) %>% 
    select(accession,Tile_color,col_pos) %>% dplyr::rename(label = "accession")
  # Load genome & taxa info
  df_genomes <- read_genomes_analyzed() %>% select(accession, Domain:Species) %>% rename(label = "accession")
  # Phylum color
  tax_level <- "Phylum"
  df_tax <- df_genomes[which(df_genomes$label %in% p_tree$data[p_tree$data$isTip,]$label), c("label", tax_level)]
  df_tax_top <- df_tax %>% group_by(.data[[tax_level]]) %>% summarize(count = dplyr::n()) %>% arrange(desc(count)) %>%
    head(11) %>% select(-count) %>% mutate(Tile_color = .data[[tax_level]])
  df_tax <- df_tax %>% left_join(df_tax_top, by = tax_level) %>% mutate(col_pos = tax_level)
  df_tax$Tile_color[is.na(df_tax$Tile_color)] <- "Others"
  # ProkAtlas habitat preference
  df_habitat_pref_whole <- read.csv("data/tsv/gtdb_ssu_prokatlas_habitat_pref.tsv", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>% as_tibble()
  df_habitat_pref <- left_join(
    tibble(accession = gtdb_node_labeled_tree$tip.label), 
    df_habitat_pref_whole, by = "accession"
  )
  df_habitat_pref[is.na(df_habitat_pref)] <- 0
  df_envs <- df_habitat_pref %>% 
    pivot_longer(cols = c(-accession), names_to = "envs", values_to = "hp") %>% 
    group_by(envs) %>% summarize(hp_sum = sum(hp)) %>% arrange(desc(hp_sum))
  df_env_color <- df_envs %>% mutate(Habitats = map2(
      hp_sum, envs, ~ if (.x < df_envs$hp_sum[18]) {"Other_envs"} else {.y}
    )) %>% unnest(Habitats)
  df_env_color$Habitats <- factor(df_env_color$Habitats, levels = unique(df_env_color$Habitats))
  df_habitat_barplot <- df_habitat_pref %>% 
    pivot_longer(cols = c(-accession), names_to = "envs", values_to = "hp") %>% 
    left_join(df_env_color, by = "envs") %>% dplyr::rename(label = "accession")
  # Make tile color data frame
  df_heatmap <- df_tax %>% select(-all_of(tax_level)) %>% rbind(df_presence)
  df_heatmap$Tile_color <- factor(df_heatmap$Tile_color, levels = c(df_tax_top$Tile_color, "Others", unique(OG_colors)))
  df_heatmap$col_pos <- factor(df_heatmap$col_pos, levels = c(tax_level, OGs))
  # Draw heatmap (taxon, OGs)
  p2 <- p_tree + 
    geom_fruit(
      data = df_heatmap, geom = geom_tile, mapping = aes(y = label, x = col_pos, fill = Tile_color), color = NA,
      axis.params = list(axis = "x", text.size = 3, text.angle = 90, hjust = 1), pwidth = p_width
    ) +
    scale_fill_manual(
      values = c(cbp12[1:length(df_tax_top$Tile_color)], cbp12[12], palette.colors(palette = "Dark 2")), 
      drop = FALSE, na.value = NA
    )
  # Draw barplot (envs)
  p3 <- p2 + new_scale_fill() +
    geom_fruit(
      data = df_habitat_barplot, geom = geom_bar, mapping = aes(y = label, x = hp, fill = Habitats), stat = "identity",
      orientation = "y", offset = off_set
    ) +
    scale_fill_manual(values = c(pals::polychrome(18), "grey10") %>% as.vector()) +
    ylim(-length(gtdb_node_labeled_tree$tip.label)/10, length(gtdb_node_labeled_tree$tip.label)) +
    theme(text = element_text(family = "Arial", size = text_size))
  return(p3)
}


species_tree_heatmap_manual_prokatlas_2 <- function (dom, ogs, og_colors, p_width, text_size, off_set) {
  #dom <- "archaea"
  #ogs <- "MgtE_TM_CBS_N;MgtE_TM_CBS_N_gt590;MgtE_TM_CBS_N_N;MgtE_TM_CBS;MgtE_TM_CBS_gt380;MgtE_TM_TM;MgtE_TM;MgtE_TM_gt280;MgtE_TM_lt170;MgtE_TM_N;MgtE_TM_N_N"
  #og_colors <- "MgtE;MgtE;MgtE;MgtE_dN;MgtE_dN;MgtE_dCBS_dN;MgtE_dCBS_dN;MgtE_dCBS_dN;MgtE_dCBS_dN;MgtE_dCBS;MgtE_dCBS"
  if (dom == "archaea") {
    og_path_suffix <- "ar"
    tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Archaea"
  } else if (dom == "bacteria") {
    og_path_suffix <- "bac"
    tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Bacteria"
  }
  # Load tree
  gtdb_node_labeled_tree <- read.tree(tree_path)
  # Draw tree
  p_tree <- ggtree(gtdb_node_labeled_tree, color = "gray", size = 0.5)
  # Load presence/absence matrix
  OGs <- strsplit(ogs, ";") %>% unlist
  OG_colors <- strsplit(og_colors, ";") %>% unlist()
  df_og_matrix <- tibble(accession = gtdb_node_labeled_tree$tip.label)
  for (i in OGs) {
    #i <- "MgtE_TM_CBS_N"
    if (length(grep("^COG", i)) == 1) {
      df_og <- read.csv(paste0("data/pastml_csv_cog2020_", og_path_suffix, "/", i, ".csv"), check.names = FALSE, header = TRUE, sep = ",", quote = "") %>% 
        as_tibble()
      df_og_matrix <- left_join(df_og_matrix, df_og, by = "accession")
    } else if (length(grep("^K", i)) == 1) {
      df_og <- read.csv(paste0("data/pastml_csv_ko_", og_path_suffix, "/", i, ".csv"), check.names = FALSE, header = TRUE, sep = ",", quote = "") %>% 
        as_tibble()
      df_og_matrix <- left_join(df_og_matrix, df_og, by = "accession")
    } else {
      df_og <- read.csv(paste0("data/pastml_csv_manual_", og_path_suffix, "/", i, ".csv"), check.names = FALSE, header = TRUE, sep = ",", quote = "") %>% 
        as_tibble()
      df_og_matrix <- left_join(df_og_matrix, df_og, by = "accession")
    }
  }
  df_presence <- df_og_matrix %>% pivot_longer(cols = c(-accession), names_to = "col_pos", values_to = "presence") %>% 
    left_join(tibble(col_pos = OGs, tl_color = OG_colors), by = "col_pos") %>% 
    mutate(Tile_color = map2(presence, tl_color, ~ if (.x == FALSE) {NA} else {.y})) %>% unnest(Tile_color) %>% 
    select(accession,Tile_color,col_pos) %>% dplyr::rename(label = "accession")
  # Load genome & taxa info
  df_genomes <- read_genomes_analyzed() %>% select(accession, Domain:Species) %>% rename(label = "accession")
  # Phylum color
  tax_level <- "Phylum"
  df_tax <- df_genomes[which(df_genomes$label %in% p_tree$data[p_tree$data$isTip,]$label), c("label", tax_level)]
  df_tax_top <- df_tax %>% group_by(.data[[tax_level]]) %>% summarize(count = dplyr::n()) %>% arrange(desc(count)) %>%
    head(11) %>% select(-count) %>% mutate(Tile_color = .data[[tax_level]])
  df_tax <- df_tax %>% left_join(df_tax_top, by = tax_level) %>% mutate(col_pos = tax_level)
  df_tax$Tile_color[is.na(df_tax$Tile_color)] <- "Others"
  # ProkAtlas habitat preference
  df_habitat_pref_whole <- read.csv("data/tsv/gtdb_ssu_prokatlas_habitat_pref.tsv", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>% as_tibble()
  df_habitat_pref <- left_join(
    tibble(accession = gtdb_node_labeled_tree$tip.label), 
    df_habitat_pref_whole, by = "accession"
  )
  df_habitat_pref[is.na(df_habitat_pref)] <- 0
  df_envs <- df_habitat_pref %>% 
    pivot_longer(cols = c(-accession), names_to = "envs", values_to = "hp") %>% 
    group_by(envs) %>% summarize(hp_sum = sum(hp)) %>% arrange(desc(hp_sum))
  df_env_color <- df_envs %>% mutate(Habitats = map2(
    hp_sum, envs, ~ if (.x < df_envs$hp_sum[18]) {"Other_envs"} else {.y}
  )) %>% unnest(Habitats)
  df_env_color$Habitats <- factor(df_env_color$Habitats, levels = unique(df_env_color$Habitats))
  df_habitat_barplot <- df_habitat_pref %>% 
    pivot_longer(cols = c(-accession), names_to = "envs", values_to = "hp") %>% 
    left_join(df_env_color, by = "envs") %>% dplyr::rename(label = "accession")
  # Make tile color data frame
  df_heatmap <- df_tax %>% select(-all_of(tax_level)) %>% rbind(df_presence)
  df_heatmap$Tile_color <- factor(df_heatmap$Tile_color, levels = c(df_tax_top$Tile_color, "Others", unique(OG_colors)))
  df_heatmap$col_pos <- factor(df_heatmap$col_pos, levels = c(tax_level, OGs))
  # Draw heatmap (taxon, OGs)
  p2 <- p_tree + 
    geom_fruit(
      data = df_heatmap, geom = geom_tile, mapping = aes(y = label, x = col_pos, fill = Tile_color), color = NA,
      axis.params = list(axis = "x", text.size = 3, text.angle = 90, hjust = 1), pwidth = p_width, offset = off_set
    ) +
    scale_fill_manual(
      values = c(cbp12[1:length(df_tax_top$Tile_color)], cbp12[12], hcl.colors(2,"Blue-Red 2"), palette.colors(palette = "Dark 2")), 
      drop = FALSE, na.value = NA
    )
  # Draw barplot (envs)
  p3 <- p2 + #new_scale_fill() +
    #geom_fruit(
    #  data = df_habitat_barplot, geom = geom_bar, mapping = aes(y = label, x = hp, fill = Habitats), stat = "identity",
    #  orientation = "y", offset = off_set
    #) +
    #scale_fill_manual(values = c(pals::polychrome(18), "grey10") %>% as.vector()) +
    ylim(-length(gtdb_node_labeled_tree$tip.label)/10, length(gtdb_node_labeled_tree$tip.label)) +
    theme(text = element_text(family = "Arial", size = text_size))
  return(p3)
}


ogs_prokatlas_sum_hp <- function (og, dom) {
  #og <- "MgtE_TM_CBS_N"
  # Load presence/absence matrix
  if (dom == "archaea") {
    og_path_suffix <- "ar"
    Dom <- "Archaea"
    color_num <- 1:18
  } else if (dom == "bacteria") {
    og_path_suffix <- "bac"
    Dom <- "Bacteria"
    color_num <- 19:36
  }
  if (length(grep("^COG", og)) == 1) {
    df_og <- read.csv(paste0("data/pastml_csv_cog2020_", og_path_suffix, "/", og, ".csv"), header = TRUE, sep = ",", quote = "") %>% 
      as_tibble() %>% mutate(Domain = Dom, og_id = og)
  } else if (length(grep("^K", og)) == 1) {
    df_og <- read.csv(paste0("data/pastml_csv_ko_", og_path_suffix, "/", og, ".csv"), header = TRUE, sep = ",", quote = "") %>% 
      as_tibble() %>% mutate(Domain = Dom, og_id = og)
  } else {
    df_og <- read.csv(paste0("data/pastml_csv_manual_", og_path_suffix, "/", og, ".csv"), header = TRUE, sep = ",", quote = "") %>% 
      as_tibble() %>% mutate(Domain = Dom, og_id = og)
  }
  # ProkAtlas habitat preference
  df_habitat_pref_whole <- read.csv("data/tsv/gtdb_ssu_prokatlas_habitat_pref.tsv", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>% as_tibble()
  df_og_hp <- left_join(
    df_og, 
    df_habitat_pref_whole, by = "accession"
  )
  df_og_hp[is.na(df_og_hp)] <- 0
  df_og_habitat <- df_og_hp %>% 
    pivot_longer(cols = c(-accession:-og_id), names_to = "envs", values_to = "hp") %>% filter(hp > 0, .data[[og]] == TRUE) %>% 
    group_by(og_id,envs,Domain) %>% summarize(sum_hp = sum(hp)) %>% ungroup()
  # Bar plot
  df_og_habitat <- df_og_habitat %>% filter(Domain == Dom) %>% arrange(desc(sum_hp))
  df_habitat_color <- df_og_habitat %>% 
    mutate(Habitat = map2(sum_hp, envs, ~ if (.x < df_og_habitat$sum_hp[18]) {"Others"} else {.y})) %>% 
    unnest(Habitat)
  df_habitat_color$Habitat <- factor(df_habitat_color$Habitat, levels = unique(df_habitat_color$Habitat))
  p_og_habitat <- left_join(
      df_og_habitat %>% filter(Domain == Dom),
      df_habitat_color %>% select(envs,Habitat), by = "envs"
    ) %>% 
    ggplot(aes(x = og_id, y = sum_hp, fill = Habitat)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = c(pals::polychrome(36)[color_num], "grey10") %>% as.vector()) +
    facet_grid(. ~ Domain)
  return(p_og_habitat)
}


ogs_prokatlas_sum_hp_multi <- function (ogs, dom) {
  #ogs <- c("MgtE_TM_CBS_N","MgtE_TM_CBS")
  #dom <- "archaea"
  # Load presence/absence matrix
  if (dom == "archaea") {
    og_path_suffix <- "ar"
    Dom <- "Archaea"
    color_num <- 1:18
  } else if (dom == "bacteria") {
    og_path_suffix <- "bac"
    Dom <- "Bacteria"
    color_num <- 19:36
  }
  df_og <- tibble()
  for (og in ogs) {
    if (length(grep("^COG", og)) == 1) {
      df_og <- rbind(df_og, read.csv(paste0("data/pastml_csv_cog2020_", og_path_suffix, "/", og, ".csv"), header = TRUE, sep = ",", quote = "") %>% 
        as_tibble() %>% mutate(Domain = Dom, og_id = og) %>% filter(.data[[og]] == TRUE) %>% select(-all_of(og)))
    } else if (length(grep("^K", og)) == 1) {
      df_og <- rbind(df_og, read.csv(paste0("data/pastml_csv_ko_", og_path_suffix, "/", og, ".csv"), header = TRUE, sep = ",", quote = "") %>% 
        as_tibble() %>% mutate(Domain = Dom, og_id = og) %>% filter(.data[[og]] == TRUE) %>% select(-all_of(og)))
    } else {
      df_og <- rbind(df_og, read.csv(paste0("data/pastml_csv_manual_", og_path_suffix, "/", og, ".csv"), header = TRUE, sep = ",", quote = "") %>% 
        as_tibble() %>% mutate(Domain = Dom, og_id = og) %>% filter(.data[[og]] == TRUE) %>% select(-all_of(og)))
    }
  }
  # ProkAtlas habitat preference
  df_n_genoms <- df_og %>% group_by(og_id) %>% summarize(n_genomes = dplyr::n())
  df_habitat_pref_whole <- read.csv("data/tsv/gtdb_ssu_prokatlas_habitat_pref.tsv", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>% as_tibble()
  df_og_hp <- left_join(
    df_og, 
    df_habitat_pref_whole, by = "accession"
  )
  df_og_hp[is.na(df_og_hp)] <- 0
  df_og_habitat <- df_og_hp %>% 
    pivot_longer(cols = c(-accession:-og_id), names_to = "envs", values_to = "hp") %>% filter(hp > 0) %>% 
    group_by(og_id,envs,Domain) %>% summarize(sum_hp = sum(hp)) %>% ungroup() %>% 
    left_join(df_n_genoms, by = "og_id") %>% mutate(norm_sum_hp = sum_hp / n_genomes)
  # Bar plot
  df_habitat <- df_og_habitat %>% group_by(envs) %>% summarize(total_norm_sum_hp = sum(norm_sum_hp)) %>% 
    arrange(desc(total_norm_sum_hp))
  df_habitat_color <- df_habitat %>% 
    mutate(Habitat = map2(total_norm_sum_hp, envs, ~ if (.x < df_habitat$total_norm_sum_hp[18]) {"Others"} else {.y})) %>% 
    unnest(Habitat)
  df_habitat_color$Habitat <- factor(df_habitat_color$Habitat, levels = unique(df_habitat_color$Habitat))
  p_og_habitat <- left_join(
      df_og_habitat,
      df_habitat_color %>% select(envs,Habitat), by = "envs"
    ) %>% 
    ggplot(aes(x = og_id, y = sum_hp, fill = Habitat)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = c(pals::polychrome(36)[color_num], "grey10") %>% as.vector()) +
    facet_grid(. ~ Domain)
  return(p_og_habitat)
}


interproscan_range_num <- function (i) {
  rng <- df_interproscan_nest[i,]$data[[1]] %>% as_iranges()
  rng_reduced <- reduce_ranges(rng) %>% mutate(range_num = seq(1, length(start)))
  bind_cols(
    df_interproscan_nest[i,c("prot_acc", "interpro_acc")],
    join_overlap_left(rng, rng_reduced) %>% as.data.frame() %>% as_tibble()
  ) %>% ungroup()
}


interpro_acc_upset <- function(df_interpro, interpro_ids) {
  p_interpro <- df_interpro %>% filter(interpro_acc %in% interpro_ids) %>% group_by(prodigal_acc) %>% 
    distinct(interpro_acc,interpro_desc) %>% ungroup() %>% 
    mutate(interpro_name = paste(interpro_desc, interpro_acc,sep = "; ")) %>% 
    group_by(prodigal_acc) %>% summarize(interpro_names = list(interpro_name)) %>% 
    ggplot(aes(x = interpro_names)) + geom_bar() + scale_x_upset()
  return(p_interpro)
}


interpro_acc_upset_2 <- function(df_interpro, interpro_ids) {
  p_interpro <- df_interpro %>% filter(interpro_acc %in% interpro_ids) %>% group_by(prodigal_acc) %>% 
    distinct(interpro_acc, short_name,interpro_desc) %>% ungroup() %>% 
    mutate(interpro_name = paste(interpro_desc, short_name, interpro_acc,sep = "; ")) %>% 
    group_by(prodigal_acc) %>% summarize(interpro_names = list(interpro_name)) %>% 
    ggplot(aes(x = interpro_names)) + geom_bar() + scale_x_upset()
  return(p_interpro)
}

interpro_acc_upset_3 <- function(df_interpro, interpro_ids) {
  p_interpro <- df_interpro %>% filter(interpro_acc %in% interpro_ids) %>% group_by(prodigal_acc) %>% 
    distinct(interpro_acc, short_name,interpro_desc) %>% ungroup() %>% 
    mutate(interpro_name = paste(short_name, interpro_acc,sep = "; ")) %>% 
    group_by(prodigal_acc) %>% summarize(interpro_names = list(interpro_name)) %>% 
    ggplot(aes(x = interpro_names)) + geom_bar() + scale_x_upset()
  return(p_interpro)
}

amino_acid_pal <- function(name) {
  
  # aliview order
  # levels = c("K","R","H","D","E","S","T","N","Q","A","V","L","I","P","M","F","W","Y","G","C","X")
  
  aliview <- c(
    # Basic
    c("#E63319","#F6442C"),
    # Basic (Imidazole)
    "#19B3B3",
    # Acidic
    c("#CC4DCC","#984097"),
    # Neutral
    c("#029602","#45C945","#19CC19","#5CED5C"),
    # Hydrophobic
    c("#276EB7","#047DF9","#78A6D5","#4EA0F3"),
    # Hydrophobic (Imino)
    "#CDCD00",
    # Aromatic
    c("#0F549B","#1980E6","#0355A9","#14C6C8"),
    # Glycine
    "#E6994D",
    # Cysteine
    "#E68080",
    # X
    "grey70"
  )
  names(aliview) <- c(
    # Basic
    c("K","R"),
    # Basic (Imidazole)
    "H",
    # Acidic
    c("D","E"),
    # Neutral
    c("S","T","N","Q"),
    # Hydrophobic
    c("A","V","L","I"),
    # Hydrophobic (Imino)
    "P",
    # Aromatic
    c("M","F","W","Y"),
    # Glycine
    "G",
    # Cysteine
    "C",
    # X
    "X"
  )
  
    clustal <- c(rep("orange", 4),
                 rep("red", 3),
                 rep("blue", 3),
                 rep("green", 4))
    names(clustal) <- c("G", "P", "S", "T",
                        "H", "K", "R",
                        "F", "W", "Y",
                        "I", "L", "M", "V")

    lesk <- c(rep("orange", 4),
              rep("green", 10),
              rep("magenta", 3),
              rep("red", 2),
              rep("blue", 2))
    names(lesk) <- c("G", "A", "S", "T",
                     "C", "V", "I", "L", "P", "Y", "F", "Y", "M", "W",
                     "N", "Q", "H",
                     "D", "E",
                     "K", "R")

#    maeditor <- c(rep("lightgreen", 4),
#                  rep("green", 9),
#                  rep("darkgreen", 3),
#                  rep("blue", 2),
#                  rep("lilac", 2),
#                  "darkblue",
#                  rep("orange", 2),
#                  "pink",
#                  rep("red", 2))
#    names(maeditor) <- c("A", "G",
#                         "C",
#                         "S", "T",
#                         "D", "E", "N", "Q",
#                         "I", "L", "M", "V",
#                         "F", "W", "Y",
#                         "H",
#                         "K", "R",
#                         "P",
#                         "S", "T")

    cinema <- c(rep("blue", 3),
                rep("red", 2),
                rep("green", 4),
                rep("white", 5),
                rep("magenta", 3),
                rep("brown", 2),
                "yellow")
    names(cinema) <- c("H", "K", "R",
                       "D", "E",
                       "S", "T", "N", "Q",
                       "A", "V", "L", "I", "M",
                       "F", "W", "Y",
                       "P", "G",
                       "C")

#    shapely <- c(rep("E60A0A", 2),
#                 rep("E6E600", 2),
#                 rep("145AFF", 2),
#                 rep("FA9600", 2),
#                 rep("3232AA", 2),
#                 "00DCDC",
#                 rep("0F820F", 3),
#                 "C8C8C8",
#                 "B45AB4",
#                 "8282D2",
#                 "DC9682")
#    names(shapely) <- c("D", "E",
#                        "C", "M",
#                        "K", "R",
#                        "S", "T",
#                        "F", "Y",
#                        "N", "Q",
#                        "G",
#                        "L", "V", "I",
#                        "A",
#                        "W",
#                        "H",
#                        "P")

    heliquest <- c("gray", "yellow",
                   "red", "red",
                   "yellow", "gray", "lightblue", "yellow", "blue",
                   "yellow", "yellow",
                   "pink", "green", "pink",
                   "blue", "purple", "purple",
                   "yellow", "yellow", "yellow")
    names(heliquest) <- c("A", "C",
                          "D", "E",
                          "F", "G", "H", "I", "K",
                          "L", "M",
                          "N", "P", "Q",
                          "R", "S", "T",
                          "V", "W", "Y")

    switch(name,
           "aliview" = aliview,
           "clustal" = clustal,
           "lesk" = lesk,
           "cinema" = cinema,
           #"maeditor" = maeditor,
           "heliquest" = heliquest)#,
           #shapely)
}

