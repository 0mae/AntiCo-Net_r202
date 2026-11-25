# Load libraries

library(RCy3)
library(dplyr)
library(tibble)
library(tidyr)
library(igraph)
library(purrr)
library(viridis)
library(viridisLite)
library(grDevices)
library(vegan)
library(philentropy)


# Set threshold

edge_thre4og_replacement_cooccu <- 6
edge_thre4og_replacement_contra <- 4
edge_thre4og_distance_cooccu <- 50
edge_thre4og_distance_contra <- 7
edge_thre_contra <- 8
edge_thre_cooccu <- 80
replace_thre <- 0

# Read

read_cog2020_def_fixed_funcat_text <- function () {
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
		unnest(cog_fun_cat_text)
	return(df_cog2020_def)
}

read_ko_def_na <- function () {
	df_ko_def <- read.csv("db/kegg_fixed/kegg_ko_def.tsv", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
		as_tibble()
	return(df_ko_def)
}

read_network_edge_stats_general <- function(OG, edge_thre, edgetype, dt_network) {
	if (OG == "cog") {
		df_network <- dt_network %>%
			as_tibble() %>%
			select(-cog_fun_cat_A:-n_cog_fun_cats_B) %>%
			filter(n_tax >= edge_thre, edge_type == edgetype)
	} else if (OG == "kegg") {
		df_network <- dt_network %>%
			as_tibble() %>%
			select(-ko_symbol_A:-n_ko_cogs_B, -n_ko_ecs_A, -n_ko_ecs_B, -n_ko_rns_A, -n_ko_rns_B, -n_ko_paths_num_specific_A, -n_ko_paths_num_specific_B, -n_ko_brs_specific_A, -n_ko_brs_specific_B) %>%
			filter(n_tax >= edge_thre, edge_type == edgetype)
	}
	return(df_network)
}

read_network_data2 <- function (OG) {
	if (OG == "cog") {
		dt_net_edge_og <- read.csv("data/tsv/anticorrelogome_network_cog.tsv.gz", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>% as_tibble()
	} else if (OG == "kegg") {
		dt_net_edge_og <- read.csv("data/tsv/anticorrelogome_network_kegg.tsv.gz", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>% as_tibble() %>% 
		  mutate(og_id_A_old = sub("^","ko:",og_id_A), og_id_B_old = sub("^","ko:",og_id_B))
	}
	dt_net_edge_og$edge_type <- factor(dt_net_edge_og$edge_type, levels = c("Anticorrelating", "Correlating"))
	return(dt_net_edge_og)
}

read_network_data_each_taxon <- function (tax_x, edgetype, OG, g_whole) {
	#tax_x <- "p__Patescibacteria"
	#edgetype <- "Anticorrelating"
	#OG <- "cog"
	#g_whole <- g_cog_0_contra
	if (edgetype == "Anticorrelating" & OG == "cog") {
		df_network <- read.csv(paste0("data/contrapositive_pairs_cog2020/contrapositive_pairs_cog2020_", tax_x, ".tsv"), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("", "NULL")) %>%
		  as_tibble() %>%
		  mutate(
		    edge_type = "Anticorrelating"
		  )
	} else if (edgetype == "Correlating" & OG == "cog") {
		df_network <- read.csv(paste0("data/cooccurring_pairs_cog2020/cooccurring_pairs_cog2020_", tax_x, ".tsv"), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("", "NULL")) %>%
		  as_tibble() %>%
		  mutate(
		    edge_type = "Anticorrelating"
		  )
	} else if (edgetype == "Anticorrelating" & OG == "kegg") {
		df_network <- read.csv(paste0("data/contrapositive_pairs_ko/contrapositive_pairs_ko_", tax_x, ".tsv"), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("", "NULL")) %>%
		  as_tibble() %>%
		  rename(
		    og_id_A_old = og_id_A,
		    og_id_B_old = og_id_B,
		  ) %>%
		  mutate(
		    og_id_A = sub("^ko:", "", og_id_A_old),
		    og_id_B = sub("^ko:", "", og_id_B_old),
		    edge_type = "Anticorrelating"
		  ) %>%
		  unnest(edge_type) %>%
		  select(og_id_A, og_id_B, everything())
	} else if (edgetype == "Correlating" & OG == "kegg") {
		df_network <- read.csv(paste0("data/cooccurring_pairs_ko/cooccurring_pairs_ko_", tax_x, ".tsv"), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("", "NULL")) %>%
		  as_tibble() %>%
		  rename(
		    og_id_A_old = og_id_A,
		    og_id_B_old = og_id_B,
		  ) %>%
		  mutate(
		    og_id_A = sub("^ko:", "", og_id_A_old),
		    og_id_B = sub("^ko:", "", og_id_B_old),
		    edge_type = "Correlating"
		  ) %>%
		  unnest(edge_type) %>%
		  select(og_id_A, og_id_B, everything())
	}
	if (OG == "cog") {
		df_def <- read_cog2020_def_fixed_funcat_text() %>%
			right_join(
				tibble(cog_id = unique(c(df_network$og_id_A, df_network$og_id_B))),
				by = "cog_id"
			)
	} else if (OG == "kegg") {
		df_def <- read_ko_def_na() %>%
			rename(ko_id_old = "ko_id") %>%
			mutate(ko_id = sub("^ko:", "", ko_id_old)) %>%
			select(ko_id, everything()) %>%
			right_join(
				tibble(
					ko_id = unique(c(df_network$og_id_A, df_network$og_id_B)),
					ko_id_old = unique(c(df_network$og_id_A_old, df_network$og_id_B_old))
				),
				by = "ko_id_old"
			)
	}
	g_network <- df_network %>% select(og_id_A:tax) %>% graph_from_data_frame(directed = FALSE, vertices = df_def)
	V(g_network)$degree <- igraph::degree(g_network)
	name_custom <- V(g_whole)$name_custom
	names(name_custom) <- V(g_whole)$name
	V(g_network)$name_custom <- name_custom[V(g_network)$name]
	igraph::E(g_network)$n_tax <- igraph::E(g_whole)[as_ids(igraph::E(g_network))]$n_tax
	igraph::E(g_network)$n_tax_pseudo <- igraph::E(g_whole)[as_ids(igraph::E(g_network))]$n_tax
	igraph::E(g_network)$n_replace_prok <- igraph::E(g_whole)[as_ids(igraph::E(g_network))]$n_replace_prok
	igraph::E(g_network)$n_replace_prok_pseudo <- igraph::E(g_whole)[as_ids(igraph::E(g_network))]$n_replace_prok_pseudo
	igraph::E(g_network)$n_tax_pseudo[igraph::E(g_network)$n_tax_pseudo >= 10] <- 10
	return(g_network)
}

read_network_data_each_taxon_mini <- function (tax_x, edgetype, OG, g_whole) {
	#tax_x <- "p__Patescibacteria"
	#edgetype <- "Anticorrelating"
	#OG <- "cog"
	#g_whole <- g_cog_0_contra
	if (edgetype == "Anticorrelating" & OG == "cog") {
		df_network <- read.csv(paste0("data/contrapositive_pairs_cog2020/contrapositive_pairs_cog2020_", tax_x, ".tsv"), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("", "NULL")) %>%
		  as_tibble() %>%
		  mutate(
		    edge_type = "Anticorrelating"
		  )
	} else if (edgetype == "Correlating" & OG == "cog") {
		df_network <- read.csv(paste0("data/cooccurring_pairs_cog2020/cooccurring_pairs_cog2020_", tax_x, ".tsv"), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("", "NULL")) %>%
		  as_tibble() %>%
		  mutate(
		    edge_type = "Anticorrelating"
		  )
	} else if (edgetype == "Anticorrelating" & OG == "kegg") {
		df_network <- read.csv(paste0("data/contrapositive_pairs_ko/contrapositive_pairs_ko_", tax_x, ".tsv"), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("", "NULL")) %>%
		  as_tibble() %>%
		  rename(
		    og_id_A_old = og_id_A,
		    og_id_B_old = og_id_B,
		  ) %>%
		  mutate(
		    og_id_A = sub("^ko:", "", og_id_A_old),
		    og_id_B = sub("^ko:", "", og_id_B_old),
		    edge_type = "Anticorrelating"
		  ) %>%
		  unnest(edge_type) %>%
		  select(og_id_A, og_id_B, everything())
	} else if (edgetype == "Correlating" & OG == "kegg") {
		df_network <- read.csv(paste0("data/cooccurring_pairs_ko/cooccurring_pairs_ko_", tax_x, ".tsv"), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("", "NULL")) %>%
		  as_tibble() %>%
		  rename(
		    og_id_A_old = og_id_A,
		    og_id_B_old = og_id_B,
		  ) %>%
		  mutate(
		    og_id_A = sub("^ko:", "", og_id_A_old),
		    og_id_B = sub("^ko:", "", og_id_B_old),
		    edge_type = "Correlating"
		  ) %>%
		  unnest(edge_type) %>%
		  select(og_id_A, og_id_B, everything())
	}
	if (OG == "cog") {
		df_def <- read_cog2020_def_fixed_funcat_text() %>%
			right_join(
				tibble(cog_id = unique(c(df_network$og_id_A, df_network$og_id_B))),
				by = "cog_id"
			)
	} else if (OG == "kegg") {
		df_def <- read_ko_def_na() %>%
			rename(ko_id_old = "ko_id") %>%
			mutate(ko_id = sub("^ko:", "", ko_id_old)) %>%
			select(ko_id, everything()) %>%
			right_join(
				tibble(
					ko_id = unique(c(df_network$og_id_A, df_network$og_id_B)),
					ko_id_old = unique(c(df_network$og_id_A_old, df_network$og_id_B_old))
				),
				by = "ko_id_old"
			)
	}
	g_network <- df_network %>% select(og_id_A:tax) %>% graph_from_data_frame(directed = FALSE, vertices = df_def)
	V(g_network)$degree <- igraph::degree(g_network)
	name_custom <- V(g_whole)$name_custom
	names(name_custom) <- V(g_whole)$name
	V(g_network)$name_custom <- name_custom[V(g_network)$name]
	igraph::E(g_network)$n_tax <- igraph::E(g_whole)[as_ids(igraph::E(g_network))]$n_tax
	igraph::E(g_network)$n_tax_pseudo <- igraph::E(g_whole)[as_ids(igraph::E(g_network))]$n_tax
	igraph::E(g_network)$n_replace_prok <- igraph::E(g_whole)[as_ids(igraph::E(g_network))]$n_replace_prok
	igraph::E(g_network)$n_replace_prok_pseudo <- igraph::E(g_whole)[as_ids(igraph::E(g_network))]$n_replace_prok_pseudo
	igraph::E(g_network)$n_tax_pseudo[igraph::E(g_network)$n_tax_pseudo >= 10] <- 10
	assign_functional_categories <- function(ko_pathway) {
	  # Generate a canonical representation of pathway sets
	  categories <- ko_pathway %>%
	    filter(!is.na(pathway_id)) %>%
	    distinct(ko_id, pathway_id) %>%
	    group_by(ko_id) %>%
	    summarise(
	      pathway_set = paste(
	        sort(unique(pathway_id)),
	        collapse = ";"
	      ),
	      .groups = "drop"
	    ) %>%
	    arrange(pathway_set, ko_id) %>%
	    mutate(
	      category_id = sprintf(
	        "FC%05d",
	        match(pathway_set, unique(pathway_set))
	      )
	    )
	  return(categories %>%
	           select(ko_id, category_id))
	}
	df_path_num_FC <- 
	  assign_functional_categories(
	    tibble(
	      ko_id = V(g_network)$ko_id_old, 
	      pathway_id = V(g_network)$ko_path_num
	    )
	  )
	names(df_path_num_FC$category_id) <- df_path_num_FC$ko_id
	V(g_network)$ko_path_num <-
	  df_path_num_FC$category_id[V(g_network)$ko_id_old]
	attr_remove <- 
	  c("ko_id_old", "ko_id.y", "ko_cog", "ko_ec", "ko_path", "ko_rn", "ko_br", 
	    "ko_path_num_specific", "ko_br_specific")
	for (attr in attr_remove) {
	  g_network <- delete_vertex_attr(g_network, attr)
	}
	for (attr in c("og_id_A_old", "og_id_B_old", "cooccurring_ratio")) {
	  g_network <- delete_edge_attr(g_network, attr)
	}
	return(g_network)
}

read_gmgc10_env_gene_abundance <- function (OG) {
	if (OG == "cog") {
		df_gmgc <- read.csv("data/tsv/gmgc10_cog_colsum.tsv.gz", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
		  select(-amplicon,-isolate) %>% as_tibble()
	} else if (OG == "kegg") {
		df_gmgc <- read.csv("data/tsv/gmgc10_kegg_colsum.tsv.gz", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
		  select(-amplicon,-isolate) %>% as_tibble()
	}
	return(df_gmgc)
}

gmgc10_env_gene_relative_abundance <- function (df_gmgc_summary) {
	df_gmgc_summary_norm <- df_gmgc_summary %>%
    mutate(total_unigenes = rowSums(df_gmgc_summary %>% select(-OGs))) %>% 
    pivot_longer(cols = c(-OGs, -total_unigenes), names_to = "habitat", values_to = "n_unigenes") %>% 
    mutate(ratio_unigenes = n_unigenes / total_unigenes) %>% 
    select(-n_unigenes) %>% 
    pivot_wider(names_from = "habitat", values_from = "ratio_unigenes")
  return(df_gmgc_summary_norm)
}

# igraph

create_displacement_network_object <- function (OG, edge_thre, edgetype, dt_network, df_og_replacement_prok, replace_thre) {
	df_network <- read_network_edge_stats_general(OG, edge_thre, edgetype, dt_network) %>%
		mutate(edge_id = paste(og_id_A, og_id_B, sep = "_")) %>%
		left_join(
			df_og_replacement_prok %>% mutate(edge_id = paste(og_id_A, og_id_B, sep = "_")) %>% select(-og_id_A, -og_id_B),
			by = "edge_id"
		) %>%
		filter(n_replace_prok >= replace_thre)
	if (OG == "cog") {
		df_def <- read_cog2020_def_fixed_funcat_text() %>%
			right_join(
				tibble(cog_id = unique(c(df_network$og_id_A, df_network$og_id_B))),
				by = "cog_id"
			)
		g_network <<- df_network %>% graph_from_data_frame(directed = FALSE, vertices = df_def)
		V(g_network)$degree <- igraph::degree(g_network)
		V(g_network)$name_custom <- tibble(
		    name = V(g_network)$name,
		    cog_gene = V(g_network)$cog_gene,
		    cog_fun_cat = V(g_network)$cog_fun_cat
		  ) %>%
		  mutate(
		    name_2 = map2(name, cog_gene, ~ if (is.na(.y)) {.x} else {.y}),
		    name_custom = paste(name_2, cog_fun_cat, sep = "\n")
		  ) %>%
		  unnest(c(name_custom)) %>% select(name_custom) %>% pull()
		V(g_network)$name_custom_2 <- tibble(
		    name = V(g_network)$name,
		    cog_gene = V(g_network)$cog_gene,
		    cog_fun_cat = V(g_network)$cog_fun_cat
		  ) %>%
		  mutate(
		    name_2 = map2(name, cog_gene, ~ if (is.na(.y)) {.x} else {.y})
		  ) %>%
		  unnest(c(name_2)) %>% select(name_2) %>% pull()
		igraph::E(g_network)$n_replace_prok_pseudo <- tibble(
		    n_replace_prok = igraph::E(g_network)$n_replace_prok
		  ) %>%
		  mutate(
		    n_replace_prok_pseudo = map(n_replace_prok, ~ if (.x >= 100) {100} else {.x})
		  ) %>%
		  unnest(c(n_replace_prok_pseudo)) %>% select(n_replace_prok_pseudo) %>% pull()
	} else if (OG == "kegg") {
		df_def <- read_ko_def_na() %>%
			rename(ko_id_old = "ko_id") %>%
			mutate(ko_id = sub("^ko:", "", ko_id_old)) %>%
			select(ko_id, everything()) %>%
			right_join(
				tibble(
					ko_id = unique(c(df_network$og_id_A, df_network$og_id_B)),
					ko_id_old = unique(c(df_network$og_id_A_old, df_network$og_id_B_old))
				),
				by = "ko_id_old"
			)
		g_network <<- df_network %>% graph_from_data_frame(directed = FALSE, vertices = df_def)
		V(g_network)$degree <- igraph::degree(g_network)
		V(g_network)$name_custom <- tibble(
		    name = V(g_network)$name,
		    ko_symbol = V(g_network)$ko_symbol
		  ) %>%
		  mutate(
		    name_custom = map2(name, ko_symbol, ~ if (is.na(.y)) {.x} else {.y})
		  ) %>%
		  unnest(c(name_custom)) %>% select(name_custom) %>% pull()
		V(g_network)$name_custom_2 <- V(g_network)$name_custom
		igraph::E(g_network)$n_replace_prok_pseudo <- tibble(
		    n_replace_prok = igraph::E(g_network)$n_replace_prok
		  ) %>%
		  mutate(
		    n_replace_prok_pseudo = map(n_replace_prok, ~ if (.x >= 100) {100} else {.x})
		  ) %>%
		  unnest(c(n_replace_prok_pseudo)) %>% select(n_replace_prok_pseudo) %>% pull()
	}
	return(g_network)
}

create_displacement_network_object_kegg_mini <- function (OG, edge_thre, edgetype, dt_network, df_og_replacement_prok, replace_thre) {
	df_network <- read_network_edge_stats_general(OG, edge_thre, edgetype, dt_network) %>%
		mutate(edge_id = paste(og_id_A, og_id_B, sep = "_")) %>%
		left_join(
			df_og_replacement_prok %>% mutate(edge_id = paste(og_id_A, og_id_B, sep = "_")) %>% select(-og_id_A, -og_id_B),
			by = "edge_id"
		) %>%
		filter(n_replace_prok >= replace_thre)
	if (OG == "cog") {
		df_def <- read_cog2020_def_fixed_funcat_text() %>%
			right_join(
				tibble(cog_id = unique(c(df_network$og_id_A, df_network$og_id_B))),
				by = "cog_id"
			)
		g_network <<- df_network %>% graph_from_data_frame(directed = FALSE, vertices = df_def)
		V(g_network)$degree <- igraph::degree(g_network)
		V(g_network)$name_custom <- tibble(
		    name = V(g_network)$name,
		    cog_gene = V(g_network)$cog_gene,
		    cog_fun_cat = V(g_network)$cog_fun_cat
		  ) %>%
		  mutate(
		    name_2 = map2(name, cog_gene, ~ if (is.na(.y)) {.x} else {.y}),
		    name_custom = paste(name_2, cog_fun_cat, sep = "\n")
		  ) %>%
		  unnest(c(name_custom)) %>% select(name_custom) %>% pull()
		V(g_network)$name_custom_2 <- tibble(
		    name = V(g_network)$name,
		    cog_gene = V(g_network)$cog_gene,
		    cog_fun_cat = V(g_network)$cog_fun_cat
		  ) %>%
		  mutate(
		    name_2 = map2(name, cog_gene, ~ if (is.na(.y)) {.x} else {.y})
		  ) %>%
		  unnest(c(name_2)) %>% select(name_2) %>% pull()
		igraph::E(g_network)$n_replace_prok_pseudo <- tibble(
		    n_replace_prok = igraph::E(g_network)$n_replace_prok
		  ) %>%
		  mutate(
		    n_replace_prok_pseudo = map(n_replace_prok, ~ if (.x >= 100) {100} else {.x})
		  ) %>%
		  unnest(c(n_replace_prok_pseudo)) %>% select(n_replace_prok_pseudo) %>% pull()
	} else if (OG == "kegg") {
		df_def <- read_ko_def_na() %>%
			rename(ko_id_old = "ko_id") %>%
			mutate(ko_id = sub("^ko:", "", ko_id_old)) %>%
			select(ko_id, ko_symbol, ko_name, ko_id_old, ko_path_num) %>%
			right_join(
				tibble(
					ko_id = unique(c(df_network$og_id_A, df_network$og_id_B)),
					ko_id_old = unique(c(df_network$og_id_A_old, df_network$og_id_B_old))
				),
				by = "ko_id_old"
			)
		g_network <<- df_network %>% graph_from_data_frame(directed = FALSE, vertices = df_def)
		V(g_network)$degree <- igraph::degree(g_network)
		V(g_network)$name_custom <- tibble(
		    name = V(g_network)$name,
		    ko_symbol = V(g_network)$ko_symbol
		  ) %>%
		  mutate(
		    name_custom = map2(name, ko_symbol, ~ if (is.na(.y)) {.x} else {.y})
		  ) %>%
		  unnest(c(name_custom)) %>% select(name_custom) %>% pull()
		V(g_network)$name_custom_2 <- V(g_network)$name_custom
		igraph::E(g_network)$n_replace_prok_pseudo <- tibble(
		    n_replace_prok = igraph::E(g_network)$n_replace_prok
		  ) %>%
		  mutate(
		    n_replace_prok_pseudo = map(n_replace_prok, ~ if (.x >= 100) {100} else {.x})
		  ) %>%
		  unnest(c(n_replace_prok_pseudo)) %>% select(n_replace_prok_pseudo) %>% pull()
	}
	return(g_network)
}

gmgc10_env_gene_abundance_network_bray2 <- function (g_network, df_gmgc_summary) {
	df_gmgc_summary_norm <- df_gmgc_summary %>%
		group_by(OGs) %>%
		nest() %>%
		mutate(total_unigenes = map(data, ~ sum(.x$sum_unigenes))) %>%
		unnest(c(data,total_unigenes)) %>%
		mutate(norm_unigenes = sum_unigenes / total_unigenes) %>%
		select(-sum_unigenes,-total_unigenes) %>%
		pivot_wider(names_from = envs, values_from = norm_unigenes) %>%
		ungroup()
	ogs_disappear_in_gmgc <- setdiff(as_ids(V(g_network)), df_gmgc_summary_norm$OGs)
	df_gmgc_summary_norm <- bind_rows(
		df_gmgc_summary_norm,
		tibble(
			OGs = ogs_disappear_in_gmgc
		)
	)
	df_gmgc_summary_norm[is.na(df_gmgc_summary_norm)] <- 0
	df_gmgc_summary_norm_t <- df_gmgc_summary_norm %>%
		select(-OGs) %>%
		as.matrix() %>%
		t() %>%
		as.data.frame() %>%
		rownames_to_column() %>%
		as_tibble()
	colnames(df_gmgc_summary_norm_t) <- c("envs", df_gmgc_summary_norm$OGs)
	paired_ogs_bray <- function (og_a, og_b) {
		bray_pair <- matrix(
			c(df_gmgc_summary_norm_t[,og_a] %>% pull(), df_gmgc_summary_norm_t[,og_b] %>% pull()),
			nrow = 2,
			byrow = TRUE
		) %>%
		vegdist(method = "bray") %>%
		as.vector()
		return(bray_pair)
	}
	df_network_bray <- as_data_frame(g_network, what = "edges") %>%
		as_tibble() %>%
		select(from, to) %>%
		mutate(bray = map2(from, to, ~ paired_ogs_bray(.x, .y))) %>%
		unnest(c(bray))
	df_network_bray[df_network_bray$from %in% ogs_disappear_in_gmgc,"bray"] <- NA
	df_network_bray[df_network_bray$to %in% ogs_disappear_in_gmgc,"bray"] <- NA
	igraph::E(g_network)$bray <- df_network_bray$bray
	return(g_network)
}

gmgc10_network_bray <- function (g_network, df_gmgc_summary) {
  l_gmgc_summary <- pmap(df_gmgc_summary %>% select(-OGs,-total_unigenes), list)
  names(l_gmgc_summary) <- df_gmgc_summary$OGs
  paired_ogs_bray <- function (og_a, og_b) {
    bray_pair <- matrix(
      c(l_gmgc_summary[[og_a]] %>% unlist(), l_gmgc_summary[[og_b]] %>% unlist()),
      nrow = 2,
      byrow = TRUE
    ) %>%
      vegdist(method = "bray") %>%
      as.vector()
    return(bray_pair)
  }
  #l_gmgc_summary[["COG4985"]] %>% is.null() #%>% unlist()
  #l_gmgc_summary[["COG1851"]] %>% unlist() %>% is.nan() %>% sum()
  #paired_ogs_bray("COG0001","COG4985")
  # A few minutes step
  df_network_bray <- as_data_frame(g_network, what = "edges") %>%
    as_tibble() %>%
    select(from, to) %>%
    mutate(bray = map2(from, to, ~ 
                         # Skip OGs that are not found in eggnog (e.g. "COG4985")
                         if (is.null(l_gmgc_summary[[.x]]) || is.null(l_gmgc_summary[[.y]])) {
                           NA
                           # Skip OGs that have 0 abundance (e.g. "COG1851")
                         } else if (sum(is.nan(unlist(l_gmgc_summary[[.x]]))) > 0 || sum(is.nan(unlist(l_gmgc_summary[[.y]]))) > 0) {
                           NA
                         } else {
                           paired_ogs_bray(.x, .y)
                         }
    )) %>%
    unnest(c(bray))
  #df_network_bray[df_network_bray$from %in% ogs_disappear_in_gmgc,"bray"] <- NA
  #df_network_bray[df_network_bray$to %in% ogs_disappear_in_gmgc,"bray"] <- NA
  igraph::E(g_network)$bray <- df_network_bray$bray
  return(g_network)
}

# Cytoscape

cytoscape_save <- function (g_whole_displace, node_color_fun, net_title) {
	object_name_whole <- deparse(substitute(g_whole_displace))
	createNetworkFromIgraph(g_whole_displace, net_title)
	set.seed(1234)
	setLayoutProperties('kamada-kawai', list(singlePartition=1, randomize=0))
	layoutNetwork('kamada-kawai')
	setNodeShapeDefault('ellipse')
	setNodeBorderWidthDefault(4)
	setNodeBorderColorDefault("#FFFAFA")
	setEdgeLineWidthMapping('n_tax', widths=c(10,40))
	setNodeSizeMapping('degree', sizes=c(50,100))
	if (node_color_fun == "cog_fun_cat") {
		df_cog_node_color <- tibble(cog_fun_cat = V(g_whole_displace)$cog_fun_cat) %>%
		  distinct() %>% arrange(cog_fun_cat) %>%
		  mutate(
		    node_color = gsub("FF$","",viridis(length(unique((V(g_whole_displace)$cog_fun_cat)))))
		  )
		setNodeColorMapping(
		    table.column = 'cog_fun_cat',
		    table.column.values = unique(V(g_whole_displace)$cog_fun_cat),
		    colors = tibble(cog_fun_cat = unique(V(g_whole_displace)$cog_fun_cat)) %>%
		        left_join(df_cog_node_color, by = "cog_fun_cat") %>%
		        select(node_color) %>% pull(),
		    mapping.type='d'
		  )
	} else if (node_color_fun == "ko_path_num") {
		ko_path_nums <- V(g_whole_displace)$ko_path_num
		df_kegg_node_color <- tibble(ko_path_num = ko_path_nums[!is.na(ko_path_nums)]) %>%
		  distinct() %>% arrange(ko_path_num) %>%
		  mutate(
		    node_color = gsub("FF$","",viridis(length(unique((ko_path_nums[!is.na(ko_path_nums)])))))
		  )
		setNodeColorMapping(
		    table.column = 'ko_path_num',
		    table.column.values = unique(ko_path_nums[!is.na(ko_path_nums)]),
		    colors = tibble(ko_path_num = unique(ko_path_nums[!is.na(ko_path_nums)])) %>%
		        left_join(df_kegg_node_color, by = "ko_path_num") %>%
		        select(node_color) %>% pull(),
		    mapping.type='d',
		    default.color = "#7E7E7E"
		  )
	}
	setNodeLabelMapping('name_custom')
	edgetype <- sub("_.*", "", sub(".*[0-9]_", "", object_name_whole))
	OG <- sub("^g_", "", sub("_[0-9].*", "", object_name_whole))
	if (edgetype == "contra") {
		nettype <- "Anticorrelating"
		plasma_pseudo <- function (n, alpha = 1, begin = 0, end = 0.9, direction = 1, option = "C") {
		    if (begin < 0 | begin > 1 | end < 0 | end > 1) {
		        stop("begin and end must be in [0,1]")
		    }
		    if (abs(direction) != 1) {
		        stop("direction must be 1 or -1")
		    }
		    if (n == 0) {
		        return(character(0))
		    }
		    if (direction == -1) {
		        tmp <- begin
		        begin <- end
		        end <- tmp
		    }
		    option <- switch(EXPR = option, A = "A", magma = "A", B = "B", 
		        inferno = "B", C = "C", plasma = "C", D = "D", viridis = "D", 
		        E = "E", cividis = "E", F = "F", rocket = "F", G = "G", 
		        mako = "G", H = "H", turbo = "H", {
		            warning(paste0("Option '", option, "' does not exist. Defaulting to 'viridis'."))
		            "D"
		        })
		    map <- viridisLite::viridis.map[viridisLite::viridis.map$opt == 
		        option, ]
		    map_cols <- grDevices::rgb(map$R, map$G, map$B)
		    fn_cols <- grDevices::colorRamp(map_cols, space = "Lab", 
		        interpolate = "spline")
		    cols <- fn_cols(seq(begin, end, length.out = n))/255
		    sub("FF$", "", grDevices::rgb(cols[, 1], cols[, 2], cols[, 3], alpha = alpha))
		  }
		setEdgeColorMapping(
		    table.column = 'n_replace_prok_pseudo',
		    colors = plasma_pseudo,
		    mapping.type = "c"
		  )
	} else if (edgetype == "cooccu") {
		nettype <- "Correlating"
		if (OG == "cog") {
			df_edge_color <- tibble(V = igraph::E(g_whole_displace)$share_cog_fun_cat) %>%
				left_join(tibble(V = c(TRUE, FALSE), edge_color = c("#F79CD4", "#0FCFC0")), by = "V")
			setEdgeColorMapping(
			    table.column = "share_cog_fun_cat",
			    table.column.values = c(TRUE, FALSE),
			    colors = c("#0FCFC0", "#F79CD4"),
			    mapping.type = "d"
			  )
		} else if (OG == "kegg") {
			df_edge_color <- tibble(V = igraph::E(g_whole_displace)$share_ko_path_num_specific) %>%
				left_join(tibble(V = c(TRUE, FALSE), edge_color = c("#F79CD4", "#0FCFC0")), by = "V")
			setEdgeColorMapping(
			    table.column = "share_ko_path_num_specific",
			    table.column.values = c(TRUE, FALSE),
			    colors = c("#0FCFC0", "#F79CD4"),
			    mapping.type = "d"
			  )
		}
	}
	df_ndex_net <- tibble(
		`Network type` = nettype,
		OG = OG,
		`Network title (links)` = net_title,
		`# edges` = igraph::E(g_whole_displace) %>% length(),
		`# nodes` = igraph::V(g_whole_displace) %>% length()
	)
	return(df_ndex_net)
}

cytoscape_save2 <- function (g_whole_displace, node_color_fun, net_title, nettype, OG) {
	object_name_whole <- deparse(substitute(g_whole_displace))
	createNetworkFromIgraph(g_whole_displace, net_title)
	set.seed(1234)
	setLayoutProperties('kamada-kawai', list(singlePartition=1, randomize=0))
	layoutNetwork('kamada-kawai')
	setNodeShapeDefault('ellipse')
	setNodeBorderWidthDefault(4)
	setNodeBorderColorDefault("#FFFAFA")
	setEdgeLineWidthMapping('n_tax_pseudo', widths=c(10,40))
	setNodeSizeMapping('degree', sizes=c(50,100))
	if (node_color_fun == "cog_fun_cat") {
		df_cog_node_color <- tibble(cog_fun_cat = V(g_whole_displace)$cog_fun_cat) %>%
		  distinct() %>% arrange(cog_fun_cat) %>%
		  mutate(
		    node_color = gsub("FF$","",viridis(length(unique((V(g_whole_displace)$cog_fun_cat)))))
		  )
		setNodeColorMapping(
		    table.column = 'cog_fun_cat',
		    table.column.values = unique(V(g_whole_displace)$cog_fun_cat),
		    colors = tibble(cog_fun_cat = unique(V(g_whole_displace)$cog_fun_cat)) %>%
		        left_join(df_cog_node_color, by = "cog_fun_cat") %>%
		        select(node_color) %>% pull(),
		    mapping.type='d'
		  )
	} else if (node_color_fun == "ko_path_num") {
		ko_path_nums <- V(g_whole_displace)$ko_path_num
		df_kegg_node_color <- tibble(ko_path_num = ko_path_nums[!is.na(ko_path_nums)]) %>%
		  distinct() %>% arrange(ko_path_num) %>%
		  mutate(
		    node_color = gsub("FF$","",viridis(length(unique((ko_path_nums[!is.na(ko_path_nums)])))))
		  )
		setNodeColorMapping(
		    table.column = 'ko_path_num',
		    table.column.values = unique(ko_path_nums[!is.na(ko_path_nums)]),
		    colors = tibble(ko_path_num = unique(ko_path_nums[!is.na(ko_path_nums)])) %>%
		        left_join(df_kegg_node_color, by = "ko_path_num") %>%
		        select(node_color) %>% pull(),
		    mapping.type='d',
		    default.color = "#7E7E7E"
		  )
	}
	setNodeLabelMapping('name_custom')
	if (nettype == "Anticorrelating") {
		plasma_pseudo <- function (n, alpha = 1, begin = 0, end = 0.9, direction = 1, option = "C") {
		    if (begin < 0 | begin > 1 | end < 0 | end > 1) {
		        stop("begin and end must be in [0,1]")
		    }
		    if (abs(direction) != 1) {
		        stop("direction must be 1 or -1")
		    }
		    if (n == 0) {
		        return(character(0))
		    }
		    if (direction == -1) {
		        tmp <- begin
		        begin <- end
		        end <- tmp
		    }
		    option <- switch(EXPR = option, A = "A", magma = "A", B = "B", 
		        inferno = "B", C = "C", plasma = "C", D = "D", viridis = "D", 
		        E = "E", cividis = "E", F = "F", rocket = "F", G = "G", 
		        mako = "G", H = "H", turbo = "H", {
		            warning(paste0("Option '", option, "' does not exist. Defaulting to 'viridis'."))
		            "D"
		        })
		    map <- viridisLite::viridis.map[viridisLite::viridis.map$opt == 
		        option, ]
		    map_cols <- grDevices::rgb(map$R, map$G, map$B)
		    fn_cols <- grDevices::colorRamp(map_cols, space = "Lab", 
		        interpolate = "spline")
		    cols <- fn_cols(seq(begin, end, length.out = n))/255
		    sub("FF$", "", grDevices::rgb(cols[, 1], cols[, 2], cols[, 3], alpha = alpha))
		  }
		setEdgeColorMapping(
		    table.column = 'n_replace_prok_pseudo',
		    colors = plasma_pseudo,
		    mapping.type = "c"
		  )
	} else if (nettype == "Correlating") {
		if (OG == "cog") {
			df_edge_color <- tibble(V = igraph::E(g_whole_displace)$share_cog_fun_cat) %>%
				left_join(tibble(V = c(TRUE, FALSE), edge_color = c("#F79CD4", "#0FCFC0")), by = "V")
			setEdgeColorMapping(
			    table.column = "share_cog_fun_cat",
			    table.column.values = c(TRUE, FALSE),
			    colors = c("#0FCFC0", "#F79CD4"),
			    mapping.type = "d"
			  )
		} else if (OG == "kegg") {
			df_edge_color <- tibble(V = igraph::E(g_whole_displace)$share_ko_path_num_specific) %>%
				left_join(tibble(V = c(TRUE, FALSE), edge_color = c("#F79CD4", "#0FCFC0")), by = "V")
			setEdgeColorMapping(
			    table.column = "share_ko_path_num_specific",
			    table.column.values = c(TRUE, FALSE),
			    colors = c("#0FCFC0", "#F79CD4"),
			    mapping.type = "d"
			  )
		}
	}
	df_ndex_net <- tibble(
		OG = OG,
		`Network title` = net_title,
		`# edges` = igraph::E(g_whole_displace) %>% length(),
		`# nodes` = igraph::V(g_whole_displace) %>% length()
	)
	return(df_ndex_net)
}

subgraph_by_n_tax <- function (g_network, i) {
	g_subgraph <- subgraph.edges(g_network, as_ids(igraph::E(g_network)[igraph::E(g_network)$n_tax >= i]), delete.vertices = TRUE)
	return(g_subgraph)
}

cytoscape_to_ndex_whole <- function (g_whole_displace, node_color_fun, usr, password) {
	object_name_whole <- deparse(substitute(g_whole_displace))
	if (object_name_whole == "g_cog_8_contra_displace") {
		net_title <- "anticorrelating-cog-whole"
		createNetworkFromIgraph(g_whole_displace, net_title)
	} else if (object_name_whole == "g_kegg_8_contra_displace") {
		net_title <- "anticorrelating-kegg-whole"
		createNetworkFromIgraph(g_whole_displace, net_title)
	} else if (object_name_whole == "g_cog_80_cooccu") {
		net_title <- "correlating-cog-whole"
		createNetworkFromIgraph(g_whole_displace, net_title)
	} else if (object_name_whole == "g_kegg_80_cooccu") {
		net_title <- "correlating-kegg-whole"
		createNetworkFromIgraph(g_whole_displace, net_title)
	}
	setLayoutProperties('kamada-kawai', list(singlePartition=1, randomize=0))
	layoutNetwork('kamada-kawai')
	setNodeShapeDefault('ellipse')
	setNodeBorderWidthDefault(4)
	setNodeBorderColorDefault("#FFFAFA")
	setEdgeLineWidthMapping('n_tax', widths=c(10,40))
	setNodeSizeMapping('degree', sizes=c(50,100))
	if (node_color_fun == "cog_fun_cat") {
		df_cog_node_color <- tibble(cog_fun_cat = V(g_whole_displace)$cog_fun_cat) %>%
		  distinct() %>% arrange(cog_fun_cat) %>%
		  mutate(
		    node_color = gsub("FF$","",viridis(length(unique((V(g_whole_displace)$cog_fun_cat)))))
		  )
		setNodeColorMapping(
		    table.column = 'cog_fun_cat',
		    table.column.values = unique(V(g_whole_displace)$cog_fun_cat),
		    colors = tibble(cog_fun_cat = unique(V(g_whole_displace)$cog_fun_cat)) %>%
		        left_join(df_cog_node_color, by = "cog_fun_cat") %>%
		        select(node_color) %>% pull(),
		    mapping.type='d'
		  )
	} else if (node_color_fun == "ko_path_num") {
		ko_path_nums <- V(g_whole_displace)$ko_path_num
		df_kegg_node_color <- tibble(ko_path_num = ko_path_nums[!is.na(ko_path_nums)]) %>%
		  distinct() %>% arrange(ko_path_num) %>%
		  mutate(
		    node_color = gsub("FF$","",viridis(length(unique((ko_path_nums[!is.na(ko_path_nums)])))))
		  )
		setNodeColorMapping(
		    table.column = 'ko_path_num',
		    table.column.values = unique(ko_path_nums[!is.na(ko_path_nums)]),
		    colors = tibble(ko_path_num = unique(ko_path_nums[!is.na(ko_path_nums)])) %>%
		        left_join(df_kegg_node_color, by = "ko_path_num") %>%
		        select(node_color) %>% pull(),
		    mapping.type='d',
		    default.color = "#7E7E7E"
		  )
	}
	setNodeLabelMapping('name_custom')
	edgetype <- sub("_.*", "", sub(".*[0-9]_", "", object_name_whole))
	OG <- sub("^g_", "", sub("_[0-9].*", "", object_name_whole))
	if (edgetype == "contra") {
		nettype <- "Anticorrelating"
		plasma_pseudo <- function (n, alpha = 1, begin = 0, end = 0.9, direction = 1, option = "C") {
		    if (begin < 0 | begin > 1 | end < 0 | end > 1) {
		        stop("begin and end must be in [0,1]")
		    }
		    if (abs(direction) != 1) {
		        stop("direction must be 1 or -1")
		    }
		    if (n == 0) {
		        return(character(0))
		    }
		    if (direction == -1) {
		        tmp <- begin
		        begin <- end
		        end <- tmp
		    }
		    option <- switch(EXPR = option, A = "A", magma = "A", B = "B", 
		        inferno = "B", C = "C", plasma = "C", D = "D", viridis = "D", 
		        E = "E", cividis = "E", F = "F", rocket = "F", G = "G", 
		        mako = "G", H = "H", turbo = "H", {
		            warning(paste0("Option '", option, "' does not exist. Defaulting to 'viridis'."))
		            "D"
		        })
		    map <- viridisLite::viridis.map[viridisLite::viridis.map$opt == 
		        option, ]
		    map_cols <- grDevices::rgb(map$R, map$G, map$B)
		    fn_cols <- grDevices::colorRamp(map_cols, space = "Lab", 
		        interpolate = "spline")
		    cols <- fn_cols(seq(begin, end, length.out = n))/255
		    sub("FF$", "", grDevices::rgb(cols[, 1], cols[, 2], cols[, 3], alpha = alpha))
		  }
		setEdgeColorMapping(
		    table.column = 'n_replace_prok_pseudo',
		    colors = plasma_pseudo,
		    mapping.type = "c"
		  )
	} else if (edgetype == "cooccu") {
		nettype <- "Correlating"
		if (OG == "cog") {
			df_edge_color <- tibble(V = igraph::E(g_whole_displace)$share_cog_fun_cat) %>%
				left_join(tibble(V = c(TRUE, FALSE), edge_color = c("#F79CD4", "#0FCFC0")), by = "V")
			setEdgeColorMapping(
			    table.column = "share_cog_fun_cat",
			    table.column.values = c(TRUE, FALSE),
			    colors = c("#0FCFC0", "#F79CD4"),
			    mapping.type = "d"
			  )
		} else if (OG == "kegg") {
			df_edge_color <- tibble(V = igraph::E(g_whole_displace)$share_ko_path_num_specific) %>%
				left_join(tibble(V = c(TRUE, FALSE), edge_color = c("#F79CD4", "#0FCFC0")), by = "V")
			setEdgeColorMapping(
			    table.column = "share_ko_path_num_specific",
			    table.column.values = c(TRUE, FALSE),
			    colors = c("#0FCFC0", "#F79CD4"),
			    mapping.type = "d"
			  )
		}
	}
	uuid <- exportNetworkToNDEx(usr, password, FALSE)
	deleteAllNetworks()
	df_ndex_net <- tibble(
		`Network type` = nettype,
		OG = OG,
		`Network title (links)` = net_title,
		`# edges` = igraph::E(g_whole_displace) %>% length(),
		`# nodes` = igraph::V(g_whole_displace) %>% length(),
		`NDEx UUID` = uuid
	)
	return(df_ndex_net)
}

cytoscape_save_gmgc_bray <- function (g_graph, df_gmgc_summary_norm, network_title) {
  environ <- c(
    c("marine"),
    c("human.oral","human.skin","human.nose","human.vagina"),
    c("soil"),
    c("human.gut","pig.gut","mouse.gut","cat.gut","dog.gut"),
    c("wastewater","freshwater","built.environment")
  )
	df_node_def <- as_data_frame(g_graph, what = "vertices") %>%
		as_tibble() %>% 
		left_join(df_gmgc_summary_norm %>% select(OGs, total_unigenes, all_of(environ)) %>% dplyr::rename(name = "OGs"), by = "name")
	g_network <- graph_from_data_frame(as_data_frame(g_graph, what = "edges"), directed = FALSE, vertices = df_node_def)
	createNetworkFromIgraph(g_network, network_title)
	set.seed(1234)
	setLayoutProperties('kamada-kawai', list(singlePartition=1, randomize=0))
	layoutNetwork('kamada-kawai')
	setNodeShapeDefault('ellipse')
	setNodeCustomPieChart(
    environ,
		colors = c("#393b79","#637939","#8ca252","#b5cf6b","#cedb9c","#8c6d31","#843c39","#ad494a","#d6616b","#e7969c","#7b4173","#a55194","#ce6dbd")
	)
	setNodeBorderWidthDefault(4)
	setNodeBorderColorDefault("#FFFAFA")
	#setEdgeLineWidthDefault(10)
	setEdgeLineWidthMapping('n_tax', widths=c(10,40))
	viridis_pseudo <- function (n, alpha = 1, begin = 0, end = 0.9, direction = 1, option = "viridis") {
	    if (begin < 0 | begin > 1 | end < 0 | end > 1) {
	        stop("begin and end must be in [0,1]")
	    }
	    if (abs(direction) != 1) {
	        stop("direction must be 1 or -1")
	    }
	    if (n == 0) {
	        return(character(0))
	    }
	    if (direction == -1) {
	        tmp <- begin
	        begin <- end
	        end <- tmp
	    }
	    option <- switch(EXPR = option, A = "A", magma = "A", B = "B", 
	        inferno = "B", C = "C", plasma = "C", D = "D", viridis = "D", 
	        E = "E", cividis = "E", F = "F", rocket = "F", G = "G", 
	        mako = "G", H = "H", turbo = "H", {
	            warning(paste0("Option '", option, "' does not exist. Defaulting to 'viridis'."))
	            "D"
	        })
	    map <- viridisLite::viridis.map[viridisLite::viridis.map$opt == 
	        option, ]
	    map_cols <- grDevices::rgb(map$R, map$G, map$B)
	    fn_cols <- grDevices::colorRamp(map_cols, space = "Lab", 
	        interpolate = "spline")
	    cols <- fn_cols(seq(begin, end, length.out = n))/255
	    sub("FF$", "", grDevices::rgb(cols[, 1], cols[, 2], cols[, 3], alpha = alpha))
	  }
	setEdgeColorMapping(
	    table.column = 'bray',
	    colors = viridis_pseudo,
	    mapping.type = "c"
	  )
	setNodeSizeMapping('degree', sizes=c(50,100))
	#setNodeSizeDefault(70)
	setNodeColorDefault("#cccccc")
	setNodeLabelMapping('name_custom')
}
