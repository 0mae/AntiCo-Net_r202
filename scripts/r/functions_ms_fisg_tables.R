# Define functions used in the analysis

# Load libraries

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
library(ggpointdensity)
library(coRdon)
library(png)
library(ggrepel)
library(scatterpie)
library(vegan)
library(philentropy)
library(pals)
#library(bio3d)
#library(phytools)

# ggplot setting

## Chenge default setting

theme_set(
    theme_bw() +
    theme(
        text = element_text(family = "Arial", size = 12),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(face = "bold", hjust = 0, size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 18),
        plot.caption = element_text(hjust = 0.5, size = 10),
	plot.tag = element_text(face = "bold", size = 18),
	strip.text.x = element_text(size = 13),
	strip.text.y = element_text(size = 13),
    )
)

scientific_notation <- function(x) {
     x <- format(x, scientific = TRUE)
     x <- gsub("^(.*)e", "'\\1'e", x)
     x <- gsub("e", "%*%10^", x)
     x <- gsub('\\+', '', x)
     parse(text = x)
}

## Color blind Palette
### 8 colors
cbp8 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
### 12 colors
cbp12 <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

edge_thre4og_replacement_cooccu <- 6
edge_thre4og_replacement_contra <- 4
edge_thre4og_distance_cooccu <- 50
edge_thre4og_distance_contra <- 7
edge_thre_contra <- 8
edge_thre_cooccu <- 80
replace_thre <- 0

# For prepare COG category texts

df_fun20tab <- read.csv("db/cog2020/fun-20.tab", header = FALSE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>% as_tibble()
colnames(df_fun20tab) <- c("cog_fun_cat","cog_fun_cat_color","cog_fun_cat_text")
v_fun20tab <- df_fun20tab$cog_fun_cat_text
names(v_fun20tab) <- df_fun20tab$cog_fun_cat

# Create general objects

read_network_data <- function (OG) {
	if (OG == "cog") {
		dt_net_edge_og <- fread("data/tsv/network_edge_stats_cog.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("", "NULL")) %>%
		  as_tibble() %>%
		  rename(
		    edge_type_old = edge_type
		  ) %>%
		  mutate(
		    edge_type = map(edge_type_old, ~ if (.x == "contrapositive") {"Anticorrelating"} else if (.x == "cooccurring") {"Correlating"})
		  ) %>%
		  unnest(edge_type)
	} else if (OG == "kegg") {
		dt_net_edge_og <- fread("data/tsv/network_edge_stats_kegg.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("", "NULL")) %>%
		  as_tibble() %>%
		  rename(
		    og_id_A_old = og_id_A,
		    og_id_B_old = og_id_B,
		    #edge_type_old = edge_type,
		  ) %>%
		  mutate(
		    og_id_A = sub("^ko:", "", og_id_A_old),
		    og_id_B = sub("^ko:", "", og_id_B_old)
		    #edge_type = map(edge_type_old, ~ if (.x == "contrapositive") {"Anticorrelating"} else if (.x == "cooccurring") {"Correlating"})
		  ) %>%
		  unnest(edge_type) %>%
		  select(og_id_A, og_id_B, everything())
	}
	dt_net_edge_og$edge_type <- factor(dt_net_edge_og$edge_type, levels = c("Anticorrelating", "Correlating"))
	return(dt_net_edge_og)
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

read_cog2020_def_fixed <- function () {
	df_cog2020_def <- read.csv("db/cog2020_fixed/cog-20.def.utf8.tsv", header = FALSE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>% as_tibble()
	colnames(df_cog2020_def) <- c("cog_id", "cog_fun_cat", "cog_name", "cog_gene", "cog_fun_path", "cog_pubmed_id", "cog_pdb_id")
	return(df_cog2020_def)
}

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

read_ko_map_title <- function () {
	df_ko_map_title <- read.csv(
			"db/kegg_fixed/pathway_mod.list", header = FALSE, sep = "\t", quote = "",
			colClasses="character",	na.strings=c("", "NULL")
		) %>%
		rename(
			fun_single = "V1",
			fun_single_title = "V2",
			fun_single_higher_2 = "V3",
			fun_single_higher_1 = "V4",
		) %>%
		as_tibble()
	return(df_ko_map_title)
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

create_original_network_object <- function (OG, edge_thre, edgetype, dt_network) {
	df_network <- read_network_edge_stats_general(OG, edge_thre, edgetype, dt_network)
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
	}
	return(g_network)
}

create_original_network_object2 <- function (OG, edge_thre, edgetype, dt_network) {
  df_network <- read_network_edge_stats_general(OG, edge_thre, edgetype, dt_network)
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
  } else if (OG == "kegg") {
    df_def <- read_ko_def_na() %>%
      rename(ko_id_old = "ko_id") %>%
      mutate(ko_id = sub("^ko:", "", ko_id_old)) %>%
      select(ko_id, everything()) %>%
      right_join(
        tibble(
          ko_id = unique(c(df_network$og_id_A, df_network$og_id_B))
          #ko_id_old = unique(c(df_network$og_id_A_old, df_network$og_id_B_old))
        ),
        #by = "ko_id_old"
        by = "ko_id"
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
  }
  return(g_network)
}

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

edgeid2function <- function (g_graph, eids, OG) {
	froms <- sub("\\|.*", "", eids)
	tos <- sub(".*\\|", "", eids)
	g_subgraph <- subgraph.edges(g_graph, eids)
	if (OG == "cog"){
		fun_def <-  paste(
			V(g_subgraph)$cog_fun_cat,
			V(g_subgraph)$cog_gene,
			V(g_subgraph)$cog_name,
			sep = "; "
		)
		names(fun_def) <- V(g_subgraph)$name
	} else if (OG == "kegg") {
		fun_def <-  paste(
			V(g_subgraph)$ko_symbol,
			V(g_subgraph)$ko_name,
			sep = "; "
		)
		names(fun_def) <- V(g_subgraph)$name
	}
	df_og_table <- tibble(
		`OG A` = froms,
		`Name A` = fun_def[froms],
		`OG B` = tos,
		`Name B` = fun_def[tos]
	)
	return(df_og_table)
}

ggnetwork_fr <- function (g_graph) {
	set.seed(1)
	df_g_graph <- g_graph %>%
		ggnetwork(layout = igraph::layout_with_fr(g_graph)) %>%
		as_tibble() %>%
		mutate(
			x_same = x == xend,
			y_same = y == yend,
			v_e = map2(x_same, y_same, ~ if (.x == TRUE & .y == TRUE) {"vertex"} else {"edge"})
		) %>%
		unnest(c(v_e))
	return(df_g_graph)
}

ggnetwork_fr_bipartite <- function (sub_network) {
	set.seed(1)
	df_sub_network <- sub_network %>%
		ggnetwork(layout = igraph::layout_as_bipartite(sub_network)) %>%
		as_tibble() %>%
		mutate(
			x_same = x == xend,
			y_same = y == yend,
			v_e = map2(x_same, y_same, ~ if (.x == TRUE & .y == TRUE) {"vertex"} else {"edge"})
		) %>%
		unnest(c(v_e))
	return(df_sub_network)
}

kable_print <- function (df_kable, OG) {
	if (OG == "cog") {
		df_kable_print <- df_kable %>%
			select(name,cog_fun_cat,cog_gene,cog_name,cog_fun_path)
		colnames(df_kable_print) <- c("COG ID","Category","Symbol","Name","Pathway")
	} else if (OG == "kegg") {
		df_kable_print <- df_kable %>%
			#mutate(renamed = sub("^ko:","",name)) %>%
			#select(renamed,ko_symbol,ko_name)
			select(name,ko_symbol,ko_name)
		colnames(df_kable_print) <- c("KO ID","Symbol","Name")
	}
	df_kable_print <- df_kable_print %>% knitr::kable()
	return(df_kable_print)
}

kable_print_additional <- function (df_kable, OG, add_col) {
	if (OG == "cog") {
		df_kable_print <- df_kable %>%
			select(.data[[add_col]], name,cog_fun_cat,cog_gene,cog_name,cog_fun_path) %>%
			arrange(desc(.data[[add_col]]))
		colnames(df_kable_print) <- c(add_col, "COG ID","Category","Symbol","Name","Pathway")
	} else if (OG == "kegg") {
		df_kable_print <- df_kable %>%
			#mutate(renamed = sub("^ko:","",name)) %>%
			#select(.data[[add_col]],renamed,ko_symbol,ko_name) %>%
			select(.data[[add_col]],name,ko_symbol,ko_name) %>%
			arrange(desc(.data[[add_col]]))
		colnames(df_kable_print) <- c(add_col, "KO ID","Symbol","Name")
	}
	df_kable_print <- df_kable_print %>% knitr::kable()
	return(df_kable_print)
}

# Number of edges

network_n_edges <- function (dt_net_edge) {
	df_net_edge <- dt_net_edge %>% as_tibble()
	df_n_edge <- df_net_edge %>%
		group_by(edge_type, n_tax) %>%
		summarize(n_edge = dplyr::n()) %>%
		ungroup()
	p_n_edge <- df_n_edge %>%
		ggplot(aes(x = n_tax, y = n_edge, color = edge_type)) +
		geom_point() +
		annotation_logticks(sides = "l") +
		scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
		scale_color_manual(values = hcl.colors(2,"TealRose"), name = "") +
		xlim(0, 296) +
		geom_vline(xintercept = edge_thre_contra, size = 0.2, color = "black") +
		labs(x = "Number of taxa", y = "Number of edges") +
		guides(
			color = guide_legend(
				override.aes = list(size=3),
				ncol = 1
			),
		) +
		theme(
			legend.position = "bottom",
			legend.title=element_blank(),
			legend.text = element_text(size = 12),
		)
	return(p_n_edge)
}

network_n_edges_anti <- function (dt_net_edge_cog, dt_net_edge_kegg) {
	df_n_edge <- rbind(
			dt_net_edge_cog %>% filter(edge_type == "Anticorrelating") %>% select(n_tax) %>% mutate(OG = "COG"),
			dt_net_edge_kegg %>% filter(edge_type == "Anticorrelating") %>% select(n_tax) %>% mutate(OG = "KEGG")
		) %>%
		group_by(n_tax, OG) %>%
		summarize(n_edge = dplyr::n()) %>%
		ungroup()
	p_n_edge <- df_n_edge %>%
		ggplot(aes(x = n_tax, y = n_edge, color = OG, shape = OG)) +
		geom_point() +
		annotation_logticks(sides = "l") +
		scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
		scale_color_manual(values = hcl.colors(2,"Earth")) +
		geom_vline(xintercept = edge_thre_contra, size = 0.2, color = "red") +
		#facet_grid(OG ~ .) +
		labs(x = "Number of taxa", y = "Number of edges") +
		theme(
			#legend.position = "bottom",
			#legend.position = c(1,1),
			#legend.justification = c(1.1,1.1),
			legend.title=element_blank(),
			legend.text = element_text(size = 12),
		)
	return(p_n_edge)
}


network_n_edges_anti_log <- function (dt_net_edge_cog, dt_net_edge_kegg) {
  df_n_edge <- rbind(
    dt_net_edge_cog %>% filter(edge_type == "Anticorrelating") %>% select(n_tax) %>% mutate(OG = "COG"),
    dt_net_edge_kegg %>% filter(edge_type == "Anticorrelating") %>% select(n_tax) %>% mutate(OG = "KEGG")
  ) %>%
    group_by(n_tax, OG) %>%
    summarize(n_edge = dplyr::n()) %>%
    ungroup()
  p_n_edge <- df_n_edge %>%
    ggplot(aes(x = n_tax, y = n_edge, color = OG, shape = OG)) +
    geom_point() +
    annotation_logticks(sides = "lb") +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_color_manual(values = hcl.colors(2,"Earth")) +
    #geom_vline(xintercept = edge_thre_contra, size = 0.2, color = "red") +
    #facet_grid(OG ~ .) +
    labs(x = "Number of taxa", y = "Number of edges") +
    theme(
      #legend.position = "bottom",
      #legend.position = c(1,1),
      #legend.justification = c(1.1,1.1),
      legend.title=element_blank(),
      legend.text = element_text(size = 12),
    )
  return(p_n_edge)
}

# Assortativity of edges

network_structure_boxplot <- function(net_indicator) {
	df_net_struct_og <- read.csv("data/tsv/network_structure_cog.tsv", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
		as_tibble() %>%
		filter(tax_level != "Kingdom")
	df_net_struct_og$tax_level <- factor(df_net_struct_og$tax_level, levels = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))
	p_net_struct_og <- ggplot(
			data = df_net_struct_og,
			aes(x = assortativity_degree, y = edge_type, color = edge_type),
		) +
		geom_jitter() +
		geom_boxplot(
			fill = NA,
			color = "grey20",
			alpha = 0.5,
			outlier.colour = NA
		) +
		scale_color_manual(values = cbp12, na.value = "black") +
		theme(
			axis.title.y = element_blank(),
			axis.text.y = element_blank(),
			legend.position = "none",
			axis.ticks.y = element_blank(),
		) +
		labs(x = "Assortativity of degree")
	return(p_net_struct_og)
}

# OG distance and displacement

load_og_replacement <- function (OG, dom) {
	df_og_replacement <- read.table(
		paste0("data/tsv/og_displacement_", dom, "_", OG, ".tsv"),
		header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")
	) %>% as_tibble()
	return(df_og_replacement)
}

load_og_replacement2 <- function (OG, dom) {
  df_og_replacement <- read.table(
    paste0("data/tsv/og_displacement_", dom, "_", OG, ".tsv.gz"),
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

combine_with_net_edge_stats2 <- function (dt_net_edge_og, df_for_edge) {
  df_net_edge_og_mod <- right_join(
    as_tibble(dt_net_edge_og) %>% mutate(edge_id = paste(og_id_A, og_id_B, sep = "_")),
    df_for_edge %>% mutate(edge_id = paste(og_id_A, og_id_B, sep = "_")) %>% select(-og_id_A, -og_id_B),
    by = "edge_id"
  )
  return(df_net_edge_og_mod)
}

network_n_replace <- function (dt_net_edge_og, OG, df_og_replacement_prok) {
	dt_net_edge_og <<- dt_net_edge_og
	df_net_edge_og_replace <- combine_with_net_edge_stats(df_og_replacement_prok)
	p_n_replace <- df_net_edge_og_replace %>%
		filter(edge_type == "Anticorrelating" & n_tax >= edge_thre4og_replacement_contra) %>%
		ggplot(aes(x = n_tax, y = n_replace_prok)) +
		geom_pointdensity() +
		geom_vline(xintercept = edge_thre_contra, size = 0.2, color = "red") +
		scale_color_viridis() +
		theme(
			legend.position = "right",
			#legend.position = c(1,0),
			#legend.justification = c(1,0),
			legend.key.height= unit(12, 'pt')
		) +
		labs(x = "Number of taxa", y = "Displacement", color = "Density")
	return(p_n_replace)
}

network_n_replace2 <- function (dt_net_edge_og, OG, df_og_replacement_prok) {
  dt_net_edge_og <<- dt_net_edge_og
  df_net_edge_og_replace <- combine_with_net_edge_stats(df_og_replacement_prok)
  p_n_replace <- df_net_edge_og_replace %>%
    filter(n_tax >= 8 | n_replace_prok >= 50) %>%
    ggplot(aes(x = n_tax, y = n_replace_prok)) +
    geom_pointdensity(size = 0.01) +
    geom_vline(xintercept = 8, size = 0.1, color = "red") +
    geom_hline(yintercept = 50, size = 0.1, color = "red") +
    scale_color_viridis() +
    theme(
      legend.position = "right",
      #legend.position = c(1,0),
      #legend.justification = c(1,0),
      legend.key.height= unit(12, 'pt')
    ) +
    labs(x = "Number of taxa", y = "Displacement", color = "Density")
  return(p_n_replace)
}

network_n_replace_table <- function (dt_net_edge_og, OG, df_og_replacement_prok) {
	dt_net_edge_og <<- dt_net_edge_og
	df_net_edge_og_replace <- combine_with_net_edge_stats(df_og_replacement_prok)
	df_kable <- df_net_edge_og_replace %>%
		filter(edge_type == "Anticorrelating" & n_tax >= edge_thre_contra) %>%
		mutate(n_replace_bin = cut(n_replace_prok, c(1,100,max(df_net_edge_og_replace$n_replace_prok)), right = FALSE, include.lowest = TRUE)) %>%
		group_by(n_replace_bin) %>%
		summarize(n_edges = dplyr::n()) %>%
		ungroup() %>% mutate(OG = OG)
	return(df_kable)
}

network_n_replace_shared_fun <- function (dt_net_edge_og, OG, df_og_replacement_prok, n_tax_thre, fun_profile) {
	dt_net_edge_og <<- dt_net_edge_og
	df_net_edge_og_replace <- combine_with_net_edge_stats(df_og_replacement_prok)
	p_n_replace <- df_net_edge_og_replace %>%
		filter(edge_type == "Anticorrelating" & n_tax >= n_tax_thre) %>%
		ggplot(aes(x = n_tax, y = n_replace_prok)) +
		geom_pointdensity() +
		geom_hline(yintercept = 100, size = 0.2, color = "red") +
		scale_color_viridis() +
		facet_grid(. ~ .data[[fun_profile]]) +
		theme(
			legend.position = "right",
			legend.key.height= unit(12, 'pt')
		) +
		labs(x = "Number of taxa", y = "Displacement")
	return(p_n_replace)
}

network_distance_cog <- function () {
	df_og_distance <- read.csv("data/tsv/og_distance_cog.tsv", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
		as_tibble()
	df_net_edge_og_distance <- combine_with_net_edge_stats(df_og_distance)
	p_og_distance <- df_net_edge_og_distance %>%
		filter(edge_type == "cooccurring" & n_tax >= edge_thre4og_distance_cooccu) %>%
		ggplot(aes(x = n_tax, y = min_intergenic_mean)) +
		geom_pointdensity() +
		scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
		annotation_logticks(sides = "l") +
		geom_vline(xintercept = edge_thre_cooccu, size = 0.2, color = "red") +
		#geom_hline(yintercept = 10000, size = 0.1, color = "red") +
		scale_color_viridis() +
		theme(legend.key.height= unit(12, 'pt')) +
		labs(x = "Number of taxa", y = "Distance (bp)", color = "Density")
	return(p_og_distance)
}

# Contrapositive & cooccurring networks

network_cog_graph2df <- function (g_network, col_name) {
	set.seed(1)
	df_g_network <- g_network %>%
		ggnetwork(layout = igraph::layout_with_fr(g_network)) %>%
		as_tibble() %>%
		mutate(
			x_same = x == xend,
			y_same = y == yend,
			v_e = map2(x_same, y_same, ~ if (.x == TRUE & .y == TRUE) {"vertex"} else {"edge"}),
			v_shape = map(.data[[col_name]], ~ if (length(grep("(R|S)", .x)) >= 1) {"23"} else {"21"}),
			v_color = .data[[col_name]]
		) %>%
		unnest(c(v_e, v_shape))
	return(df_g_network)
}

network_cog_graph2df_bipartite <- function (g_network, col_name) {
  set.seed(1)
  df_g_network <- g_network %>%
    ggnetwork(layout = igraph::layout_as_bipartite(g_network)) %>%
    as_tibble() %>%
    mutate(
      x_same = x == xend,
      y_same = y == yend,
      v_e = map2(x_same, y_same, ~ if (.x == TRUE & .y == TRUE) {"vertex"} else {"edge"}),
      v_shape = map(.data[[col_name]], ~ if (length(grep("(R|S)", .x)) >= 1) {"23"} else {"21"}),
      v_color = .data[[col_name]]
    ) %>%
    unnest(c(v_e, v_shape))
  return(df_g_network)
}

ggplot_network <- function(df_g_network, col_name, category) {
	if (category == "discrete") {
		p_network <- ggplot() +
			geom_segment(
				aes(x = x, y = y, xend = xend, yend = yend, color = .data[[col_name]], size = n_tax),
				data = df_g_network %>% filter(v_e == "edge")
			) +
			scale_size_continuous(range = c(1,2.8)) +
			guides(size = "none") +
			new_scale("size") +
			geom_point(
				#aes(x = x, y = y, fill = v_color, shape = v_shape),
				aes(x = x, y = y, fill = v_color, size = degree),
				shape = 21,
				#size = 1.8,
				color = "white",
				data = df_g_network %>% filter(v_e == "vertex")
			) +
			scale_fill_viridis_d(na.value = "grey50") +
			scale_size_continuous(range = c(1.8,3.5)) +
			#scale_shape_manual(values = c(21, 23)) +
			#scale_color_manual(values = cbp8, na.value = "grey50") +
			scale_color_manual(values = rev(hcl.colors(2,"Cyan-Magenta")), na.value = "grey50") +
			#scale_color_manual(values = rev(hcl.colors(2,"Berlin")), na.value = "grey50") +
			coord_fixed(ratio = 1) +
			guides(fill = "none", shape = "none")
	} else if (category == "continuous") {
		p_network <- ggplot() +
			geom_segment(
				aes(x = x, y = y, xend = xend, yend = yend, color = .data[[col_name]], size = n_tax),
				data = df_g_network %>% filter(v_e == "edge")
			) +
			scale_size_continuous(range = c(1,2.8)) +
			guides(size = "none") +
			new_scale("size") +
			geom_point(
				#aes(x = x, y = y, fill = v_color, shape = v_shape),
				aes(x = x, y = y, fill = v_color, size = degree),
				shape = 21,
				#size = 1.8,
				color = "white",
				data = df_g_network %>% filter(v_e == "vertex")
			) +
			scale_fill_viridis_d(na.value = "grey50") +
			scale_size_continuous(range = c(1.8,3.5)) +
			#scale_shape_manual(values = c(21, 23)) +
			scale_color_viridis_c(option = "plasma", begin = 0, end = 0.9) +
			coord_fixed(ratio = 1) +
			guides(fill = "none", shape = "none")
	}
	return(p_network)
}

contrapositive_network_show <- function (g_network, fun_profile) {
	df_g_network <- network_cog_graph2df(g_network, fun_profile)
	df_g_network$n_replace_prok[df_g_network$n_replace_prok >= 100] <- 100
	p_network <- ggplot_network(df_g_network, "n_replace_prok", "continuous") +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
        		plot.title = element_text(face = "bold", hjust = 0, size = 24),
			plot.subtitle = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			legend.position = "bottom"
		) +
		guides(size = "none") +
		labs(subtitle = "Anticorrelating", color = "Displacement")
	return(p_network)
}

cooccurring_network_show <- function (g_network, fun_profile, share_fun_profile) {
	df_g_network <- network_cog_graph2df(g_network, fun_profile)
	p_network <- ggplot_network(df_g_network, share_fun_profile, "discrete") +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
        		plot.title = element_text(face = "bold", hjust = 0, size = 24),
			plot.subtitle = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			legend.position = "bottom"
		) +
		guides(size = "none") +
		labs(subtitle = "Correlating")
	return(p_network)
}

add_geom_text_repel <- function (p_network, g_network, fun_profile, OGs) {
	p_network_labeled <- p_network +
		geom_text_repel(
			aes(x = x, y = y, label = name_custom_2),
			min.segment.length = 0,
			segment.size = 0.1,
			segment.linetype = 1,
			segment.color = "black",
			size = 2.5,
			data = network_cog_graph2df(g_network, fun_profile) %>%
				filter(v_e == "vertex") %>%
				distinct() %>%
				inner_join(
					tibble(name = OGs),
					by = "name"
				)
		)
	return(p_network_labeled)
}

add_geom_text_repel_align <- function (p_network, g_network, fun_profile, OGs) {
	p_network_labeled <- p_network +
		geom_text_repel(
			aes(x = x, y = y, label = name_custom_2),
			min.segment.length = 0,
			segment.size = 0.1,
			segment.linetype = 1,
			segment.color = "black",
			nudge_x = 1.02,
			direction = "y",
			hjust = 1,
			size = 3,
			data = network_cog_graph2df(g_network, fun_profile) %>%
				filter(v_e == "vertex") %>%
				distinct() %>%
				inner_join(
					tibble(name = OGs),
					by = "name"
				)
		)
	return(p_network_labeled)
}

read_genomes_analyzed <- function() {
	df_genomes_analyzed <- rbind(
		read.csv("data/tsv/genomes_analyzed_ar.tsv", header = TRUE, sep = "\t", quote = "") %>% as_tibble(),
		read.csv("data/tsv/genomes_analyzed_bac.tsv", header = TRUE, sep = "\t", quote = "") %>% as_tibble()
	)
	return(df_genomes_analyzed)
}

plot_tree_legend_taxa <- function (dom, tax_level) {
	if (dom == "archaea") {
		tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Archaea"
	} else if (dom == "bacteria") {
		tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Bacteria"
	}
	gtdb_node_labeled_tree <- read.tree(tree_path)
	p_tree <- ggtree(gtdb_node_labeled_tree, branch.length = 'none', layout = 'circular', color = "gray")
	df_genomes <- read_genomes_analyzed() %>% select(accession, Domain:Species) %>% rename(label = "accession")
	df_tax <- df_genomes[which(df_genomes$label %in% p_tree$data[p_tree$data$isTip,]$label), c("label", tax_level)]
	df_tax_top <- df_tax %>% group_by(.data[[tax_level]]) %>% summarize(count = dplyr::n()) %>% arrange(desc(count)) %>%
		head(11) %>% select(-count) %>% mutate(tax_color = .data[[tax_level]]) %>% bind_rows(tibble(tax_color = "Others")) %>%
		cbind(tibble(x = 0.2, y_pseudo = seq(12,1)) %>% mutate(y = y_pseudo / 12)) %>% as_tibble()
	df_tax_top$tax_color <- factor(df_tax_top$tax_color, levels = c(df_tax_top$tax_color))
	return(df_tax_top)
}

visualize_og_displacement <- function (OG, dom, og_a, og_b, gene_a, gene_b, tax_level) {
	if (OG == "cog" & dom == "archaea") {
		og_a_path <- paste0("data/pastml_cog2020_ar/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_cog2020_ar/", og_b, ".txt")
		tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Archaea"
		node_plot_size <- 1.4
	} else if (OG == "cog" & dom == "bacteria") {
		og_a_path <- paste0("data/pastml_cog2020_bac/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_cog2020_bac/", og_b, ".txt")
		tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Bacteria"
		node_plot_size <- 0.8
	} else if (OG == "kegg" & dom == "archaea") {
		og_a_path <- paste0("data/pastml_ko_ar/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_ko_ar/", og_b, ".txt")
		tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Archaea"
		node_plot_size <- 1.4
	} else if (OG == "kegg" & dom == "bacteria") {
		og_a_path <- paste0("data/pastml_ko_bac/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_ko_bac/", og_b, ".txt")
		tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Bacteria"
		node_plot_size <- 0.8
	}
	df_og_a <- read.csv(og_a_path, header = TRUE, sep = "\t", quote = "") %>% as_tibble()
	df_og_b <- read.csv(og_b_path, header = TRUE, sep = "\t", quote = "") %>% as_tibble()
	gtdb_node_labeled_tree <- read.tree(tree_path)
	df_genomes <- read_genomes_analyzed() %>% select(accession, Domain:Species) %>% rename(label = "accession")
	df_og_a_pivot <- df_og_a %>% mutate(presence = 1) %>% group_by(node) %>%
		pivot_wider(names_from = all_of(og_a), values_from = presence, values_fill = 0) %>% rename(label = "node", og_a_TRUE = "TRUE", og_a_FALSE = "FALSE")
	df_og_b_pivot <- df_og_b %>% mutate(presence = 1) %>% group_by(node) %>%
		pivot_wider(names_from = all_of(og_b), values_from = presence, values_fill = 0) %>% rename(label = "node", og_b_TRUE = "TRUE", og_b_FALSE = "FALSE")
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
		unnest(c(og_a_presence, og_b_presence))
	p_tree <- ggtree(gtdb_node_labeled_tree, branch.length = 'none', layout = 'circular', color = "gray", size = 0.27)
	df_tip_data <- df_tree_mapped[which(df_tree_mapped$label %in% gtdb_node_labeled_tree$tip.label), c("label", "og_a_presence", "og_b_presence")] %>%
		pivot_longer(-label, names_to = "ogs_tip", values_to = "pastml_tip") %>%
		mutate(
			og_presence = pastml_tip >= 1,
			OGs_tip = map(ogs_tip, ~ if (.x == "og_a_presence") {paste0(og_a, "\n", gene_a)} else if (.x == "og_b_presence") {paste0(og_b, "\n", gene_b)})
			#OGs_tip = map(ogs_tip, ~ if (.x == "og_a_presence") {gene_a} else if (.x == "og_b_presence") {gene_b})
		) %>%
		unnest(c(OGs_tip)) %>%
		left_join(p_tree$data, by = "label") %>%
		rename(ID = "label")
	df_tip_data$OGs_tip <- factor(df_tip_data$OGs_tip, levels = c(paste0(og_a, "\n", gene_a), paste0(og_b, "\n", gene_b)))
	#df_tip_data$OGs_tip <- factor(df_tip_data$OGs_tip, levels = c(paste0(gene_a), paste0(gene_b)))
	df_tax <- df_genomes[which(df_genomes$label %in% p_tree$data[p_tree$data$isTip,]$label), c("label", tax_level)]
	df_tax_top <- df_tax %>% group_by(.data[[tax_level]]) %>% summarize(count = dplyr::n()) %>% arrange(desc(count)) %>%
		head(11) %>% select(-count) %>% mutate(tax_color = .data[[tax_level]])
	df_tax <- df_tax %>% left_join(df_tax_top, by = tax_level)
	df_tax$tax_color[is.na(df_tax$tax_color)] <- "Others"
	df_tax$tax_color <- factor(df_tax$tax_color, levels = c(df_tax_top$tax_color, "Others"))
	df_tax_data <- left_join(p_tree$data, df_tax, by = "label") %>% filter(isTip == TRUE)
	df_node_data <- df_tree_mapped[which(df_tree_mapped$label %in% gtdb_node_labeled_tree$node.label), c("label", "og_a_presence", "og_b_presence")] %>%
		pivot_longer(-label, names_to = "ogs_node", values_to = "pastml_node") %>%
		mutate(
			node_size = factor(pastml_node),
			OGs_node = map(ogs_node, ~ if (.x == "og_a_presence") {paste0(og_a, "\n", gene_a)} else if (.x == "og_b_presence") {paste0(og_b, "\n", gene_b)})
		) %>%
		unnest(c(OGs_node)) %>%
		left_join(p_tree$data, by = "label")
	df_node_data$node_size[df_node_data$pastml_node == 0] <- NA
	df_node_data$OGs_node <- factor(df_node_data$OGs_node, levels = c(paste0(og_a, "\n", gene_a), paste0(og_b, "\n", gene_b)))
	p_tree_plotted <- p_tree +
		geom_point(
			data = df_node_data,
			aes(x = x, y = y, color = OGs_node, size = node_size),
			shape = 16,
			alpha = 0.5
		) +
		guides(size = "none") +
		scale_size_manual(values = c(0, node_plot_size/2, node_plot_size), na.translate = FALSE, drop = FALSE) +
		scale_color_manual(values = hcl.colors(2,"Blue-Red 2"), guide = "none") +
		new_scale_fill() +
		geom_fruit(
			data = df_tip_data %>% filter(og_presence == TRUE) %>% select(ID, OGs_tip),
			geom = geom_tile,
			mapping = aes(y = ID, x = OGs_tip, fill = OGs_tip),
			color = NA,
			#alpha = 0.8,
			pwidth = 0.25,
			offset = 0.08
		) +
		labs(fill = "OGs") +
		scale_fill_manual(values = hcl.colors(2,"Blue-Red 2")) +
		new_scale_fill() +
		geom_fruit(
			data = df_tax_data %>% rename(ID = "label") %>% select(ID, tax_color) %>% mutate(x_pos = 1),
			geom = geom_star,
			mapping = aes(y = ID, x = x_pos, fill = tax_color),
			starshape = 25,
			color = NA,
			size = 1.5,
			#alpha = 0.8,
			#pwidth = 0.04,
			offset = 0.04
		) +
		guides(fill = "none") +
		scale_fill_manual(values = cbp12, na.value = cbp12[12]) +
		labs(size = "Presence\nprobability", fill = tax_level) +
    		theme_void() +
		#coord_fixed(ratio = 1) +
    		theme(
			plot.background = element_rect(fill = "white",color = NA),
    			text = element_text(family = "Arial", size = 12),
    			axis.title.x = element_text(size = 13),
    			axis.title.y = element_text(size = 13),
    			plot.title = element_text(face = "bold", hjust = 0, size = 18),
    			plot.subtitle = element_text(hjust = 0.5, size = 12),
    			plot.caption = element_text(hjust = 0.5, size = 10),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13),
			legend.position = "top",
			legend.title = element_blank()
		)
}

visualize_og_displacement_rect <- function (OG, dom, og_a, og_b, gene_a, gene_b, tax_level) {
	if (OG == "cog" & dom == "archaea") {
		og_a_path <- paste0("data/pastml_cog2020_ar/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_cog2020_ar/", og_b, ".txt")
		tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Archaea"
		node_plot_size <- 1.4
	} else if (OG == "cog" & dom == "bacteria") {
		og_a_path <- paste0("data/pastml_cog2020_bac/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_cog2020_bac/", og_b, ".txt")
		tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Bacteria"
		node_plot_size <- 0.8
	} else if (OG == "kegg" & dom == "archaea") {
		og_a_path <- paste0("data/pastml_ko_ar/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_ko_ar/", og_b, ".txt")
		tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Archaea"
		node_plot_size <- 1.4
	} else if (OG == "kegg" & dom == "bacteria") {
		og_a_path <- paste0("data/pastml_ko_bac/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_ko_bac/", og_b, ".txt")
		tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Bacteria"
		node_plot_size <- 0.8
	}
	df_og_a <- read.csv(og_a_path, header = TRUE, sep = "\t", quote = "") %>% as_tibble()
	df_og_b <- read.csv(og_b_path, header = TRUE, sep = "\t", quote = "") %>% as_tibble()
	gtdb_node_labeled_tree <- read.tree(tree_path)
	df_genomes <- read_genomes_analyzed() %>% select(accession, Domain:Species) %>% rename(label = "accession")
	df_og_a_pivot <- df_og_a %>% mutate(presence = 1) %>% group_by(node) %>%
		pivot_wider(names_from = all_of(og_a), values_from = presence, values_fill = 0) %>% rename(label = "node", og_a_TRUE = "TRUE", og_a_FALSE = "FALSE")
	df_og_b_pivot <- df_og_b %>% mutate(presence = 1) %>% group_by(node) %>%
		pivot_wider(names_from = all_of(og_b), values_from = presence, values_fill = 0) %>% rename(label = "node", og_b_TRUE = "TRUE", og_b_FALSE = "FALSE")
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
		unnest(c(og_a_presence, og_b_presence))
	p_tree <- ggtree(gtdb_node_labeled_tree, branch.length = 'none', layout = 'circular', color = "gray", size = 0.5)
	df_tip_data <- df_tree_mapped[which(df_tree_mapped$label %in% gtdb_node_labeled_tree$tip.label), c("label", "og_a_presence", "og_b_presence")] %>%
		pivot_longer(-label, names_to = "ogs_tip", values_to = "pastml_tip") %>%
		mutate(
			og_presence = pastml_tip >= 1,
			OGs_tip = map(ogs_tip, ~ if (.x == "og_a_presence") {paste0(og_a, "\n", gene_a)} else if (.x == "og_b_presence") {paste0(og_b, "\n", gene_b)})
			#OGs_tip = map(ogs_tip, ~ if (.x == "og_a_presence") {gene_a} else if (.x == "og_b_presence") {gene_b})
		) %>%
		unnest(c(OGs_tip)) %>%
		left_join(p_tree$data, by = "label") %>%
		rename(ID = "label")
	df_tip_data$OGs_tip <- factor(df_tip_data$OGs_tip, levels = c(paste0(og_a, "\n", gene_a), paste0(og_b, "\n", gene_b)))
	#df_tip_data$OGs_tip <- factor(df_tip_data$OGs_tip, levels = c(paste0(gene_a), paste0(gene_b)))
	df_tax <- df_genomes[which(df_genomes$label %in% p_tree$data[p_tree$data$isTip,]$label), c("label", tax_level)]
	df_tax_top <- df_tax %>% group_by(.data[[tax_level]]) %>% summarize(count = dplyr::n()) %>% arrange(desc(count)) %>%
		head(11) %>% select(-count) %>% mutate(tax_color = .data[[tax_level]])
	df_tax <- df_tax %>% left_join(df_tax_top, by = tax_level)
	df_tax$tax_color[is.na(df_tax$tax_color)] <- "Others"
	df_tax$tax_color <- factor(df_tax$tax_color, levels = c(df_tax_top$tax_color, "Others"))
	df_tax_data <- left_join(p_tree$data, df_tax, by = "label") %>% filter(isTip == TRUE)
	df_node_data <- df_tree_mapped[which(df_tree_mapped$label %in% gtdb_node_labeled_tree$node.label), c("label", "og_a_presence", "og_b_presence")] %>%
		pivot_longer(-label, names_to = "ogs_node", values_to = "pastml_node") %>%
		mutate(
			node_size = factor(pastml_node),
			OGs_node = map(ogs_node, ~ if (.x == "og_a_presence") {paste0(og_a, "\n", gene_a)} else if (.x == "og_b_presence") {paste0(og_b, "\n", gene_b)})
		) %>%
		unnest(c(OGs_node)) %>%
		left_join(p_tree$data, by = "label")
	df_node_data$node_size[df_node_data$pastml_node == 0] <- NA
	df_node_data$OGs_node <- factor(df_node_data$OGs_node, levels = c(paste0(og_a, "\n", gene_a), paste0(og_b, "\n", gene_b)))
	p_tree_plotted <- p_tree +
		geom_point(
			data = df_node_data,
			aes(x = x, y = y, color = OGs_node, size = node_size),
			shape = 16,
			alpha = 0.5
		) +
		guides(size = "none") +
		scale_size_manual(values = c(0, node_plot_size/2, node_plot_size), na.translate = FALSE, drop = FALSE) +
		scale_color_manual(values = hcl.colors(2,"Blue-Red 2"), guide = "none") +
		new_scale_fill() +
		geom_fruit(
			data = df_tip_data %>% filter(og_presence == TRUE) %>% select(ID, OGs_tip),
			geom = geom_tile,
			mapping = aes(y = ID, x = OGs_tip, fill = OGs_tip),
			color = NA,
			#alpha = 0.8,
			pwidth = 0.4,
			offset = 0.14
		) +
		labs(fill = "OGs") +
		scale_fill_manual(values = hcl.colors(2,"Blue-Red 2")) +
		new_scale_fill() +
		geom_fruit(
			#data = df_tax_data %>% rename(ID = "label") %>% select(ID, tax_color) %>% mutate(x_pos = 1),
			data = rbind(
				df_tax_data %>% rename(ID = "label") %>% select(ID, tax_color) %>% mutate(x_pos = 1),
				df_tax_data %>% rename(ID = "label") %>% select(ID, tax_color) %>% mutate(x_pos = 2)
			),
			#geom = geom_star,
			geom = geom_tile,
			mapping = aes(y = ID, x = x_pos, fill = tax_color),
			#axis.params = list(axis = "x"),
			#starshape = 15,
			color = NA,
			#size = 2.5,
			#alpha = 0.8,
			pwidth = 0.12,
			offset = 0.30
		) +
		guides(fill = "none") +
		scale_fill_manual(values = cbp12, na.value = cbp12[12]) +
		labs(size = "Presence\nprobability", fill = tax_level) +
    		theme_void() +
  		layout_rectangular() +
		coord_flip() +
		scale_y_reverse() +
		#coord_fixed(ratio = 1) +
    		theme(
			plot.background = element_rect(fill = "white",color = NA),
    			text = element_text(family = "Arial", size = 12),
    			axis.title.x = element_text(size = 13),
    			axis.title.y = element_text(size = 13),
    			plot.title = element_text(face = "bold", hjust = 0, size = 18),
    			plot.subtitle = element_text(hjust = 0.5, size = 12),
    			plot.caption = element_text(hjust = 0.5, size = 10),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13),
			legend.position = "right",
			legend.title = element_blank()
		)
}

visualize_og_displacement_rect_tax <- function (OG, dom, og_a, og_b, gene_a, gene_b, tax_level) {
	if (OG == "cog" & dom == "archaea") {
		og_a_path <- paste0("data/pastml_cog2020_ar/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_cog2020_ar/", og_b, ".txt")
		tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Archaea"
		node_plot_size <- 1.4
	} else if (OG == "cog" & dom == "bacteria") {
		og_a_path <- paste0("data/pastml_cog2020_bac/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_cog2020_bac/", og_b, ".txt")
		tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Bacteria"
		node_plot_size <- 0.8
	} else if (OG == "kegg" & dom == "archaea") {
		og_a_path <- paste0("data/pastml_ko_ar/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_ko_ar/", og_b, ".txt")
		tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Archaea"
		node_plot_size <- 1.4
	} else if (OG == "kegg" & dom == "bacteria") {
		og_a_path <- paste0("data/pastml_ko_bac/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_ko_bac/", og_b, ".txt")
		tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Bacteria"
		node_plot_size <- 0.8
	}
	df_og_a <- read.csv(og_a_path, header = TRUE, sep = "\t", quote = "") %>% as_tibble()
	df_og_b <- read.csv(og_b_path, header = TRUE, sep = "\t", quote = "") %>% as_tibble()
	gtdb_node_labeled_tree <- read.tree(tree_path)
	df_genomes <- read_genomes_analyzed() %>% select(accession, Domain:Species) %>% rename(label = "accession")
	df_og_a_pivot <- df_og_a %>% mutate(presence = 1) %>% group_by(node) %>%
		pivot_wider(names_from = all_of(og_a), values_from = presence, values_fill = 0) %>% rename(label = "node", og_a_TRUE = "TRUE", og_a_FALSE = "FALSE")
	df_og_b_pivot <- df_og_b %>% mutate(presence = 1) %>% group_by(node) %>%
		pivot_wider(names_from = all_of(og_b), values_from = presence, values_fill = 0) %>% rename(label = "node", og_b_TRUE = "TRUE", og_b_FALSE = "FALSE")
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
		unnest(c(og_a_presence, og_b_presence))
	p_tree <- ggtree(gtdb_node_labeled_tree, branch.length = 'none', layout = 'circular', color = "gray", size = 0.5)
	df_tip_data <- df_tree_mapped[which(df_tree_mapped$label %in% gtdb_node_labeled_tree$tip.label), c("label", "og_a_presence", "og_b_presence")] %>%
		pivot_longer(-label, names_to = "ogs_tip", values_to = "pastml_tip") %>%
		mutate(
			og_presence = pastml_tip >= 1,
			OGs_tip = map(ogs_tip, ~ if (.x == "og_a_presence") {paste0(og_a, "\n", gene_a)} else if (.x == "og_b_presence") {paste0(og_b, "\n", gene_b)})
			#OGs_tip = map(ogs_tip, ~ if (.x == "og_a_presence") {gene_a} else if (.x == "og_b_presence") {gene_b})
		) %>%
		unnest(c(OGs_tip)) %>%
		left_join(p_tree$data, by = "label") %>%
		rename(ID = "label")
	df_tip_data$OGs_tip <- factor(df_tip_data$OGs_tip, levels = c(paste0(og_a, "\n", gene_a), paste0(og_b, "\n", gene_b)))
	#df_tip_data$OGs_tip <- factor(df_tip_data$OGs_tip, levels = c(paste0(gene_a), paste0(gene_b)))
	df_tax <- df_genomes[which(df_genomes$label %in% p_tree$data[p_tree$data$isTip,]$label), c("label", tax_level)]
	df_tax_top <- df_tax %>% group_by(.data[[tax_level]]) %>% summarize(count = dplyr::n()) %>% arrange(desc(count)) %>%
		head(11) %>% select(-count) %>% mutate(tax_color = .data[[tax_level]])
	df_tax <- df_tax %>% left_join(df_tax_top, by = tax_level)
	df_tax$tax_color[is.na(df_tax$tax_color)] <- "Others"
	df_tax$tax_color <- factor(df_tax$tax_color, levels = c(df_tax_top$tax_color, "Others"))
	df_tax_data <- left_join(p_tree$data, df_tax, by = "label") %>% filter(isTip == TRUE)
	df_node_data <- df_tree_mapped[which(df_tree_mapped$label %in% gtdb_node_labeled_tree$node.label), c("label", "og_a_presence", "og_b_presence")] %>%
		pivot_longer(-label, names_to = "ogs_node", values_to = "pastml_node") %>%
		mutate(
			node_size = factor(pastml_node),
			OGs_node = map(ogs_node, ~ if (.x == "og_a_presence") {paste0(og_a, "\n", gene_a)} else if (.x == "og_b_presence") {paste0(og_b, "\n", gene_b)})
		) %>%
		unnest(c(OGs_node)) %>%
		left_join(p_tree$data, by = "label")
	df_node_data$node_size[df_node_data$pastml_node == 0] <- NA
	df_node_data$OGs_node <- factor(df_node_data$OGs_node, levels = c(paste0(og_a, "\n", gene_a), paste0(og_b, "\n", gene_b)))
	p_tree_plotted <- p_tree +
		geom_point(
			data = df_node_data,
			aes(x = x, y = y, color = OGs_node, size = node_size),
			shape = 16,
			alpha = 0.5
		) +
		guides(size = "none") +
		scale_size_manual(values = c(0, node_plot_size/2, node_plot_size), na.translate = FALSE, drop = FALSE) +
		scale_color_manual(values = hcl.colors(2,"Blue-Red 2"), guide = "none") +
		new_scale_fill() +
		geom_fruit(
			data = df_tip_data %>% filter(og_presence == TRUE) %>% select(ID, OGs_tip),
			geom = geom_tile,
			mapping = aes(y = ID, x = OGs_tip, fill = OGs_tip),
			color = NA,
			#alpha = 0.8,
			pwidth = 0.4,
			offset = 0.14
		) +
		labs(fill = "OGs") +
		scale_fill_manual(values = hcl.colors(2,"Blue-Red 2")) +
		new_scale_fill() +
		geom_fruit(
			#data = df_tax_data %>% rename(ID = "label") %>% select(ID, tax_color) %>% mutate(x_pos = 1),
			data = rbind(
				df_tax_data %>% rename(ID = "label") %>% select(ID, tax_color) %>% mutate(x_pos = 1),
				df_tax_data %>% rename(ID = "label") %>% select(ID, tax_color) %>% mutate(x_pos = 2)
			),
			#geom = geom_star,
			geom = geom_tile,
			mapping = aes(y = ID, x = x_pos, fill = tax_color),
			#axis.params = list(axis = "x"),
			#starshape = 15,
			color = NA,
			#size = 2.5,
			#alpha = 0.8,
			pwidth = 0.12,
			offset = 0.30
		) +
		#guides(fill = "none") +
		scale_fill_manual(values = cbp12, na.value = cbp12[12]) +
		labs(size = "Presence\nprobability", fill = tax_level) +
    		theme_void() +
  		layout_rectangular() +
		coord_flip() +
		scale_y_reverse() +
		#coord_fixed(ratio = 1) +
    		theme(
			plot.background = element_rect(fill = "white",color = NA),
    			text = element_text(family = "Arial", size = 12),
    			axis.title.x = element_text(size = 13),
    			axis.title.y = element_text(size = 13),
    			plot.title = element_text(face = "bold", hjust = 0, size = 18),
    			plot.subtitle = element_text(hjust = 0.5, size = 12),
    			plot.caption = element_text(hjust = 0.5, size = 10),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13),
			legend.position = "bottom",
			legend.title = element_blank()
		)
}

vis_edge_stats_divide_n_edge <- function (dt_net_edge_og, edge_thre_contra, edge_thre_cooccu, OG, direct) {
	df_net_edge_og <- dt_net_edge_og %>%
		select(og_id_A,og_id_B,edge_type,n_tax) %>%
		as_tibble()
	if (direct == "h") {
		p_n_edges_div <- df_net_edge_og %>%
			mutate(n_tax_bin = cut(n_tax, c(1,8,max(df_net_edge_og$n_tax)), right = FALSE, include.lowest = TRUE)) %>%
			group_by(edge_type,n_tax_bin) %>%
			summarize(n_edges = dplyr::n()) %>%
			mutate(OG = OG) %>%
			ggplot(aes(x = n_edges, y = n_tax_bin, fill = edge_type)) +
			geom_bar(stat = "identity", position = "dodge") +
			scale_fill_manual(values = cbp12) +
			scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
			annotation_logticks(sides = "b") +
			labs(y = "Number of taxa", x = "Number of edges") +
			facet_grid(. ~ OG) +
			theme(legend.position = "none")
	} else if (direct == "v") {
		p_n_edges_div <- df_net_edge_og %>%
			mutate(n_tax_bin = cut(n_tax, c(1,8,max(df_net_edge_og$n_tax)), right = FALSE, include.lowest = TRUE)) %>%
			group_by(edge_type,n_tax_bin) %>%
			summarize(n_edges = dplyr::n()) %>%
			mutate(OG = OG) %>%
			ggplot(aes(y = n_edges, x = n_tax_bin, fill = edge_type)) +
			geom_bar(stat = "identity", position = "dodge") +
			scale_fill_manual(values = cbp12) +
			scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
			annotation_logticks(sides = "l") +
			labs(y = "Number of taxa", x = "Number of edges") +
			facet_grid(OG ~ .) +
			theme(legend.position = "none")
	}
	return(p_n_edges_div)
}

dataframe_edge_stats_div <- function(dt_net_edge_og, shared_profile) {
	if (shared_profile == "share_cog_fun_cat") {
		shared_profile_name = "COG Category"
	} else if (shared_profile == "share_cog_fun_path") {
		shared_profile_name = "COG Pathway"
	} else if (shared_profile == "share_ko_path_num_specific") {
		shared_profile_name = "KEGG Pathway"
	} else if (shared_profile == "share_ko_rn") {
		shared_profile_name = "KEGG Reaction"
	} else if (shared_profile == "share_ko_ec1") {
		shared_profile_name = "EC (1st digit)"
	} else if (shared_profile == "share_ko_ec2") {
		shared_profile_name = "EC (2 digits)"
	} else if (shared_profile == "share_ko_ec3") {
		shared_profile_name = "EC (3 digits)"
	} else if (shared_profile == "share_ko_ec") {
		#shared_profile_name = "EC (full id)"
		shared_profile_name = "KEGG Enzyme"
	}
	df_net_edge_og <- dt_net_edge_og %>% as_tibble()
	df_net_edge_og_div <- df_net_edge_og %>%
		mutate(n_tax_bin = cut(n_tax, c(1,8,max(df_net_edge_og$n_tax)), right = FALSE, include.lowest = TRUE)) %>%
		group_by(edge_type,n_tax_bin,.data[[shared_profile]]) %>%
		summarize(n_edges = dplyr::n()) %>%
		mutate(function_type = shared_profile_name) %>%
		rename(shared_profile = all_of(shared_profile)) %>%
		ungroup() %>%
		group_by(edge_type,n_tax_bin) %>%
		nest() %>%
		mutate(n_edges_total = map(data, ~ sum(.x$n_edges))) %>%
		unnest(c(data,n_edges_total)) %>%
		ungroup() %>%
		mutate(edge_ratio = n_edges / n_edges_total) %>%
		filter(shared_profile == TRUE)
	return(df_net_edge_og_div)
}

edge_tax_summary <- function (OG, df_taxa_summarized_50) {
	if (OG == "cog") {
		df_edge_tax <- fread("data/tsv/network_edge_stats_cog_taxa.tsv", stringsAsFactors=FALSE,header = TRUE, sep = "\t", quote = "") %>%
			as_tibble()
	} else if (OG == "kegg") {
		df_edge_tax <- fread("data/tsv/network_edge_stats_kegg_taxa.tsv", stringsAsFactors=FALSE,header = TRUE, sep = "\t", quote = "") %>%
			as_tibble()
	}
	df_edge_tax_summary <- df_edge_tax %>%
		filter(n_tax >= 8) %>%
		mutate(tax = strsplit(taxa, ";")) %>%
		unnest(c(tax)) %>%
		left_join(df_taxa_summarized_50, by = "tax")
	return(df_edge_tax_summary)
}

vis_edge_stats_div <- function(df_edge_og_div, ali) {
	if (ali == "v") {
		p_edge_kegg_div <- df_edge_og_div %>%
			ggplot(aes(y = edge_ratio, x = n_tax_bin, fill = edge_type)) +
			geom_bar(stat = "identity", position = "dodge") +
			facet_grid(function_type ~ .) +
			scale_fill_manual(values = hcl.colors(2,"TealRose")) +
			labs(y = "Number of taxa", x = "Ratio of edges") +
			theme(legend.position = "none")
	} else if (ali == "h") {
		p_edge_kegg_div <- df_edge_og_div %>%
			ggplot(aes(x = edge_ratio, y = n_tax_bin, fill = edge_type)) +
			geom_bar(stat = "identity", position = "dodge") +
			facet_grid(. ~ function_type) +
			scale_fill_manual(values = hcl.colors(2,"TealRose")) +
			labs(y = "Number of taxa", x = "Ratio of edges") +
			theme(legend.position = "none")
	}
	return(p_edge_kegg_div)
}

vis_edge_fun_sharing <- function(df_edge_og, shared_profile) {
	p_edge <- df_edge_og %>%
		filter(!is.na(.data[[shared_profile]])) %>%
		mutate(edge_type_short = sub("Correlating", "Corr.", sub("Anticorrelating", "Anti.", edge_type))) %>%
		group_by(n_tax,edge_type_short,.data[[shared_profile]]) %>%
		summarize(n_edges = dplyr::n()) %>%
		ggplot(aes(y = n_edges, x = n_tax, fill = .data[[shared_profile]])) +
		geom_bar(stat = "identity", position = "fill") +
		geom_vline(xintercept = edge_thre_contra, size = 0.2, color = "black") +
		facet_grid(edge_type_short ~ .) +
		scale_fill_manual(values = hcl.colors(2,"TealRose")) +
		labs(x = "Number of taxa", y = "Ratio of edges", fill = "Sharing\nfunction") +
		theme(strip.text.y = element_text(size = 13))
	return(p_edge)
}

vis_edge_fun_sharing_anti <- function(df_edge_og, shared_profiles) {
	df_plot <- data.frame()
	for (shared_profile in shared_profiles) {
		df_plot <- rbind(
			df_plot,
			df_edge_og %>%
				filter(edge_type == "Anticorrelating", !is.na(.data[[shared_profile]])) %>%
				rename(share_tf = all_of(shared_profile)) %>%
				group_by(n_tax,share_tf) %>%
				summarize(n_edges = dplyr::n()) %>%
				mutate(fun_profile = shared_profile)
		)
	}
	df_plot <- df_plot %>%
		mutate(
			fun_profile_short = map(fun_profile, ~ 
				if (.x == "share_cog_fun_cat") {
					"COG Category"
				} else if (.x == "share_cog_fun_path") {
					"COG Pathway"
				} else if (.x == "share_ko_path_num_specific") {
					"KEGG\nPathway"
				} else if (.x == "share_ko_rn") {
					"KEGG\nReaction"
				} else if (.x == "share_ko_ec1") {
					"EC\n1 digit"
				} else if (.x == "share_ko_ec2") {
					"EC\n2 digits"
				} else if (.x == "share_ko_ec3") {
					"EC\n3 digits"
				} else if (.x == "share_ko_ec") {
					"EC number"
				}
			)
		) %>%
		unnest(c(fun_profile_short))
	df_plot$fun_profile_short <- factor(df_plot$fun_profile_short, levels = c("COG Category","COG Pathway","KEGG\nPathway","KEGG\nReaction","EC\n1 digit","EC\n2 digits","EC\n3 digits","EC number"))
	p_edge <- df_plot %>%
		ggplot(aes(y = n_edges, x = n_tax, fill = share_tf)) +
		geom_bar(stat = "identity", position = "fill") +
		geom_vline(xintercept = edge_thre_contra, size = 0.2, color = "black") +
		scale_fill_manual(values = hcl.colors(2,"TealRose")) +
		facet_grid(fun_profile_short ~ .) +
		labs(x = "Number of taxa", y = "Ratio of edges", fill = "Sharing\nfunction") +
		theme(strip.text.y = element_text(size = 13))
	return(p_edge)
}

vis_edge_fun_sharing_anti_displacement <- function(df_edge_og, shared_profiles) {
	df_plot <- data.frame()
	for (shared_profile in shared_profiles) {
		df_plot <- rbind(
			df_plot,
			df_edge_og %>%
				filter(edge_type == "Anticorrelating", !is.na(.data[[shared_profile]])) %>%
				rename(share_tf = all_of(shared_profile)) %>%
				group_by(n_replace_prok,share_tf) %>%
				summarize(n_edges = dplyr::n()) %>%
				mutate(fun_profile = shared_profile)
		)
	}
	df_plot <- df_plot %>%
		mutate(
			fun_profile_short = map(fun_profile, ~ 
				if (.x == "share_cog_fun_cat") {
					"COG Category"
				} else if (.x == "share_cog_fun_path") {
					"COG Pathway"
				} else if (.x == "share_ko_path_num_specific") {
					"KEGG\nPathway"
				} else if (.x == "share_ko_rn") {
					"KEGG\nReaction"
				} else if (.x == "share_ko_ec1") {
					"EC\n1 digit"
				} else if (.x == "share_ko_ec2") {
					"EC\n2 digits"
				} else if (.x == "share_ko_ec3") {
					"EC\n3 digits"
				} else if (.x == "share_ko_ec") {
					"EC number"
				}
			)
		) %>%
		unnest(c(fun_profile_short))
	df_plot$fun_profile_short <- factor(df_plot$fun_profile_short, levels = c("COG Category","COG Pathway","KEGG\nPathway","KEGG\nReaction","EC\n1 digit","EC\n2 digits","EC\n3 digits","EC number"))
	p_edge <- df_plot %>%
		ggplot(aes(y = n_edges, x = n_replace_prok, fill = share_tf)) +
		geom_bar(stat = "identity", position = "fill") +
		#geom_vline(xintercept = 90, size = 0.2, color = "black") +
		scale_fill_manual(values = hcl.colors(2,"TealRose")) +
		facet_grid(fun_profile_short ~ .) +
		labs(x = "Displacement frequency", y = "Ratio of edges", fill = "Sharing\nfunction") +
		theme(strip.text.y = element_text(size = 13))
	return(p_edge)
}

vis_edge_fun_sharing_anti_displacement2 <- function(df_edge_og, shared_profiles) {
  df_plot <- data.frame()
  for (shared_profile in shared_profiles) {
    df_plot <- rbind(
      df_plot,
      df_edge_og %>%
        filter(edge_type == "Anticorrelating", !is.na(.data[[shared_profile]])) %>%
        rename(share_tf = all_of(shared_profile)) %>%
        group_by(n_tax,share_tf) %>%
        summarize(n_edges = dplyr::n()) %>%
        mutate(fun_profile = shared_profile)
    )
  }
  df_plot <- df_plot %>%
    mutate(
      fun_profile_short = map(fun_profile, ~ 
                                if (.x == "share_cog_fun_cat") {
                                  "COG Category"
                                } else if (.x == "share_cog_fun_path") {
                                  "COG Pathway"
                                } else if (.x == "share_ko_path_num_specific") {
                                  "KEGG\nPathway"
                                } else if (.x == "share_ko_rn") {
                                  "KEGG\nReaction"
                                } else if (.x == "share_ko_ec1") {
                                  "EC\n1 digit"
                                } else if (.x == "share_ko_ec2") {
                                  "EC\n2 digits"
                                } else if (.x == "share_ko_ec3") {
                                  "EC\n3 digits"
                                } else if (.x == "share_ko_ec") {
                                  "EC number"
                                }
      )
    ) %>%
    unnest(c(fun_profile_short))
  df_plot$fun_profile_short <- factor(df_plot$fun_profile_short, levels = c("COG Category","COG Pathway","KEGG\nPathway","KEGG\nReaction","EC\n1 digit","EC\n2 digits","EC\n3 digits","EC number"))
  p_edge <- df_plot %>%
    ggplot(aes(y = n_edges, x = n_tax, fill = share_tf)) +
    geom_bar(stat = "identity", position = "fill") +
    geom_vline(xintercept = edge_thre_contra, size = 0.2, color = "black") +
    scale_fill_manual(values = hcl.colors(2,"TealRose")) +
    facet_grid(fun_profile_short ~ .) +
    labs(x = "Number of taxa", y = "Ratio of edges", fill = "Sharing\nfunction") +
    theme(strip.text.y = element_text(size = 13))
  return(p_edge)
}

vis_edge_bray_n_tax <- function(g_network) {
  df_plot <- as_data_frame(g_network, "edges") %>% as_tibble() %>% filter(!(is.na(bray))) %>% 
    mutate(bray07 = bray >= 0.7) %>% group_by(n_tax, bray07) %>% summarize(n_edges = dplyr::n()) %>% ungroup()
  p_edge <- df_plot %>%
    ggplot(aes(y = n_edges, x = n_tax, fill = bray07)) +
    geom_bar(stat = "identity", position = "fill") +
    geom_vline(xintercept = edge_thre_contra, size = 0.2, color = "black") +
    scale_fill_manual(values = c("grey70", hcl.colors(2,"Zissou 1")[2])) +
    labs(x = "Number of taxa", y = "Ratio of edges") #+
    #theme(strip.text.y = element_text(size = 13))
  return(p_edge)
}

vis_edge_fun_sharing_anti_single <- function(df_edge_og, shared_profile) {
	df_plot <- df_edge_og %>%
		filter(edge_type == "Anticorrelating", !is.na(.data[[shared_profile]])) %>%
		rename(share_tf = all_of(shared_profile)) %>%
		group_by(n_tax,share_tf) %>%
		summarize(n_edges = dplyr::n()) %>%
		mutate(fun_profile = shared_profile)
	df_plot <- df_plot %>%
		mutate(
			fun_profile_short = map(fun_profile, ~ 
				if (.x == "share_cog_fun_cat") {
					"COG Category"
				} else if (.x == "share_cog_fun_path") {
					"COG Pathway"
				} else if (.x == "share_ko_path_num_specific") {
					"KEGG\nPathway"
				} else if (.x == "share_ko_rn") {
					"KEGG\nReaction"
				} else if (.x == "share_ko_ec1") {
					"EC\n1 digit"
				} else if (.x == "share_ko_ec2") {
					"EC\n2 digits"
				} else if (.x == "share_ko_ec3") {
					"EC\n3 digits"
				} else if (.x == "share_ko_ec") {
					"EC number"
				}
			)
		) %>%
		unnest(c(fun_profile_short))
	df_plot$fun_profile_short <- factor(df_plot$fun_profile_short, levels = c("COG Category","COG Pathway","KEGG\nPathway","KEGG\nReaction","EC\n1 digit","EC\n2 digits","EC\n3 digits","EC number"))
	p_edge <- df_plot %>%
		ggplot(aes(y = n_edges, x = n_tax, fill = share_tf)) +
		geom_bar(stat = "identity", position = "fill") +
		geom_vline(xintercept = edge_thre_contra, size = 0.2, color = "black") +
		scale_fill_manual(values = hcl.colors(2,"TealRose")) +
		labs(x = "Number of taxa", y = "Ratio of edges", fill = "Sharing\nfunction") +
		theme(strip.text.y = element_text(size = 13))
	return(p_edge)
}

vis_edge_fun_sharing_binary <- function(df_edge_og, shared_profile) {
	n_tax_max <- max(df_edge_og$n_tax)
	p_edge <- df_edge_og %>%
		filter(!is.na(.data[[shared_profile]])) %>%
		mutate(edge_type_short = sub("Correlating", "Corr.", sub("Anticorrelating", "Anti.", edge_type))) %>%
		group_by(n_tax,edge_type_short,.data[[shared_profile]]) %>%
		summarize(n_edges = dplyr::n()) %>%
		mutate(n_tax_bin = cut(n_tax, c(1,8,n_tax_max), right = FALSE, include.lowest = TRUE)) %>%
		ggplot(aes(y = n_edges, x = n_tax_bin, fill = .data[[shared_profile]])) +
		geom_bar(stat = "identity", position = "fill") +
		geom_vline(xintercept = edge_thre_contra, size = 0.2, color = "black") +
		facet_grid(edge_type_short ~ .) +
		scale_fill_manual(values = hcl.colors(2,"TealRose")) +
		labs(x = "Number of taxa", y = "Ratio of edges", fill = "Sharing\nfunction") +
		theme(strip.text.y = element_text(size = 13))
	return(p_edge)
}

fun_color_df_generator <- function (g_network, fun_profile) {
	V(g_network)$degree <- igraph::degree(g_network)
	df_fun_color <- g_network %>%
		as_data_frame(what = "vertices") %>%
		as_tibble() %>%
		mutate(
			fun_single = strsplit(.data[[fun_profile]], ";")
		) %>%
		unnest(fun_single) %>%
		group_by(fun_single) %>%
		summarize(degree_weighted_fun = sum(degree)) %>%
		ungroup() %>%
		arrange(desc(degree_weighted_fun)) %>%
		filter(!is.na(fun_single)) %>%
		mutate(fun_color = fun_single) %>%
		head(11)
	df_fun_color$fun_color <- factor(df_fun_color$fun_color, levels = df_fun_color$fun_color)
	return(df_fun_color)
}

vis_stats_category_div <- function(g_contra_0, g_cooccu_0, fun_profile) {
	n_tax_max <- max(igraph::E(g_cooccu_0)$n_tax)
	df_contra_0_vertices <- g_contra_0 %>%
		as_data_frame(what = "vertices") %>%
		as_tibble()
	df_contra_0 <- g_contra_0 %>%
		as_data_frame(what = "edges") %>%
		as_tibble() %>%
		#filter(.data[[share_profile]] == TRUE) %>%
		mutate(n_tax_bin = cut(n_tax, c(1,8,n_tax_max), right = FALSE, include.lowest = TRUE)) %>%
		select(from,to,edge_type,n_tax_bin) %>%
		pivot_longer(cols = c(from,to), names_to = "vertex_group", values_to = "name") %>%
		select(-vertex_group) %>%
		group_by(edge_type,n_tax_bin,name) %>%
		summarize(n_edges = dplyr::n()) %>%
		ungroup() %>%
		left_join(
			df_contra_0_vertices,
			by = "name"
		) %>%
		mutate(
			fun_single = strsplit(.data[[fun_profile]], ";")
		) %>%
		unnest(fun_single) %>%
		group_by(edge_type,n_tax_bin,fun_single) %>%
		summarize(edges_weighted_fun = sum(n_edges)) %>%
		ungroup()
	df_cooccu_0_vertices <- g_cooccu_0 %>%
		as_data_frame(what = "vertices") %>%
		as_tibble()
	df_cooccu_0 <- g_cooccu_0 %>%
		as_data_frame(what = "edges") %>%
		as_tibble() %>%
		#filter(.data[[share_profile]] == TRUE) %>%
		mutate(n_tax_bin = cut(n_tax, c(1,8,n_tax_max), right = FALSE, include.lowest = TRUE)) %>%
		select(from,to,edge_type,n_tax_bin) %>%
		pivot_longer(cols = c(from,to), names_to = "vertex_group", values_to = "name") %>%
		select(-vertex_group) %>%
		group_by(edge_type,n_tax_bin,name) %>%
		summarize(n_edges = dplyr::n()) %>%
		ungroup() %>%
		left_join(
			df_cooccu_0_vertices,
			by = "name"
		) %>%
		mutate(
			fun_single = strsplit(.data[[fun_profile]], ";")
		) %>%
		unnest(fun_single) %>%
		group_by(edge_type,n_tax_bin,fun_single) %>%
		summarize(edges_weighted_fun = sum(n_edges)) %>%
		ungroup()
	df_fun_color <- df_contra_0 %>%
		filter(
			!is.na(fun_single),
			n_tax_bin == paste0("[", edge_thre_contra, ",", n_tax_max, "]")
		) %>%
		arrange(desc(edges_weighted_fun)) %>%
		head(11) %>%
		select(fun_single) %>%
		mutate(fun_color = fun_single)
	df_fun_color$fun_color <- factor(df_fun_color$fun_color, levels = c(df_fun_color$fun_color, "Others"))
	df_network_summary <- rbind(df_contra_0, df_cooccu_0) %>%
		filter(!is.na(fun_single)) %>%
		left_join(df_fun_color, by = "fun_single")
	df_network_summary$fun_color[is.na(df_network_summary$fun_color)] <- "Others"
	p_edge <- df_network_summary %>%
		ggplot(aes(x = n_tax_bin, y = edges_weighted_fun, fill = fun_color)) +
		geom_bar(stat = "identity", position = "fill") +
		facet_grid(edge_type ~ .) +
		scale_fill_manual(values = cbp12, na.value = "black") +
		scale_color_manual(values = cbp12, na.value = "black") +
		labs(x = "Number of taxa", y = "Degree-weighted counts")
	return(p_edge)
}

vis_stats_category_div_cog_fun_cat<- function(g_contra_0, g_contra_8) {
	V(g_contra_0)$degree <- igraph::degree(g_contra_0)
	V(g_contra_8)$degree <- igraph::degree(g_contra_8)
	df_network_summary <- rbind(
		g_contra_0 %>%
			as_data_frame(what = "vertices") %>%
			as_tibble() %>%
			select(name,cog_fun_cat_text,degree) %>%
			mutate(
				fun_single = strsplit(cog_fun_cat_text, ";")
			) %>%
			unnest(fun_single) %>%
			group_by(fun_single) %>%
			summarize(`Node count` = dplyr::n(), `Degree\nweighted` = sum(degree)) %>%
			mutate(Network = "Whole") %>%
			pivot_longer(cols = c(`Node count`, `Degree\nweighted`), names_to = "count_type", values_to = "count"),
		g_contra_8 %>%
			as_data_frame(what = "vertices") %>%
			as_tibble() %>%
			select(name,cog_fun_cat_text,degree) %>%
			mutate(
				fun_single = strsplit(cog_fun_cat_text, ";")
			) %>%
			unnest(fun_single) %>%
			group_by(fun_single) %>%
			summarize(`Node count` = dplyr::n(), `Degree\nweighted` = sum(degree)) %>%
			mutate(Network = "≥8 taxa") %>%
			pivot_longer(cols = c(`Node count`, `Degree\nweighted`), names_to = "count_type", values_to = "count")
	)
	df_fun_color <- df_network_summary %>%
		filter(Network == "≥8 taxa", count_type == "Degree\nweighted") %>%
		arrange(desc(count)) %>%
		head(11) %>%
		select(fun_single) %>%
		mutate(fun_color = fun_single)
	df_fun_color$fun_color <- factor(df_fun_color$fun_color, levels = c(df_fun_color$fun_color, "Others"))
	df_network_summary <- df_network_summary %>%
		left_join(df_fun_color, by = "fun_single")
	df_network_summary$fun_color[is.na(df_network_summary$fun_color)] <- "Others"
	df_network_summary$Network <- factor(df_network_summary$Network, levels = c("Whole","≥8 taxa"))
	df_network_summary$count_type <- factor(df_network_summary$count_type, levels = c("Node count","Degree\nweighted"))
	p_edge <- df_network_summary %>%
		ggplot(aes(x = Network, y = count, fill = fun_color)) +
		geom_bar(stat = "identity", position = "fill", color = "black") +
		facet_grid(count_type ~ .) +
		scale_fill_manual(values = cbp12, na.value = "black") +
		labs(y = "Ratio", fill = "COG Category") #+
		#theme(axis.text.x = element_text(angle = 0))
	return(p_edge)
}

node_degree_table <- function (g_graph) {
	V(g_graph)$degree <- igraph::degree(g_graph)
	df_node_degree <- g_graph %>%
		as_data_frame(what = "vertices") %>%
		as_tibble() %>%
		select(name,degree) %>%
		arrange(desc(degree))
	return(df_node_degree)
}

nodes_function_ranking <- function (g_graph, fun_profile) {
	V(g_graph)$degree <- igraph::degree(g_graph)
	if (fun_profile == "ko_path_num") {
		df_fun_title <- read_ko_map_title()
		df_fun_rank <- g_graph %>%
			as_data_frame(what = "vertices") %>%
			as_tibble() %>%
			select(.data[[fun_profile]],degree) %>%
			mutate(fun_single = strsplit(.data[[fun_profile]], ";")) %>%
			unnest(c(fun_single)) %>%
			left_join(
				df_fun_title,
				by = "fun_single"
			) %>%
			mutate(fun_title = paste(fun_single, fun_single_title, sep = ", ")) %>%
			filter(!is.na(fun_single)) %>%
			group_by(fun_single_higher_1, fun_title) %>%
			summarize(
				n_nodes = dplyr::n(),
				degree_weighted = sum(degree)
			) %>%
			ungroup() %>%
			arrange(desc(n_nodes))
		df_fun_rank$fun_title <- factor(df_fun_rank$fun_title, levels = rev(df_fun_rank$fun_title))
		df_fun_rank$fun_single_higher_1 <- factor(df_fun_rank$fun_single_higher_1, levels = unique(df_fun_rank$fun_single_higher_1))
	}
	return(df_fun_rank)
}

network_cog_graph2df_node_fun_cat <- function (cog_fun_name) {
	cog_fun <- "cog_fun_cat"
	g_network_tmp <- g_network
	set.seed(1)
	vertex_color_ogs <- g_network_tmp %>%
		as_data_frame(what = "vertices") %>%
		as_tibble() %>%
		rename(profile_category = .data[[cog_fun]]) %>%
		select(name,profile_category) %>%
		mutate(vertex_color = strsplit(profile_category, "")) %>%
		unnest(c(vertex_color)) %>%
		filter(vertex_color == cog_fun_name) %>%
		select(name) %>%
		pull()
	V(g_network_tmp)[vertex_color_ogs]$vertex_color <- cog_fun_name
	g_network_tmp %>%
		ggnetwork(layout = igraph::layout_with_fr(g_network_tmp)) %>%
		as_tibble() %>%
		mutate(
			x_same = x == xend,
			y_same = y == yend,
			v_e = map2(x_same, y_same, ~ if (.x == TRUE & .y == TRUE) {"vertex"} else {"edge"}),
			v_shape = map(.data[[cog_fun]], ~ if (length(grep("(R|S)", .x)) >= 1) {"23"} else {"21"}),
		) %>%
		unnest(c(v_e, v_shape)) %>%
		mutate(profile_category = cog_fun_name)
}

network_cog_graph2df_node_fun_path <- function (cog_fun_name) {
	cog_fun <- "cog_fun_path"
	g_network_tmp <- g_network
	set.seed(1)
	vertex_color_ogs <- g_network_tmp %>%
		as_data_frame(what = "vertices") %>%
		as_tibble() %>%
		rename(profile_category = .data[[cog_fun]]) %>%
		select(name,profile_category) %>%
		mutate(vertex_color = strsplit(profile_category, ";")) %>%
		unnest(c(vertex_color)) %>%
		filter(vertex_color == cog_fun_name) %>%
		select(name) %>%
		pull()
	V(g_network_tmp)[vertex_color_ogs]$vertex_color <- cog_fun_name
	g_network_tmp %>%
		ggnetwork(layout = igraph::layout_with_fr(g_network_tmp)) %>%
		as_tibble() %>%
		mutate(
			x_same = x == xend,
			y_same = y == yend,
			v_e = map2(x_same, y_same, ~ if (.x == TRUE & .y == TRUE) {"vertex"} else {"edge"}),
			v_shape = map(.data[[cog_fun]], ~ if (length(grep("(R|S)", .x)) >= 1) {"23"} else {"21"}),
		) %>%
		unnest(c(v_e, v_shape)) %>%
		mutate(profile_category = cog_fun_name)
}

contrapositive_network_cog_node_fun <- function (cog_fun, top_i) {
	profile_A <- paste0(cog_fun, "_A")
	profile_B <- paste0(cog_fun, "_B")
	df_og_replacement_cog_ar <- load_og_replacement("cog", "archaea")
	df_og_replacement_cog_bac <- load_og_replacement("cog", "bacteria")
	df_og_replacement_cog_prok <- og_replacement_prokaryotes(df_og_replacement_cog_ar, df_og_replacement_cog_bac)
	df_network <- read_network_edge_stats("cog", edge_thre_contra, "contrapositive") %>%
		mutate(edge_id = paste(og_id_A, og_id_B, sep = "_")) %>%
		left_join(
			df_og_replacement_cog_prok %>% mutate(edge_id = paste(og_id_A, og_id_B, sep = "_")) %>% select(-og_id_A, -og_id_B),
			by = "edge_id"
		) %>%
		filter(n_replace_prok >= replace_thre)
	df_def <- read_cog2020_def_fixed() %>%
		right_join(
			tibble(cog_id = unique(c(df_network$og_id_A, df_network$og_id_B))),
			by = "cog_id"
		) %>%
		rename(og_id = "cog_id")
	#na_nadh_dehydro <- df_def$cog_fun_path == "Na+-translocating NADH dehydrogenase"
	#na_nadh_dehydro[is.na(na_nadh_dehydro)] <- FALSE
	#df_def$cog_fun_path[na_nadh_dehydro] <- "Na+-translocating\nNADH dehydrogenase"
	g_network <<- df_network %>% graph_from_data_frame(directed = FALSE, vertices = df_def)
	colnames(df_def) <- sub("$", "_A", colnames(df_def))
	df_network <- df_network %>% left_join(df_def, by = "og_id_A")
	colnames(df_def) <- sub("_A$", "_B", colnames(df_def))
	df_network <- df_network %>% left_join(df_def, by = "og_id_B")
	if (cog_fun == "cog_fun_cat") {
		df_network_cog_fun <- df_network %>%
			select(n_tax, edge_type, .data[[profile_A]], .data[[profile_B]]) %>%
			pivot_longer(cols = c(-n_tax,-edge_type), names_to = "profile", values_to = "profile_categories") %>%
			mutate(profile_category = strsplit(profile_categories, "")) %>%
			unnest(c(profile_category)) %>%
			filter(!is.na(profile_category)) %>%
			group_by(profile_category) %>%
			summarize(n_edges = dplyr::n()) %>%
			arrange(desc(n_edges)) %>%
			head(top_i)
		df_g_network <- do.call(rbind, mclapply(df_network_cog_fun$profile_category, network_cog_graph2df_node_fun_cat, mc.cores = 1))
		df_g_network$profile_category <- factor(df_g_network$profile_category, levels = df_network_cog_fun$profile_category)
		df_g_network$vertex_color <- factor(df_g_network$vertex_color, levels = df_network_cog_fun$profile_category)
		df_g_network$n_replace_prok[df_g_network$n_replace_prok >= 100] <- 100
		p_network <- ggplot() +
			geom_segment(
				aes(x = x, y = y, xend = xend, yend = yend, size = n_tax),
				color = "grey70",
				data = df_g_network %>% filter(v_e == "edge")
			) +
			geom_point(
				aes(x = x, y = y, fill = vertex_color),
				color = "white",
				size = 2,
				shape = 21,
				data = df_g_network %>% filter(v_e == "vertex", !is.na(vertex_color))
			) +
			scale_fill_manual(values = cbp12, na.value = "black") +
			facet_wrap(. ~ profile_category, ncol = 4) +
			#coord_fixed(ratio = 1) +
			theme_blank() +
			theme(
				text = element_text(family = "Arial", size = 12),
				plot.title = element_text(hjust = 0.5, size = 13),
				plot.background = element_rect(fill = "white",color = NA),
				plot.tag = element_text(face = "bold", size = 18),
				strip.text.x = element_text(size = 13),
				strip.text.y = element_text(size = 13),
				legend.position = "none"
			)
	} else if (cog_fun == "cog_fun_path") {
		df_network_cog_fun <- df_network %>%
			select(n_tax, edge_type, .data[[profile_A]], .data[[profile_B]]) %>%
			pivot_longer(cols = c(-n_tax,-edge_type), names_to = "profile", values_to = "profile_categories") %>%
			mutate(profile_category = strsplit(profile_categories, ";")) %>%
			unnest(c(profile_category)) %>%
			filter(!is.na(profile_category)) %>%
			group_by(profile_category) %>%
			summarize(n_edges = dplyr::n()) %>%
			arrange(desc(n_edges)) %>%
			head(top_i)
		df_g_network <- do.call(rbind, mclapply(df_network_cog_fun$profile_category, network_cog_graph2df_node_fun_path, mc.cores = 1))
		df_g_network$profile_category <- factor(df_g_network$profile_category, levels = df_network_cog_fun$profile_category)
		df_g_network$vertex_color <- factor(df_g_network$vertex_color, levels = df_network_cog_fun$profile_category)
		df_g_network$n_replace_prok[df_g_network$n_replace_prok >= 100] <- 100
		p_network <- ggplot() +
			geom_segment(
				aes(x = x, y = y, xend = xend, yend = yend, size = n_tax),
				color = "grey70",
				data = df_g_network %>% filter(v_e == "edge")
			) +
			geom_point(
				aes(x = x, y = y, fill = vertex_color),
				color = "white",
				size = 3,
				shape = 21,
				data = df_g_network %>% filter(v_e == "vertex", !is.na(vertex_color))
			) +
			scale_fill_manual(values = cbp12, na.value = "black") +
			guides(size = "none") +
			labs(fill = "COG Pathway") +
			#facet_wrap(. ~ profile_category, ncol = 4) +
			#coord_fixed(ratio = 1) +
			theme_blank() +
			theme(
				text = element_text(family = "Arial", size = 12),
				plot.title = element_text(hjust = 0.5, size = 13),
				plot.background = element_rect(fill = "white",color = NA),
				plot.tag = element_text(face = "bold", size = 18),
				strip.text.x = element_text(size = 13),
				strip.text.y = element_text(size = 13)
				#legend.position = "none"
			)
	}
	return(p_network)
}

contrapositive_network_node_cog_fun_cat_exclusive <- function (top_i) {
	cog_fun <- "cog_fun_cat"
	profile_A <- paste0(cog_fun, "_A")
	profile_B <- paste0(cog_fun, "_B")
	df_og_replacement_cog_ar <- load_og_replacement("cog", "archaea")
	df_og_replacement_cog_bac <- load_og_replacement("cog", "bacteria")
	df_og_replacement_cog_prok <- og_replacement_prokaryotes(df_og_replacement_cog_ar, df_og_replacement_cog_bac)
	df_network <- read_network_edge_stats("cog", edge_thre_contra, "contrapositive") %>%
		mutate(edge_id = paste(og_id_A, og_id_B, sep = "_")) %>%
		left_join(
			df_og_replacement_cog_prok %>% mutate(edge_id = paste(og_id_A, og_id_B, sep = "_")) %>% select(-og_id_A, -og_id_B),
			by = "edge_id"
		) %>%
		filter(n_replace_prok >= replace_thre)
	df_def <- read_cog2020_def_fixed() %>%
		right_join(
			tibble(cog_id = unique(c(df_network$og_id_A, df_network$og_id_B))),
			by = "cog_id"
		) %>%
		rename(og_id = "cog_id")
	g_network <<- df_network %>% graph_from_data_frame(directed = FALSE, vertices = df_def)
	colnames(df_def) <- sub("$", "_A", colnames(df_def))
	df_network <- df_network %>% left_join(df_def, by = "og_id_A")
	colnames(df_def) <- sub("_A$", "_B", colnames(df_def))
	df_network <- df_network %>% left_join(df_def, by = "og_id_B")
	df_network_cog_fun <- df_network %>%
		select(n_tax, edge_type, .data[[profile_A]], .data[[profile_B]]) %>%
		pivot_longer(cols = c(-n_tax,-edge_type), names_to = "profile", values_to = "profile_categories") %>%
		mutate(cog_fun_cat = strsplit(profile_categories, "")) %>%
		unnest(c(cog_fun_cat)) %>%
		filter(!is.na(cog_fun_cat)) %>%
		group_by(cog_fun_cat) %>%
		summarize(n_edges = dplyr::n()) %>%
		arrange(desc(n_edges)) %>%
		head(top_i) %>%
		select(-n_edges) %>%
		mutate(vertex_color = cog_fun_cat)
	set.seed(1)
	df_g_network <- g_network %>%
		ggnetwork(layout = igraph::layout_with_fr(g_network)) %>%
		as_tibble() %>%
		mutate(
			x_same = x == xend,
			y_same = y == yend,
			v_e = map2(x_same, y_same, ~ if (.x == TRUE & .y == TRUE) {"vertex"} else {"edge"}),
			v_shape = map(.data[[cog_fun]], ~ if (length(grep("(R|S)", .x)) >= 1) {"23"} else {"21"}),
		) %>%
		unnest(c(v_e, v_shape)) %>%
		left_join(df_network_cog_fun, by = "cog_fun_cat")
	df_g_network$vertex_color <- factor(df_g_network$vertex_color, levels = df_network_cog_fun$vertex_color)
	df_g_network$n_replace_prok[df_g_network$n_replace_prok >= 100] <- 100
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend, size = n_tax),
			color = "grey70",
			data = df_g_network %>% filter(v_e == "edge")
		) +
		geom_point(
			aes(x = x, y = y, fill = vertex_color),
			color = "white",
			size = 3,
			shape = 21,
			data = df_g_network %>% filter(v_e == "vertex", !is.na(vertex_color))
		) +
		scale_fill_manual(values = cbp12, na.value = "black") +
		guides(size = "none") +
		#coord_fixed(ratio = 1) +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

contrapositive_network_og_bipartite <- function (g_graph, OG, As, Bs, A_color, B_color, label) {
	og_As <- strsplit(As, ";") %>% unlist()
	og_Bs <- strsplit(Bs, ";") %>% unlist()
	g_network <- g_graph
	sub_network <- subgraph(g_network, c(og_As, og_Bs))
	V(sub_network)$type <- bipartite_mapping(sub_network)$type
	df_sub_network <- ggnetwork_fr_bipartite(sub_network)
	if (label == TRUE) {
		p_network <- ggplot() +
			geom_segment(
				aes(x = x, y = y, xend = xend, yend = yend, color = n_replace_prok),
				size = 1,
				data = df_sub_network %>% filter(v_e == "edge")
			) +
			geom_point(
				aes(x = x, y = y, fill = type),
				shape = 21,
				size = 4,
				color = "white",
				data = df_sub_network %>% filter(v_e == "vertex")
			) +
			geom_text_repel(
				aes(x = x, y = y, label = name_custom),
				min.segment.length = 0,
				segment.size = 0.1,
				segment.linetype = 1,
				segment.color = "black",
				size = 2,
				data = df_sub_network %>% filter(v_e == "vertex")
			) +
			scale_fill_manual(values = c(A_color, B_color), na.value = "grey50") +
			scale_color_viridis_c(option = "plasma", begin = 0, end = 0.9) +
			xlim(-0.15, 1.15) +
			ylim(-0.4, 1.4) +
			theme_blank() +
			theme(
				plot.title = element_text(face = "bold", hjust = 0, size = 18),
				plot.background = element_rect(fill = "white",color = NA),
				plot.tag = element_text(face = "bold", size = 18),
				strip.text.x = element_text(size = 13),
				strip.text.y = element_text(size = 13),
				legend.position = "none"
			)
	} else if (label == FALSE) {
		p_network <- ggplot() +
			geom_segment(
				aes(x = x, y = y, xend = xend, yend = yend, color = n_replace_prok),
				size = 1,
				data = df_sub_network %>% filter(v_e == "edge")
			) +
			geom_point(
				aes(x = x, y = y, fill = type),
				shape = 21,
				size = 4,
				color = "white",
				data = df_sub_network %>% filter(v_e == "vertex")
			) +
			#geom_text_repel(
			#	aes(x = x, y = y, label = name_custom),
			#	min.segment.length = 0,
			#	segment.size = 0.1,
			#	segment.linetype = 1,
			#	segment.color = "black",
			#	size = 2,
			#	data = df_sub_network %>% filter(v_e == "vertex")
			#) +
			scale_fill_manual(values = c(A_color, B_color), na.value = "grey50") +
			scale_color_viridis_c(option = "plasma", begin = 0, end = 0.9) +
			xlim(-0.15, 1.15) +
			ylim(-0.4, 1.4) +
			theme_blank() +
			theme(
				plot.title = element_text(face = "bold", hjust = 0, size = 18),
				plot.background = element_rect(fill = "white",color = NA),
				plot.tag = element_text(face = "bold", size = 18),
				strip.text.x = element_text(size = 13),
				strip.text.y = element_text(size = 13),
				legend.position = "none"
			)
	}
	return(p_network)
}

visualize_og_heatmap <- function (OG, dom, As, Bs, tax_level, A_color, B_color) {
	if (OG == "cog" & dom == "archaea") {
		og_path <- "data/pastml_cog2020_ar/"
		tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Archaea"
		node_plot_size <- 1.4
	} else if (OG == "cog" & dom == "bacteria") {
		og_path <- "data/pastml_cog2020_bac/"
		tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Bacteria"
		node_plot_size <- 0.8
	} else if (OG == "kegg" & dom == "archaea") {
		og_path <- "data/pastml_ko_ar/"
		tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Archaea"
		node_plot_size <- 1.4
	} else if (OG == "kegg" & dom == "bacteria") {
		og_path <- "data/pastml_ko_bac/"
		tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Bacteria"
		node_plot_size <- 0.8
	}
	OGs_A <- strsplit(As, ";") %>% unlist()
	OGs_B <- strsplit(Bs, ";") %>% unlist()
	df_og_a <- tibble()
	df_og_b <- tibble()
	for (i in OGs_A) {
		df_og_a <- rbind(
			df_og_a,
			read.csv(paste0(og_path, all_of(i), ".txt"), header = TRUE, sep = "\t", quote = "") %>%
				as_tibble() %>% rename(label = "node", TF = all_of(i)) %>% mutate(og_id = all_of(i))
		)
	}
	for (i in OGs_B) {
		df_og_b <- rbind(
			df_og_b,
			read.csv(paste0(og_path, all_of(i), ".txt"), header = TRUE, sep = "\t", quote = "") %>%
				as_tibble() %>% rename(label = "node", TF = all_of(i)) %>% mutate(og_id = all_of(i))
		)
	}
	gtdb_node_labeled_tree <- read.tree(tree_path)
	df_genomes <- read_genomes_analyzed() %>% select(accession, Domain:Species) %>% rename(label = "accession")
	df_og_a_pivot <- df_og_a %>% mutate(presence = 1) %>%
		pivot_wider(names_from = TF, values_from = presence, values_fill = 0) %>% rename(og_TRUE = "TRUE", og_FALSE = "FALSE") %>% mutate(OG_group = "A")
	df_og_b_pivot <- df_og_b %>% mutate(presence = 1) %>%
		pivot_wider(names_from = TF, values_from = presence, values_fill = 0) %>% rename(og_TRUE = "TRUE", og_FALSE = "FALSE") %>% mutate(OG_group = "B")
	df_tree_mapped <- gtdb_node_labeled_tree %>% as_tibble() %>%
		left_join(rbind(df_og_a_pivot, df_og_b_pivot), by = "label") %>%
		mutate(
			og_presence = map2(og_TRUE, og_FALSE, ~
				if (.x == 1 & .y == 0){
					1
				} else if (.x == 1 & .y == 1) {
					0.5
				} else if (.x == 0 & .y == 1) {
					0
				}
			)
		) %>%
		unnest(c(og_presence))
	p_tree <- ggtree(gtdb_node_labeled_tree, branch.length = 'none', layout = 'circular', color = "gray")
	df_tip_data <- df_tree_mapped[which(df_tree_mapped$label %in% gtdb_node_labeled_tree$tip.label), c("label", "og_id", "og_presence", "OG_group")] %>%
		pivot_longer(c(-label,-og_id,-OG_group), names_to = "ogs_tip", values_to = "pastml_tip") %>%
		mutate(
			og_presence = pastml_tip >= 1,
			OGs_tip = og_id
		) %>%
		left_join(p_tree$data, by = "label") %>%
		rename(ID = "label")
	df_tip_data$OGs_tip <- factor(df_tip_data$OGs_tip, levels = c(OGs_A, OGs_B))
	df_tax <- df_genomes[which(df_genomes$label %in% p_tree$data[p_tree$data$isTip,]$label), c("label", tax_level)]
	df_tax_top <- df_tax %>% group_by(.data[[tax_level]]) %>% summarize(count = dplyr::n()) %>% arrange(desc(count)) %>%
		head(11) %>% select(-count) %>% mutate(tax_color = .data[[tax_level]])
	df_tax <- df_tax %>% left_join(df_tax_top, by = tax_level)
	df_tax$tax_color[is.na(df_tax$tax_color)] <- "Others"
	df_tax$tax_color <- factor(df_tax$tax_color, levels = c(df_tax_top$tax_color, "Others"))
	df_tax_data <- left_join(p_tree$data, df_tax, by = "label") %>% filter(isTip == TRUE)
	p_tree_plotted <- p_tree +
		geom_fruit(
			data = df_tip_data %>% filter(og_presence == TRUE) %>% select(ID, OGs_tip, OG_group),
			geom = geom_tile,
			mapping = aes(y = ID, x = OGs_tip, fill = OG_group),
			color = NA,
			pwidth = 0.07 * length(c(OGs_A, OGs_B)),
			offset = 0.05
		) +
		guides(fill = "none") +
		scale_fill_manual(values = c(A_color, B_color)) +
		new_scale_fill() +
		geom_fruit(
			data = df_tax_data %>% rename(ID = "label") %>% select(ID, tax_color) %>% mutate(x_pos = 1),
			geom = geom_star,
			mapping = aes(y = ID, x = x_pos, fill = tax_color),
			starshape = 15,
			color = NA,
			alpha = 0.8,
			pwidth = 0.08
		) +
		scale_fill_manual(values = cbp12, na.value = cbp12[12]) +
		labs(fill = tax_level) +
		layout_rectangular() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.tag = element_text(face = "bold", size = 18),
			plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
		)
		return(p_tree_plotted)
}

visualize_ogs_tree_binary <- function (OG, dom, As, Bs, tax_level, og_pwidth, og_offset, tax_pwidth, tax_offset) {
	if (OG == "cog" & dom == "archaea") {
		og_path <- "data/pastml_cog2020_ar/"
		tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Archaea"
		node_plot_size <- 1.4
	} else if (OG == "cog" & dom == "bacteria") {
		og_path <- "data/pastml_cog2020_bac/"
		tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Bacteria"
		node_plot_size <- 0.8
	} else if (OG == "kegg" & dom == "archaea") {
		og_path <- "data/pastml_ko_ar/"
		tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Archaea"
		node_plot_size <- 1.4
	} else if (OG == "kegg" & dom == "bacteria") {
		og_path <- "data/pastml_ko_bac/"
		tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Bacteria"
		node_plot_size <- 0.8
	}
	df_habitat_pref <-
		read.csv("data/tsv/gtdb_ssu_prokatlas_habitat_pref.tsv", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
		as_tibble()
	OGs_A <- strsplit(As, ";") %>% unlist()
	OGs_B <- strsplit(Bs, ";") %>% unlist()
	df_og_a <- tibble()
	df_og_b <- tibble()
	for (i in OGs_A) {
		df_og_a <- rbind(
			df_og_a,
			read.csv(paste0(og_path, all_of(i), ".txt"), header = TRUE, sep = "\t", quote = "") %>%
				as_tibble() %>% rename(label = "node", TF = all_of(i)) %>% mutate(og_id = all_of(i))
		)
	}
	for (i in OGs_B) {
		df_og_b <- rbind(
			df_og_b,
			read.csv(paste0(og_path, all_of(i), ".txt"), header = TRUE, sep = "\t", quote = "") %>%
				as_tibble() %>% rename(label = "node", TF = all_of(i)) %>% mutate(og_id = all_of(i))
		)
	}
	gtdb_node_labeled_tree <- read.tree(tree_path)
	df_genomes <- read_genomes_analyzed() %>% select(accession, Domain:Species) %>% rename(label = "accession")
	df_og_a_pivot <- df_og_a %>% mutate(presence = 1) %>%
		pivot_wider(names_from = TF, values_from = presence, values_fill = 0) %>% rename(og_TRUE = "TRUE", og_FALSE = "FALSE") %>% mutate(OG_group = "A OGs")
	df_og_b_pivot <- df_og_b %>% mutate(presence = 1) %>%
		pivot_wider(names_from = TF, values_from = presence, values_fill = 0) %>% rename(og_TRUE = "TRUE", og_FALSE = "FALSE") %>% mutate(OG_group = "B OGs")
	df_tree_mapped <- gtdb_node_labeled_tree %>% as_tibble() %>%
		left_join(rbind(df_og_a_pivot, df_og_b_pivot), by = "label") %>%
		mutate(
			og_presence = map2(og_TRUE, og_FALSE, ~
				if (.x == 1 & .y == 0){
					1
				} else if (.x == 1 & .y == 1) {
					0.5
				} else if (.x == 0 & .y == 1) {
					0
				}
			)
		) %>%
		unnest(c(og_presence))
	p_tree <- ggtree(gtdb_node_labeled_tree, branch.length = 'none', layout = 'circular', color = "gray")
	df_tip_data <- df_tree_mapped[which(df_tree_mapped$label %in% gtdb_node_labeled_tree$tip.label), c("label", "og_id", "og_presence", "OG_group")] %>%
		pivot_longer(c(-label,-og_id,-OG_group), names_to = "ogs_tip", values_to = "pastml_tip") %>%
		mutate(
			og_presence = pastml_tip >= 1,
			OGs_tip = og_id
		) %>%
		left_join(p_tree$data, by = "label") %>%
		rename(ID = "label")
	df_tip_data$OGs_tip <- factor(df_tip_data$OGs_tip, levels = c(OGs_A, OGs_B))
	df_tax <- df_genomes[which(df_genomes$label %in% p_tree$data[p_tree$data$isTip,]$label), c("label", tax_level)]
	df_tax_top <- df_tax %>% group_by(.data[[tax_level]]) %>% summarize(count = dplyr::n()) %>% arrange(desc(count)) %>%
		head(11) %>% select(-count) %>% mutate(tax_color = .data[[tax_level]])
	df_tax <- df_tax %>% left_join(df_tax_top, by = tax_level)
	df_tax$tax_color[is.na(df_tax$tax_color)] <- "Others"
	df_tax$tax_color <- factor(df_tax$tax_color, levels = c(df_tax_top$tax_color, "Others"))
	df_tax_data <- left_join(p_tree$data, df_tax, by = "label") %>% filter(isTip == TRUE)
	p_tree_plotted <- p_tree +
		geom_fruit(
			data = df_tip_data %>% filter(og_presence == TRUE) %>% select(ID,OGs_tip,OG_group),
			geom = geom_tile,
			mapping = aes(y = ID, x = OGs_tip, fill = OG_group),
			color = NA,
			pwidth = og_pwidth,
			offset = og_offset
		) +
		labs(fill = "OGs") +
		scale_fill_manual(values = hcl.colors(2,"Blue-Red 2")) +
		new_scale_fill() +
		geom_fruit(
			data = rbind(
				df_tax_data %>% rename(ID = "label") %>% select(ID, tax_color) %>% mutate(x_pos = 1),
				df_tax_data %>% rename(ID = "label") %>% select(ID, tax_color) %>% mutate(x_pos = 2)
			),
			geom = geom_tile,
			mapping = aes(y = ID, x = x_pos, fill = tax_color),
			color = NA,
			pwidth = tax_pwidth,
			offset = tax_offset
		) +
		scale_fill_manual(values = cbp12, na.value = cbp12[12]) +
		layout_rectangular() +
		coord_flip() +
		scale_y_reverse() +
		labs(fill = tax_level) +
		theme(
			text = element_text(family = "Arial", size = 12),
			legend.position = "right",
			plot.title = element_text(hjust = 0, face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_tree_plotted)
}

generate_large_component_obj <- function (g_network) {
	compo_g_network <- components(g_network)
	large_compo_vertices <- compo_g_network$membership[compo_g_network$membership == which(compo_g_network$csize %in% max(compo_g_network$csize))] %>% names()
	g_large_component <- subgraph(g_network, large_compo_vertices)
	V(g_large_component)$betweenness <- betweenness(g_large_component)
	V(g_large_component)$degree <- igraph::degree(g_large_component)
	g_large_component
	return(g_large_component)
}

show_network_large_component_color <- function (g_large, node_color, show_label) {
	df_g_large <- ggnetwork_fr(g_large)
	df_label_show <- g_large %>%
		as_data_frame(what = "vertices") %>%
		as_tibble() %>%
		arrange(desc(.data[[node_color]])) %>%
		head(15) %>%
		mutate(
			label_show = map2(name, .data[[show_label]], ~ if (is.na(.y)) {.x} else {.y})
		) %>%
		unnest(c(label_show))
	df_g_large <- left_join(
		df_g_large,
		df_label_show %>% select(-.data[[node_color]]),
		by = "name"
	)
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend),
			size = 1,
			color = "grey70",
			data = df_g_large %>% filter(v_e == "edge")
		) +
		geom_point(
			aes(x = x, y = y, fill = .data[[node_color]]),
			color = "white",
			shape = 21,
			data = df_g_large %>% filter(v_e == "vertex", !is.na(.data[[node_color]]))
		) +
		geom_text_repel(
			aes(x = x, y = y, label = label_show, color = .data[[node_color]]),
			min.segment.length = 0,
			segment.size = 0.1,
			segment.linetype = 1,
			segment.color = "black",
			nudge_x = 1.02,
			direction = "y",
			hjust = 1,
			size = 2,
			data = df_g_large %>% filter(v_e == "vertex")
		) +
		xlim(0,1.2) +
		scale_fill_viridis_c(na.value = "grey70") +
		scale_color_viridis_c(na.value = "grey70") +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(list(p_network, df_label_show))
}

show_network_node_fun <- function (g_network, df_fun_color, fun_profile, show_label, ogs_a, ogs_b) {
	df_ogs <- rbind(
			tibble(
				name = strsplit(ogs_a, ";") %>% unlist(),
				label_group = "og_a"
			),
			tibble(
				name = strsplit(ogs_b, ";") %>% unlist(),
				label_group = "og_b"
			)
		)
	set.seed(1)
	if (fun_profile == "cog_fun_path") {
		df_g_network <- ggnetwork_fr(g_network) %>%
			rename(fun_single = .data[[fun_profile]])
	} else if (fun_profile == "ko_path_num") {
		df_fun_title <- read_ko_map_title() %>%
			rename(
				ko_path_num_single = fun_single,
				fun_single = fun_single_higher_1
			)
		df_g_network <- ggnetwork_fr(g_network) %>%
			as_tibble() %>%
			mutate(
				ko_path_num_single = strsplit(ko_path_num, ";")
			) %>%
			unnest(ko_path_num_single) %>%
			left_join(df_fun_title, by = "ko_path_num_single")
	}
	df_label_show <- g_network %>%
		as_data_frame(what = "vertices") %>%
		as_tibble() %>%
		mutate(
			label_show = map2(name, .data[[show_label]], ~ if (is.na(.y)) {.x} else {.y})
		) %>%
		unnest(c(label_show)) %>%
		right_join(
			df_ogs,
			by = "name"
		)
	df_g_network <- left_join(
			df_g_network,
			df_fun_color,
			by = "fun_single"
		) %>%
		left_join(
			df_label_show %>% select(name,label_group,label_show),
			by = "name"
		)
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend),
			size = 1,
			color = "grey70",
			data = df_g_network %>% filter(v_e == "edge")
		) +
		geom_point(
			aes(x = x, y = y, fill = fun_color),
			size = 1.5,
			color = "white",
			shape = 21,
			data = df_g_network %>% filter(v_e == "vertex", !is.na(fun_color))
		) +
		geom_text_repel(
			aes(x = x, y = y, label = label_show, color = fun_color),
			min.segment.length = 0,
			segment.size = 0.1,
			segment.linetype = 1,
			segment.color = "black",
			nudge_x = 1.02,
			direction = "y",
			hjust = 1,
			size = 3,
			data = df_g_network %>% filter(v_e == "vertex", !is.na(fun_color))
		) +
		xlim(0,1.2) +
		scale_fill_manual(values = cbp12) +
		scale_color_manual(values = cbp12, na.value = "grey70") +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(list(p_network, df_label_show))
}

pick_nodes_in_whole_network <- function (g_network, ogs, og_color, show_label) {
	OGs <- strsplit(ogs, ";") %>% unlist()
	set.seed(1)
	df_g_network <- ggnetwork_fr(g_network)
	df_label_show <- g_network %>%
		as_data_frame(what = "vertices") %>%
		as_tibble() %>%
		mutate(
			label_show = map2(name, .data[[show_label]], ~ if (is.na(.y)) {.x} else {.y})
		) %>%
		unnest(c(label_show)) %>%
		right_join(
			tibble(
				name = OGs
			),
			by = "name"
		)
	df_g_network <- left_join(
			df_g_network,
			df_label_show %>% select(name,label_show),
			by = "name"
		)
	df_g_network$label_show <- factor(df_g_network$label_show, levels = df_label_show$label_show)
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend),
			size = 1,
			color = "grey70",
			data = df_g_network %>% filter(v_e == "edge")
		) +
		geom_point(
			aes(x = x, y = y, fill = label_show),
			size = 1.5,
			color = "white",
			shape = 21,
			data = df_g_network %>% filter(v_e == "vertex", !is.na(label_show))
		) +
		geom_text_repel(
			aes(x = x, y = y, label = label_show, color = label_show),
			min.segment.length = 0,
			segment.size = 0.1,
			segment.linetype = 1,
			segment.color = "black",
			nudge_x = 1.02,
			direction = "y",
			hjust = 1,
			size = 3,
			data = df_g_network %>% filter(v_e == "vertex", !is.na(label_show))
		) +
		xlim(0,1.2) +
		scale_fill_manual(values = og_color) +
		scale_color_manual(values = og_color, na.value = "grey70") +
		theme_blank() +
		theme(
			legend.position = "none",
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

pick_nodes_in_network <- function (g_network, ogs, og_color) {
	OGs <- strsplit(ogs, ";") %>% unlist()
	g_subgraph <- subgraph(g_network, OGs)
	set.seed(1)
	df_g_network <- ggnetwork_fr(g_subgraph)
	df_g_network$name <- factor(df_g_network$name, levels = OGs)
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend, size = n_tax, color = n_replace_prok),
			data = df_g_network %>% filter(v_e == "edge")
		) +
		scale_color_viridis_c(
			option = "plasma",
			begin = 0,
			end = 0.9,
			breaks = c(min(igraph::E(g_network)$n_replace_prok), 30, 60, max(igraph::E(g_network)$n_replace_prok))
		) +
		scale_size_continuous(
			breaks = c(min(igraph::E(g_network)$n_tax), 100, 200, 300, max(igraph::E(g_network)$n_tax))
		) +
		geom_point(
			aes(x = x, y = y, fill = name),
			size = 1.5,
			color = "white",
			shape = 21,
			data = df_g_network %>% filter(v_e == "vertex")
		) +
		geom_text_repel(
			aes(x = x, y = y, label = name),
			min.segment.length = 0,
			segment.size = 0.1,
			segment.linetype = 1,
			segment.color = "black",
			#nudge_x = 1.02,
			#direction = "y",
			#hjust = 1,
			size = 3,
			data = df_g_network %>% filter(v_e == "vertex")
		) +
		scale_fill_manual(values = og_color) +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			legend.position = "right",
			plot.title = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

pick_nodes_in_whole_network_cog_fun_cat <- function (g_network, ogs) {
	df_color <- fun_color_df_generator(g_network, "cog_fun_cat_text")
	OGs <- ogs
	#g_subgraph <- subgraph(g_network, OGs)
	set.seed(1)
	df_g_network <- ggnetwork_fr(g_network)
	#df_g_network$name <- factor(df_g_network$name, levels = OGs)
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend),
			size = 0.8,
			color = "grey80",
			data = df_g_network %>% filter(v_e == "edge")
		) +
		#scale_color_viridis_c(
		#	option = "plasma",
		#	begin = 0,
		#	end = 0.9,
		#	breaks = c(min(igraph::E(g_network)$n_replace_prok), 30, 60, max(igraph::E(g_network)$n_replace_prok))
		#) +
		#scale_size_continuous(
		#	breaks = c(min(igraph::E(g_network)$n_tax), 100, 200, 300, max(igraph::E(g_network)$n_tax))
		#) +
		geom_point(
			aes(x = x, y = y, size = degree),
			color = "white",
			fill = "grey90",
			shape = 21,
			data = df_g_network %>% filter(v_e == "vertex", !(name %in% OGs)) %>%
				mutate(fun_single = cog_fun_cat_text) %>% left_join(df_color, by = "fun_single")
		) +
		geom_point(
			aes(x = x, y = y, fill = fun_color, size = degree),
			color = "white",
			shape = 21,
			data = df_g_network %>% filter(v_e == "vertex", name %in% OGs) %>%
				mutate(fun_single = cog_fun_cat_text) %>% left_join(df_color, by = "fun_single")
		) +
		scale_size(range = c(1.2, 3)) +
		#geom_text_repel(
		#	aes(x = x, y = y, label = name),
		#	min.segment.length = 0,
		#	segment.size = 0.1,
		#	segment.linetype = 1,
		#	segment.color = "black",
		#	#nudge_x = 1.02,
		#	#direction = "y",
		#	#hjust = 1,
		#	size = 3,
		#	data = df_g_network %>% filter(v_e == "vertex")
		#) +
		coord_equal() +
		scale_fill_manual(values = cbp12) +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			legend.position = "none",
			plot.title = element_text(hjust = 0, size = 18, face = "bold"),
			#plot.title = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

load_genome_og_tidy <- function (dom, OG) {
	if (dom == "bacteria" & OG == "cog") {
		df_genome_og_tidy <- fread("data/proteins_raw_gtdbtax_cog2020_tidy/proteins_raw_gtdbtax_cog2020_tidy_d__Bacteria.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("", "NULL")) %>% as_tibble() %>% rename(og_id = "cog_id")
	} else if (dom == "archaea" & OG == "cog") {
		df_genome_og_tidy <- fread("data/proteins_raw_gtdbtax_cog2020_tidy/proteins_raw_gtdbtax_cog2020_tidy_d__Archaea.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("", "NULL")) %>% as_tibble() %>% rename(og_id = "cog_id")
	} else if (dom == "bacteria" & OG == "kegg") {
		df_genome_og_tidy <- fread("data/proteins_raw_gtdbtax_ko_tidy/proteins_raw_gtdbtax_ko_tidy_d__Bacteria.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("", "NULL")) %>% as_tibble() %>% mutate(og_id = sub("^ko:", "", ko_id)) %>% select(-ko_id)
	} else if (dom == "archaea" & OG == "kegg") {
		df_genome_og_tidy <- fread("data/proteins_raw_gtdbtax_ko_tidy/proteins_raw_gtdbtax_ko_tidy_d__Archaea.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("", "NULL")) %>% as_tibble() %>% mutate(og_id = sub("^ko:", "", ko_id)) %>% select(-ko_id)
	}
	return(df_genome_og_tidy)
}

summary_genome_habitat <- function (df_genomes_bac_og, df_genomes_ar_og, df_habitat_pref, g_large) {
	df_env <- tibble(
		habitat = c("bovine_gut","chicken_gut","fish_gut","gut","human_gut","insect_gut","mouse_gut","pig_gut","rat_gut","termite_gut",
			"soil","rhizosphere",
			"marine","seawater","marine_sediment",
			"freshwater","freshwater_sediment","lake_water"),
		env = c("gut","gut","gut","gut","gut","gut","gut","gut","gut","gut",
			"soil","rhizosphere",
			"marine","marine","marine",
			"freshwater","freshwater","freshwater")
	)
	df_env$env <- factor(df_env$env, levels = c(unique(df_env$env), "Others"))
	df_genomes_prok_env <- rbind(
			df_genomes_bac_og,
			df_genomes_ar_og
		) %>%
		inner_join(
			g_large %>% as_data_frame(what = "vertices") %>% as_tibble() %>% select(name) %>% rename(og_id = name),
			by = "og_id"
		) %>%
		inner_join(
			df_habitat_pref,
			by = "accession"
		) %>%
		left_join(
			df_env,
			by = "habitat"
		)
	df_genomes_prok_env$env[is.na(df_genomes_prok_env$env)] <- "Others"
	df_genomes_prok_env_summary <- df_genomes_prok_env %>%
		group_by(og_id,env) %>%
		summarize(hp = sum(habitat_pref)) %>%
		ungroup() %>%
		pivot_wider(names_from = env, values_from = hp, values_fill = 0)
	return(df_genomes_prok_env_summary)
}

show_network_large_component_env <- function (g_large, df_og_env, show_label) {
	df_og_env <- df_og_env %>%
		mutate(
			total = select(., -og_id) %>% rowSums(na.rm = TRUE)
		) %>%
		pivot_longer(
			cols = c(-og_id,-total),
			names_to = "env",
			values_to = "hp"
		) %>%
		mutate(hp_ratio = hp / total) %>%
		filter(env == "gut" | env == "soil" | env == "marine")
		#filter(env != "Others")
	#df_og_env$env <- factor(df_og_env$env, levels = c("gut","soil","rhizosphere","marine","freshwater"))
	df_og_env$env <- factor(df_og_env$env, levels = c("gut","soil","marine"))
	df_og_env$hp_ratio[df_og_env$hp_ratio >= 0.3] <- 0.3
	df_g_large <- ggnetwork_fr(g_large) %>%
		left_join(
			df_og_env %>% rename(name = "og_id"),
			by = "name"
		)
	df_label_show <- df_g_large %>%
		filter(v_e == "vertex", env == "gut" | env == "marine") %>%
		group_by(env) %>%
		nest() %>%
		mutate(
			top_genes = map(data, ~ .x %>% arrange(desc(hp_ratio)) %>% head(9) %>% select(name, .data[[show_label]], hp_ratio))
		) %>%
		select(-data) %>%
		unnest(top_genes) %>%
		ungroup() %>%
		group_by(name, .data[[show_label]]) %>%
		summarize(count = dplyr::n()) %>%
		mutate(
			label_show = map2(name, .data[[show_label]], ~ if (is.na(.y)) {.x} else {.y})
		) %>%
		unnest(c(label_show))
	df_g_large <- left_join(
		df_g_large,
		df_label_show %>% select(-.data[[show_label]]),
		by = "name"
	)
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend),
			size = 1,
			color = "grey80",
			data = df_g_large %>% filter(v_e == "edge")
		) +
		geom_point(
			aes(x = x, y = y, fill = hp_ratio),
			shape = 21,
			color = "white",
			data = df_g_large %>% filter(v_e == "vertex")
		) +
		geom_text_repel(
			aes(x = x, y = y, label = label_show, color = hp_ratio),
			min.segment.length = 0,
			segment.size = 0.1,
			segment.linetype = 1,
			segment.color = "black",
			nudge_x = 1.02,
			direction = "y",
			hjust = 1,
			size = 2,
			data = df_g_large %>% filter(v_e == "vertex")
		) +
		xlim(0,1.2) +
		scale_fill_viridis_c(na.value = "grey70") +
		scale_color_viridis_c(na.value = "grey70") +
		facet_wrap(. ~ env) +
		labs(fill = "Habitat\npreference\nratio", color = "Habitat\npreference\nratio") +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	#if (show_label == "cog_gene") {
	#	df_kable <- df_g_large %>%
	#		filter(v_e == "vertex", !is.na(label_show)) %>%
	#		select(name,cog_fun_cat,cog_gene.x,cog_name) %>%
	#		distinct()
	#	colnames(df_kable) <- c("COG","Category","Symbol","Name")
	#	df_kable %>% knitr::kable()
	#}
	return(list(p_network, df_g_large))
}

show_network_large_component_env_degree <- function (g_large, df_og_env, show_label) {
	df_og_env <- df_og_env %>%
		mutate(
			total = select(., -og_id) %>% rowSums(na.rm = TRUE)
		) %>%
		pivot_longer(
			cols = c(-og_id,-total),
			names_to = "env",
			values_to = "hp"
		) %>%
		mutate(hp_ratio = hp / total) %>%
		filter(env == "gut" | env == "soil" | env == "marine")
		#filter(env != "Others")
	#df_og_env$env <- factor(df_og_env$env, levels = c("gut","soil","rhizosphere","marine","freshwater"))
	df_og_env$env <- factor(df_og_env$env, levels = c("gut","soil","marine"))
	df_og_env$hp_ratio[df_og_env$hp_ratio >= 0.3] <- 0.3
	df_g_large <- ggnetwork_fr(g_large) %>%
		left_join(
			df_og_env %>% rename(name = "og_id"),
			by = "name"
		)
	df_label_show <- df_g_large %>%
		filter(v_e == "vertex", env == "gut" | env == "marine" | degree >= 16) %>%
		group_by(env) %>%
		nest() %>%
		mutate(
			top_genes = map(data, ~ .x %>% arrange(desc(hp_ratio)) %>% head(9) %>% select(name, .data[[show_label]], hp_ratio))
		) %>%
		select(-data) %>%
		unnest(top_genes) %>%
		ungroup() %>%
		group_by(name, .data[[show_label]]) %>%
		summarize(count = dplyr::n()) %>%
		mutate(
			label_show = map2(name, .data[[show_label]], ~ if (is.na(.y)) {.x} else {.y})
		) %>%
		unnest(c(label_show))
	df_g_large <- left_join(
		df_g_large,
		df_label_show %>% select(-.data[[show_label]]),
		by = "name"
	)
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend),
			size = 1,
			color = "grey80",
			data = df_g_large %>% filter(v_e == "edge")
		) +
		geom_point(
			aes(x = x, y = y, size = degree, fill = hp_ratio),
			shape = 21,
			color = "white",
			data = df_g_large %>% filter(v_e == "vertex")
		) +
		guides(size = guide_legend(override.aes = list(fill = "grey30"))) +
		scale_size(range = c(1, 3)) +
		geom_text_repel(
			aes(x = x, y = y, label = label_show, color = hp_ratio),
			min.segment.length = 0,
			segment.size = 0.1,
			segment.linetype = 1,
			segment.color = "black",
			nudge_x = 1.02,
			direction = "y",
			hjust = 1,
			size = 3,
			data = df_g_large %>% filter(v_e == "vertex")
		) +
		xlim(0,1.2) +
		scale_fill_viridis_c(na.value = "grey70") +
		scale_color_viridis_c(na.value = "grey70") +
		facet_wrap(. ~ env) +
		labs(fill = "Habitat\npreference\nratio", color = "Habitat\npreference\nratio", size = "Degree") +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	#if (show_label == "cog_gene") {
	#	df_kable <- df_g_large %>%
	#		filter(v_e == "vertex", !is.na(label_show)) %>%
	#		select(name,cog_fun_cat,cog_gene.x,cog_name) %>%
	#		distinct()
	#	colnames(df_kable) <- c("COG","Category","Symbol","Name")
	#	df_kable %>% knitr::kable()
	#}
	return(list(p_network, df_g_large))
}

show_network_large_component_env_degree_single <- function (g_large, df_og_env, show_label, env_label) {
	df_og_env <- df_og_env %>%
		mutate(
			total = select(., -og_id) %>% rowSums(na.rm = TRUE)
		) %>%
		pivot_longer(
			cols = c(-og_id,-total),
			names_to = "env",
			values_to = "hp"
		) %>%
		mutate(hp_ratio = hp / total) %>%
		filter(env == env_label)
	df_og_env$hp_ratio[df_og_env$hp_ratio >= 0.3] <- 0.3
	set.seed(1)
	df_g_large <- ggnetwork_fr(g_large) %>%
		left_join(
			df_og_env %>% rename(name = "og_id"),
			by = "name"
		) %>%
	#df_label_show <- df_g_large %>%
		#filter(v_e == "vertex", env == "gut" | env == "marine" | degree >= 16) %>%
		#group_by(env) %>%
		#nest() %>%
		#mutate(
		#	top_genes = map(data, ~ .x %>% arrange(desc(hp_ratio)) %>% head(9) %>% select(name, .data[[show_label]], hp_ratio))
		#) %>%
		#select(-data) %>%
		#unnest(top_genes) %>%
		#ungroup() %>%
		#group_by(name, .data[[show_label]]) %>%
		#summarize(count = dplyr::n()) %>%
		mutate(
			label_show = map2(name, .data[[show_label]], ~ if (is.na(.y)) {.x} else {.y})
		) %>%
		unnest(c(label_show))
	#df_g_large <- left_join(
	#	df_g_large,
	#	df_label_show %>% select(-.data[[show_label]]),
	#	by = "name"
	#)
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend),
			size = 1,
			color = "grey80",
			data = df_g_large %>% filter(v_e == "edge")
		) +
		geom_point(
			aes(x = x, y = y, size = degree, fill = hp_ratio),
			shape = 21,
			color = "white",
			data = df_g_large %>% filter(v_e == "vertex")
		) +
		guides(size = guide_legend(override.aes = list(fill = "grey30"))) +
		scale_size(range = c(1, 4)) +
		geom_text_repel(
			aes(x = x, y = y, label = label_show, color = hp_ratio),
			min.segment.length = 0,
			segment.size = 0.1,
			segment.linetype = 1,
			segment.color = "black",
			#nudge_x = 1.02,
			#direction = "y",
			#hjust = 1,
			size = 3,
			data = df_g_large %>% filter(v_e == "vertex")
		) +
		xlim(0,1.2) +
		scale_fill_viridis_c(na.value = "grey70") +
		scale_color_viridis_c(na.value = "grey70") +
		labs(fill = "Habitat\npreference\nratio", color = "Habitat\npreference\nratio", size = "Degree") +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

show_network_fun_color_label <- function (g_graph, fun_profile, show_label) {
	set.seed(1)
	df_g_network <- ggnetwork_fr(g_graph) %>%
		as_tibble() %>%
		mutate(
			fun_single = strsplit(.data[[fun_profile]], ";"),
			label_show = map2(name, .data[[show_label]], ~ if (is.na(.y)) {.x} else {.y})
		) %>%
		unnest(c(fun_single, label_show))
	df_fun_color <- df_g_network %>%
		group_by(fun_single) %>%
		summarize(count = dplyr::n()) %>%
		arrange(desc(count)) %>%
		head(11) %>%
		mutate(fun_color = fun_single)
	df_fun_color$fun_color <- factor(df_fun_color$fun_color, levels = c(df_fun_color$fun_color, "Others"))
	df_g_network <- df_g_network %>%
		left_join(df_fun_color, by = "fun_single")
	#if (is.na(df_g_network$fun_color) %>% sum() >= 1) {
		df_g_network$fun_color[is.na(df_g_network$fun_color)] <- "Others"
	#}
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend),
			size = 1,
			color = "grey80",
			data = df_g_network %>% filter(v_e == "edge")
		) +
		geom_point(
			aes(x = x, y = y, size = degree, fill = fun_color),
			shape = 21,
			color = "white",
			data = df_g_network %>% filter(v_e == "vertex")
		) +
		#guides(size = guide_legend(override.aes = list(fill = "grey30"))) +
		guides(size = "none") +
		scale_size(range = c(1, 4)) +
		geom_text_repel(
			aes(x = x, y = y, label = label_show, color = fun_color),
			max.overlaps = Inf,
			min.segment.length = 0,
			segment.size = 0.1,
			segment.linetype = 1,
			segment.color = "black",
			#nudge_x = 1.02,
			#direction = "y",
			#hjust = 1,
			size = 3,
			data = df_g_network %>% filter(v_e == "vertex")
		) +
		xlim(0,1.2) +
		scale_fill_manual(values = cbp12, na.value = "grey70") +
		scale_color_manual(values = cbp12, na.value = "grey70") +
		guides(fill = guide_legend(override.aes = list(size = 3))) +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

show_network_cluster_color <- function (g_graph, n_cluster, bray_thre, e_size, v_size) {
	set.seed(1)
	df_g_network <- ggnetwork_fr(g_graph) %>%
		as_tibble()
	df_g_network$cluster_colors <- factor(df_g_network$cluster_colors, levels = as.character(1:n_cluster))
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend, color = bray_col),
			size = e_size,
			#color = "grey80",
			#data = df_g_network %>% filter(v_e == "edge")
			data = df_g_network %>% filter(v_e == "edge") %>% mutate(bray_col = bray >= bray_thre)
		) +
		scale_color_manual(values = c("grey70", hcl.colors(2,"Zissou 1")[2]), na.value = "grey70") +
		geom_point(
			#aes(x = x, y = y, size = degree, fill = cluster_colors),
			aes(x = x, y = y, fill = cluster_colors),
			shape = 21,
			color = "white",
			size = v_size,
			data = df_g_network %>% filter(v_e == "vertex")
		) +
		scale_fill_manual(values = cols25(15), na.value = "black") +
		guides(fill = guide_legend(override.aes = list(size = 3))) +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

show_network_cluster_color2 <- function (g_graph, n_cluster, bray_thre, e_size, v_size) {
	set.seed(1)
	df_g_network <- ggnetwork_fr(g_graph) %>%
		as_tibble()
	df_g_network$cluster_colors <- factor(df_g_network$cluster_colors, levels = as.character(1:n_cluster))
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend, color = bray_col),
			size = e_size,
			#color = "grey80",
			#data = df_g_network %>% filter(v_e == "edge")
			data = df_g_network %>% filter(v_e == "edge") %>% mutate(bray_col = bray >= bray_thre)
		) +
		scale_color_manual(values = c("grey70", hcl.colors(2,"Zissou 1")[2]), na.value = "grey70") +
		geom_point(
			#aes(x = x, y = y, size = degree, fill = cluster_colors),
			aes(x = x, y = y, fill = cluster_colors),
			shape = 21,
			color = "white",
			size = v_size,
			data = df_g_network %>% filter(v_e == "vertex")
		) +
		scale_fill_manual(values = cols25(15), na.value = "grey80") +
		guides(fill = guide_legend(override.aes = list(size = 3))) +
		coord_equal() +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

show_network_cluster_color3 <- function (g_graph, n_cluster, bray_thre, e_size, v_size, df_ogs) {
  df_ogs$cluster_id <- factor(df_ogs$cluster_id, levels = c(levels(df_ogs$cluster_id), n_cluster + 1))
  df_ogs$cluster_id[is.na(df_ogs$cluster_id)] <- n_cluster + 1
  set.seed(1)
	df_g_network <- ggnetwork_fr(g_graph) %>%
		as_tibble() %>% left_join(df_ogs %>% dplyr::rename(name = "OGs"), by = "name")
	df_g_network$cluster_id <- factor(df_g_network$cluster_id, levels = as.character(1:(n_cluster + 1)))
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend, color = bray_col),
			size = e_size,
			data = df_g_network %>% filter(v_e == "edge") %>% mutate(bray_col = bray >= bray_thre)
		) +
		scale_color_manual(values = c("grey70", hcl.colors(2,"Zissou 1")[2]), na.value = "grey70") +
		geom_point(
			aes(x = x, y = y, fill = cluster_id, shape = habitat_group),
			#shape = 21,
			color = "white",
			size = v_size,
			data = df_g_network %>% filter(v_e == "vertex")
		) +
		scale_fill_manual(values = c(cols25(15)[1:n_cluster],"black"), na.value = "grey80", drop = FALSE) +
	  scale_shape_manual(values = c(23,21,24,22,25), na.value = 21, drop = FALSE) +
		guides(fill = guide_legend(override.aes = list(size = 3))) +
		coord_equal() +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

contra_network_large_component <- function (cog_fun, fun_cats, fun_colors, cogs, cog_labels) {
	fun_cat <- strsplit(fun_cats, ";") %>% unlist()
	df_og_replacement_cog_ar <- load_og_replacement("cog", "archaea")
	df_og_replacement_cog_bac <- load_og_replacement("cog", "bacteria")
	df_og_replacement_cog_prok <- og_replacement_prokaryotes(df_og_replacement_cog_ar, df_og_replacement_cog_bac)
	df_network <- read_network_edge_stats("cog", edge_thre_contra, "contrapositive") %>%
		mutate(edge_id = paste(og_id_A, og_id_B, sep = "_")) %>%
		left_join(
			df_og_replacement_cog_prok %>% mutate(edge_id = paste(og_id_A, og_id_B, sep = "_")) %>% select(-og_id_A, -og_id_B),
			by = "edge_id"
		) %>%
		filter(n_replace_prok >= replace_thre)
	df_def <- read_cog2020_def_fixed() %>%
		right_join(
			tibble(cog_id = unique(c(df_network$og_id_A, df_network$og_id_B))),
			by = "cog_id"
		) %>%
		rename(og_id = "cog_id")
	g_network <<- df_network %>% graph_from_data_frame(directed = FALSE, vertices = df_def)
	compo_g_network <- components(g_network)
	large_compo_vertices <- compo_g_network$membership[compo_g_network$membership == which(compo_g_network$csize %in% max(compo_g_network$csize))] %>% names()
	g_large_component <- subgraph(g_network, large_compo_vertices)
	V(g_large_component)$betweenness <- betweenness(g_large_component)
	V(g_large_component)[unlist(strsplit(cogs, ";"))]$plot_text <- unlist(strsplit(cog_labels, ";"))
	set.seed(1)
	df_g_large_component <- g_large_component %>%
		ggnetwork(layout = igraph::layout_with_fr(g_large_component)) %>%
		as_tibble() %>%
		mutate(
			x_same = x == xend,
			y_same = y == yend,
			v_e = map2(x_same, y_same, ~ if (.x == TRUE & .y == TRUE) {"vertex"} else {"edge"}),
			v_shape = map(.data[[cog_fun]], ~ if (length(grep("(R|S)", .x)) >= 1) {"23"} else {"21"})
		) %>%
		unnest(c(v_e, v_shape)) %>%
		left_join(tibble(cog_fun_cat = fun_cat, vertex_color = fun_cat), by = "cog_fun_cat")
	df_g_large_component$vertex_color[is.na(df_g_large_component$vertex_color)] <- "Others"
	df_g_large_component$vertex_color <- factor(df_g_large_component$vertex_color, levels = c(fun_cat, "Others"))
	df_g_large_component$n_replace_prok[df_g_large_component$n_replace_prok >= 100] <- 100
	p_network <- ggplot() +
		geom_segment(
			#aes(x = x, y = y, xend = xend, yend = yend, size = n_tax),
			aes(x = x, y = y, xend = xend, yend = yend),
			size = 1,
			color = "grey70",
			data = df_g_large_component %>% filter(v_e == "edge")
		) +
		geom_point(
			aes(x = x, y = y, fill = vertex_color, size = betweenness),
			color = "white",
			shape = 21,
			data = df_g_large_component %>% filter(v_e == "vertex", !is.na(vertex_color))
		) +
		scale_size_continuous(range = c(2,5)) +
		scale_fill_manual(values = fun_colors, na.value = "black") +
		labs(fill = "COG Category") +
		guides(size = "none") +
		geom_text_repel(
			aes(x = x, y = y, label = plot_text),
			data = df_g_large_component %>% filter(v_e == "vertex", !is.na(vertex_color))
		) +
		#coord_fixed(ratio = 1) +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(list(p_network, df_g_large_component, g_large_component))
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

habitat_analysis_prokatlas <- function (OG, g_network, show_label) {
	if (OG == "cog") {
		hp_thre <- 0
		df_contingency_habitat_hp0 <-
			fread(paste0("data/tsv/contingency_tsv_cog2020_habitat_hp", hp_thre, ".tsv"), header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble()
		hp_thre <- 5
		df_contingency_habitat_hp5 <-
			fread(paste0("data/tsv/contingency_tsv_cog2020_habitat_hp", hp_thre, ".tsv"), header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble()
		hp_thre <- 10
		df_contingency_habitat_hp10 <-
			fread(paste0("data/tsv/contingency_tsv_cog2020_habitat_hp", hp_thre, ".tsv"), header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble()
		in_file <- "data/tsv/network_edge_stats_cog.tsv"
		df_def <- read_cog2020_def_fixed() %>%
			rename(og_id = "cog_id") %>%
			mutate(node_type = "OG")
	} else if (OG == "kegg") {
		hp_thre <- 0
		df_contingency_habitat_hp0 <-
			fread(paste0("data/tsv/contingency_tsv_ko_habitat_hp", hp_thre, ".tsv"), header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>%
			rename(og_id_old = og_id) %>%
			mutate(og_id = sub("^ko:", "", og_id_old))
		hp_thre <- 5
		df_contingency_habitat_hp5 <-
			fread(paste0("data/tsv/contingency_tsv_ko_habitat_hp", hp_thre, ".tsv"), header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>%
			rename(og_id_old = og_id) %>%
			mutate(og_id = sub("^ko:", "", og_id_old))
		hp_thre <- 10
		df_contingency_habitat_hp10 <-
			fread(paste0("data/tsv/contingency_tsv_ko_habitat_hp", hp_thre, ".tsv"), header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>%
			rename(og_id_old = og_id) %>%
			mutate(og_id = sub("^ko:", "", og_id_old))
		in_file <- "data/tsv/network_edge_stats_kegg.tsv"
		df_def <- read_ko_def_na() %>%
			mutate(
				og_id = sub("^ko:","",ko_id),
				node_type = "OG"
			) %>%
			select(og_id, names(.))
	}
	node_ogs <- g_network %>% V() %>% names()
	df_network_og_habitat_hp0 <- inner_join(
		tibble(og_id = node_ogs),
		df_contingency_habitat_hp0 %>%
			filter(q_value <= 0.05, cooccurring_ratio_hp > 0.5) %>%
			mutate(FT_ratio = FT / (TF + FT), TT_ratio = TT/ (TT + FF)),
		by = "og_id"
	)
	df_network_og_habitat_hp5 <- inner_join(
		tibble(og_id = node_ogs),
		df_contingency_habitat_hp5 %>%
			filter(q_value <= 0.05, cooccurring_ratio_hp > 0.5) %>%
			mutate(FT_ratio = FT / (TF + FT), TT_ratio = TT/ (TT + FF)),
		by = "og_id"
	)
	df_network_og_habitat_hp10 <- inner_join(
		tibble(og_id = node_ogs),
		df_contingency_habitat_hp10 %>%
			filter(q_value <= 0.05, cooccurring_ratio_hp > 0.5) %>%
			mutate(FT_ratio = FT / (TF + FT), TT_ratio = TT/ (TT + FF)),
		by = "og_id"
	)
	df_network_og_habitat <- rbind(
		df_network_og_habitat_hp5 %>% mutate(hp_threshold = "HP_5"),
		df_network_og_habitat_hp0 %>% mutate(hp_threshold = "HP_0"),
		df_network_og_habitat_hp10 %>% mutate(hp_threshold = "HP_10")
	)
	df_network_og_habitat$hp_threshold <- factor(df_network_og_habitat$hp_threshold, levels = c("HP_0", "HP_5", "HP_10"))
	df_network_og_habitat_select <- df_network_og_habitat %>%
		filter(TT_ratio > 0.1, TT_ratio < 0.9, cooccurring_ratio_hp > 0.8)
	df_def_env <- bind_rows(
		df_def %>%
			inner_join(df_network_og_habitat_select %>% distinct(og_id), by = "og_id"),
		tibble(og_id = df_network_og_habitat_select$habitat %>% unique()) %>%
			mutate(node_type = "Habitat")
	)
	df_network_og_habitat_select <- left_join(
		df_network_og_habitat_select %>% mutate(pair_id = paste(og_id, habitat, sep = "--")),
		df_network_og_habitat_select %>%
			group_by(og_id,habitat) %>%
			summarize(hp_threshold_nr = paste(hp_threshold, collapse = ";")) %>%
			ungroup() %>%
	 		mutate(pair_id = paste(og_id, habitat, sep = "--")) %>%
			select(-og_id,-habitat),
		by = "pair_id"
	)
	g_network_og_habitat <- df_network_og_habitat_select %>%
		distinct(og_id,habitat,hp_threshold_nr) %>%
		graph_from_data_frame(directed = FALSE, vertices = df_def_env)
	set.seed(1)
	df_g_network <- g_network_og_habitat %>%
		ggnetwork(layout = igraph::layout_with_fr(g_network_og_habitat)) %>%
		as_tibble() %>%
		mutate(
			x_same = x == xend,
			y_same = y == yend,
			v_e = map2(x_same, y_same, ~ if (.x == TRUE & .y == TRUE) {"vertex"} else {"edge"})
		) %>%
		unnest(c(v_e)) %>%
		left_join(
			tibble(
				name = c("human_gut","mouse_gut","pig_gut","chicken_gut","human_bile","marine","freshwater","soil","anaerobic_digester","compost","fossil","glacier","biosolids","indoor","rock_porewater"),
				habitat_color = c("Guts","Guts","Guts","Guts","Host-associated","marine","freshwater","soil","anaerobic_digester","compost","Others","Others","Others","Others","Others")
			),
			by = "name"
		)
	df_g_network$habitat_color <- factor(df_g_network$habitat_color, levels = c("Guts","Host-associated","soil","freshwater","marine","anaerobic_digester","compost","Others"))
	df_label_show <- rbind(
		df_g_network %>%
			mutate(
				label_show = map2(name, .data[[show_label]], ~ if (is.na(.y)) {.x} else {.y})
			) %>%
			unnest(c(label_show)) %>%
			filter(node_type == "OG") %>%
			select(name,label_show),
		df_g_network %>%
			mutate(
				label_show = map2(name, habitat_color, ~ if (!is.na(.y) & .y == "Others") {.x} else if (.x == "human_gut") {.x} else {NA})
			) %>%
			unnest(c(label_show)) %>%
			filter(node_type == "Habitat", !is.na(label_show)) %>%
			select(name,label_show)
	)
	df_g_network <- df_g_network %>%
		left_join(
			df_label_show,
			by = "name"
		)
	df_g_network$hp_threshold_nr <- factor(df_g_network$hp_threshold_nr, levels = c("HP_0","HP_5","HP_10","HP_5;HP_10"))
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend, color = hp_threshold_nr),
			size = 1,
			data = df_g_network %>% filter(v_e == "edge")
		) +
		#scale_color_viridis_d() +
		scale_color_manual(values = rev(hcl.colors(5,"Light Grays"))) +
		guides(color = guide_legend(ncol = 2)) +
		labs(color = "Habitat preference") +
		new_scale_color() +
		geom_point(
			aes(x = x, y = y, shape = node_type, color = habitat_color),
			size = 5,
			data = df_g_network %>% filter(v_e == "vertex")
		) +
		geom_text_repel(
			#aes(x = x, y = y, label = label_show, color = habitat_color),
			aes(x = x, y = y, label = label_show),
			color = "grey10",
			size = 4,
			data = df_g_network %>% filter(v_e == "vertex") %>% distinct()
		) +
		scale_shape_manual(values = c(18,16)) +
		scale_color_manual(values = cbp8, na.value = "grey50", drop = FALSE) +
		labs(color = "Habitat", shape = "Node type") +
		guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			#legend.position = "left",
			legend.position = "right",
			plot.title = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

#visualize_og_habitat_tree <- function (OG, dom, As, Bs, habitats, tax_level) {
visualize_og_habitat_tree <- function (OG, dom, As, habitats, tax_level) {
	if (OG == "cog" & dom == "archaea") {
		og_path <- "data/pastml_cog2020_ar/"
		tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Archaea"
		node_plot_size <- 1.4
	} else if (OG == "cog" & dom == "bacteria") {
		og_path <- "data/pastml_cog2020_bac/"
		tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Bacteria"
		node_plot_size <- 0.8
	} else if (OG == "kegg" & dom == "archaea") {
		og_path <- "data/pastml_ko_ar/"
		tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Archaea"
		node_plot_size <- 1.4
	} else if (OG == "kegg" & dom == "bacteria") {
		og_path <- "data/pastml_ko_bac/"
		tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Bacteria"
		node_plot_size <- 0.8
	}
	df_habitat_pref <-
		read.csv("data/tsv/gtdb_ssu_prokatlas_habitat_pref.tsv", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
		as_tibble()
	OGs_A <- strsplit(As, ";") %>% unlist()
	#OGs_B <- strsplit(Bs, ";") %>% unlist()
	df_og_a <- tibble()
	#df_og_b <- tibble()
	for (i in OGs_A) {
		df_og_a <- rbind(
			df_og_a,
			read.csv(paste0(og_path, all_of(i), ".txt"), header = TRUE, sep = "\t", quote = "") %>%
				as_tibble() %>% rename(label = "node", TF = all_of(i)) %>% mutate(og_id = all_of(i))
		)
	}
	#for (i in OGs_B) {
	#	df_og_b <- rbind(
	#		df_og_b,
	#		read.csv(paste0(og_path, all_of(i), ".txt"), header = TRUE, sep = "\t", quote = "") %>%
	#			as_tibble() %>% rename(label = "node", TF = all_of(i)) %>% mutate(og_id = all_of(i))
	#	)
	#}
	gtdb_node_labeled_tree <- read.tree(tree_path)
	df_genomes <- read_genomes_analyzed() %>% select(accession, Domain:Species) %>% rename(label = "accession")
	df_og_a_pivot <- df_og_a %>% mutate(presence = 1) %>%
		pivot_wider(names_from = TF, values_from = presence, values_fill = 0) %>% rename(og_TRUE = "TRUE", og_FALSE = "FALSE") %>% mutate(OG_group = "A OGs")
	#df_og_b_pivot <- df_og_b %>% mutate(presence = 1) %>%
	#	pivot_wider(names_from = TF, values_from = presence, values_fill = 0) %>% rename(og_TRUE = "TRUE", og_FALSE = "FALSE") %>% mutate(OG_group = "B OGs")
	df_tree_mapped <- gtdb_node_labeled_tree %>% as_tibble() %>%
		#left_join(rbind(df_og_a_pivot, df_og_b_pivot), by = "label") %>%
		left_join(df_og_a_pivot, by = "label") %>%
		mutate(
			og_presence = map2(og_TRUE, og_FALSE, ~
				if (.x == 1 & .y == 0){
					1
				} else if (.x == 1 & .y == 1) {
					0.5
				} else if (.x == 0 & .y == 1) {
					0
				}
			)
		) %>%
		unnest(c(og_presence))
	p_tree <- ggtree(gtdb_node_labeled_tree, branch.length = 'none', layout = 'circular', color = "gray")
	df_tip_data <- df_tree_mapped[which(df_tree_mapped$label %in% gtdb_node_labeled_tree$tip.label), c("label", "og_id", "og_presence", "OG_group")] %>%
		pivot_longer(c(-label,-og_id,-OG_group), names_to = "ogs_tip", values_to = "pastml_tip") %>%
		mutate(
			og_presence = pastml_tip >= 1,
			OGs_tip = og_id
		) %>%
		left_join(p_tree$data, by = "label") %>%
		rename(ID = "label")
	#df_tip_data$OGs_tip <- factor(df_tip_data$OGs_tip, levels = c(OGs_A, OGs_B))
	df_tip_data$OGs_tip <- factor(df_tip_data$OGs_tip, levels = c(OGs_A))
	df_tip_data <- bind_rows(
		df_tip_data %>% mutate(
			habitat_pref = 100
		),
		df_habitat_pref[,c("accession", unlist(strsplit(habitats, ";")))] %>%
			pivot_longer(-accession, names_to = "habitat", values_to = "habitat_pref") %>%
			filter(habitat_pref > 0) %>%
			mutate(og_presence = TRUE) %>%
			rename(ID = "accession", OGs_tip = "habitat") %>%
			left_join(
				tibble(
					OGs_tip = c("human_gut","mouse_gut","pig_gut","chicken_gut","human_bile","marine","freshwater","soil","anaerobic_digester","compost","fossil","glacier","biosolids","indoor","rock_porewater"),
					OG_group = c("Guts","Guts","Guts","Guts","Host-associated","marine","freshwater","soil","anaerobic_digester","compost","Others","Others","Others","Others","Others")
				),
				by = "OGs_tip"
			)
		) #%>%
		#mutate(tile_alpha = cut(habitat_pref, breaks=c(0,5,10,100)))
	df_tip_data$habitat_pref[df_tip_data$habitat_pref >= 10] <- 10
	#df_tip_data$OGs_tip <- factor(df_tip_data$OGs_tip, levels = c(OGs_A, OGs_B, unlist(strsplit(habitats, ";"))))
	df_tip_data$OGs_tip <- factor(df_tip_data$OGs_tip, levels = c(OGs_A, unlist(strsplit(habitats, ";"))))
	df_tip_data$OG_group[is.na(df_tip_data$OG_group)] <- "Others"
	#df_tip_data$OG_group <- factor(df_tip_data$OG_group, levels = c("A OGs", "B OGs", c("Guts","Host-associated","soil","freshwater","marine","anaerobic_digester","compost","Others")))
	df_tip_data$OG_group <- factor(df_tip_data$OG_group, levels = c("A OGs", c("Guts","Host-associated","soil","freshwater","marine","anaerobic_digester","compost","Others")))
	df_tax <- df_genomes[which(df_genomes$label %in% p_tree$data[p_tree$data$isTip,]$label), c("label", tax_level)]
	df_tax_top <- df_tax %>% group_by(.data[[tax_level]]) %>% summarize(count = dplyr::n()) %>% arrange(desc(count)) %>%
		head(11) %>% select(-count) %>% mutate(tax_color = .data[[tax_level]])
	df_tax <- df_tax %>% left_join(df_tax_top, by = tax_level)
	df_tax$tax_color[is.na(df_tax$tax_color)] <- "Others"
	df_tax$tax_color <- factor(df_tax$tax_color, levels = c(df_tax_top$tax_color, "Others"))
	df_tax_data <- left_join(p_tree$data, df_tax, by = "label") %>% filter(isTip == TRUE)
	p_tree_plotted <- p_tree +
		geom_fruit(
			data = df_tip_data %>% filter(og_presence == TRUE) %>% select(ID,OGs_tip,OG_group,habitat_pref),
			geom = geom_tile,
			mapping = aes(y = ID, x = OGs_tip, fill = OG_group, alpha = habitat_pref),
			color = NA,
			#pwidth = 0.07 * (1 + length(c(OGs_A, OGs_B, unlist(strsplit(habitats, ";"))))),
			pwidth = 0.07 * (1 + length(c(OGs_A, unlist(strsplit(habitats, ";"))))),
			offset = 0.05
		) +
		labs(fill = "OGs and habitats") +
		#scale_fill_manual(values = c("red", "blue", cbp8[6], cbp8[7], cbp8[5], cbp8[2], cbp8[1], "grey50"), drop = FALSE) +
		scale_fill_manual(values = c("grey50", cbp8), drop = FALSE) +
		scale_alpha_continuous(range = c(0.1, 1), na.value = 0) +
		new_scale_fill() +
		geom_fruit(
			data = df_tax_data %>% rename(ID = "label") %>% select(ID, tax_color) %>% mutate(x_pos = 1),
			geom = geom_star,
			mapping = aes(y = ID, x = x_pos, fill = tax_color),
			starshape = 15,
			color = NA,
			alpha = 0.8,
			pwidth = 0.08
		) +
		scale_fill_manual(values = cbp12, na.value = cbp12[12]) +
		guides(alpha = "none") +
		layout_rectangular() +
		labs(fill = tax_level) +
		theme(
			text = element_text(family = "Arial", size = 12),
			legend.position = "left",
			plot.tag = element_text(face = "bold", size = 18),
			plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_tree_plotted)
}

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
	df_gmgc_habitats$bar_color[grep("gut", df_gmgc_habitats$habitat)] <- "Guts"
	df_gmgc_habitats$bar_color[grep("marine", df_gmgc_habitats$habitat)] <- "Marine"
	df_gmgc_habitats$bar_color[grep("soil", df_gmgc_habitats$habitat)] <- "Soil"
	df_gmgc_habitats$bar_color[grep("freshwater", df_gmgc_habitats$habitat)] <- "Freshwater"
	df_gmgc_habitats$bar_color <- factor(df_gmgc_habitats$bar_color, levels = c("Total Unigenes", "Guts", "Marine", "Freshwater", "Soil", "Others"))
	p_gmgc_habitats <- df_gmgc_habitats %>%
		ggplot(aes(x = n_unigene, y = habitat, fill = bar_color)) +
		geom_bar(stat = "identity") +
		scale_fill_manual(values = c("black", cbp8[6], cbp8[5], cbp8[2], cbp8[1], "grey50")) +
		labs(
			subtitle = sub("cog", "COG", ogid), 
			x = "Number of Unigenes",
			y = "Habitats"
		) +
		theme(
			legend.position = "none"
		)
	return(p_gmgc_habitats)
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


gmgc10_env_gene_abundance <- function (OG) {
	if (OG == "cog") {
		df_gmgc <- fread("data/tsv/gmgc10_cog.tsv.gz", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>%
			select(-unigene,-eggnog) %>%
			group_by(OGs) %>%
			summarize_all(list(sum)) %>%
			ungroup() %>%
			pivot_longer(cols = -OGs, names_to = "environment", values_to = "n_unigenes")
	} else if (OG == "kegg") {
		df_gmgc <- fread("data/tsv/gmgc10_kegg.tsv.gz", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>%
			select(-unigene,-eggnog) %>%
			group_by(OGs) %>%
			summarize_all(list(sum)) %>%
			ungroup() %>%
			pivot_longer(cols = -OGs, names_to = "environment", values_to = "n_unigenes")
	}
	df_gmgc_summary <- df_gmgc %>%
		mutate(
			envs = map(environment, ~
				if (length(grep("gut",.x)) > 0) {
					"gut"
				} else if (length(grep("human",.x)) > 0) {
					"human-associated"
				} else if (.x == "amplicon" | .x == "isolate") {
					NA
				} else {
					.x
				}
			)
		) %>%
		unnest(c(envs)) %>%
		filter(!is.na(envs)) %>%
		group_by(OGs,envs) %>%
		summarize(sum_unigenes = sum(n_unigenes)) %>%
		ungroup()
	df_gmgc_summary$envs <- factor(df_gmgc_summary$envs, levels = c("gut","human-associated","wastewater","marine","freshwater","soil","built-environment"))
	return(df_gmgc_summary)
}

gmgc10_env_norm_gene_abundance_network <- function (g_network, df_gmgc_summary) {
	set.seed(1)
	df_n_envs <- read.csv("raw/gmgc10/n_samples_envs.tsv", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
		as_tibble() %>% dplyr::rename(environment = envs) %>%
		mutate(
			envs = map(environment, ~
				if (length(grep("gut",.x)) > 0) {
					"gut"
				} else if (length(grep("human",.x)) > 0) {
					"human-associated"
				} else if (.x == "amplicon" | .x == "isolate") {
					NA
				} else {
					.x
				}
			)
		) %>%
		unnest(c(envs)) %>%
		filter(!is.na(envs)) %>%
		group_by(envs) %>%
		summarize(sum_samples = sum(n_sanples)) %>%
		ungroup()
	df_gmgc_summary_norm <- df_gmgc_summary %>%
		dplyr::rename(name = "OGs") %>%
		left_join(df_n_envs, by = "envs") %>%
		mutate(norm_unigenes = sum_unigenes / sum_samples) %>%
		select(-sum_unigenes,-sum_samples) %>%
		pivot_wider(names_from = envs, values_from = norm_unigenes)
	df_g_network <- ggnetwork_fr(g_network)
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend),
			size = 1,
			color = "grey70",
			data = df_g_network %>% filter(v_e == "edge")
		) +
		geom_scatterpie(
			aes(x = x, y = y, group = name, r = 0.01),
			cols = c("gut","human-associated","marine","freshwater","soil","built-environment"),
			color = NA,
			data = left_join(
					df_g_network,
					df_gmgc_summary_norm,
					by = "name"
				) %>% filter(v_e == "vertex")
		) +
		coord_equal() +
		#geom_text_repel(
		#	aes(x = x, y = y, label = label_show, color = label_show),
		#	min.segment.length = 0,
		#	segment.size = 0.1,
		#	segment.linetype = 1,
		#	segment.color = "black",
		#	nudge_x = 1.02,
		#	direction = "y",
		#	hjust = 1,
		#	size = 3,
		#	data = df_g_network %>% filter(v_e == "vertex", !is.na(label_show))
		#) +
		#xlim(0,1.2) +
		scale_fill_manual(values = cbp12) +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

gmgc10_env_gene_abundance_network <- function (g_network, df_gmgc_summary) {
	set.seed(1)
	df_g_network <- ggnetwork_fr(g_network)
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend),
			size = 1,
			color = "grey70",
			data = df_g_network %>% filter(v_e == "edge")
		) +
		geom_scatterpie(
			aes(x = x, y = y, group = name, r = 0.01),
			cols = c("gut","human-associated","wastewater","marine","freshwater","soil","built-environment"),
			color = NA,
			data = left_join(
					df_g_network,
					df_gmgc_summary %>% dplyr::rename(name = "OGs") %>% pivot_wider(names_from = envs, values_from = sum_unigenes),
					by = "name"
				) %>% filter(v_e == "vertex")
		) +
		coord_equal() +
		#geom_text_repel(
		#	aes(x = x, y = y, label = label_show, color = label_show),
		#	min.segment.length = 0,
		#	segment.size = 0.1,
		#	segment.linetype = 1,
		#	segment.color = "black",
		#	nudge_x = 1.02,
		#	direction = "y",
		#	hjust = 1,
		#	size = 3,
		#	data = df_g_network %>% filter(v_e == "vertex", !is.na(label_show))
		#) +
		#xlim(0,1.2) +
		scale_fill_manual(values = unname(palette.colors(palette = "Tableau 10"))) +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

gmgc10_env_gene_abundance_network_taxon <- function (g_network, df_edge_taxon_summary, taxon) {
	#g_network <- g_cog_8_contra_bray
	#df_edge_taxon_summary <- df_edge_phylum_8taxa
	#taxon <- "p__Firmicutes_A"
	df_edge_tax <- df_edge_taxon_summary %>% filter(tax == taxon)
	igraph::E(g_network)[df_edge_tax$from %--% df_edge_tax$to]
	
	igraph::E(g_network)[df_edge_tax$edge_id[3]]
	igraph::E(g_network)["COG0757|COG0710"]
	igraph::E(g_network)["COG0757" %--% c("COG0710", "COG0710")]
	igraph::E(g_network)["COG0710|COG0757"]










	df_network_tax <- as_data_frame(g_network, what = "edges") %>%
		as_tibble() #%>%
		select(from, to) #%>%
	#paired_ogs_bray <- function (og_a, og_b) {
	#	bray_pair <- matrix(
	#		c(df_edge_tax_norm_t[,og_a] %>% pull(), df_edge_tax_norm_t[,og_b] %>% pull()),
	#		nrow = 2,
	#		byrow = TRUE
	#	) %>%
	#	vegdist(method = "bray") %>%
	#	as.vector()
	#	return(bray_pair)
	#}
	#df_network_bray <- as_data_frame(g_network, what = "edges") %>%
	#	as_tibble() %>%
	#	select(from, to) #%>%
	#	#mutate(bray = map2(from, to, ~ paired_ogs_bray(.x, .y))) %>%
	#	#unnest(c(bray))
	#df_network_bray[df_network_bray$from %in% ogs_disappear_in_gmgc,"bray"] <- NA
	#df_network_bray[df_network_bray$to %in% ogs_disappear_in_gmgc,"bray"] <- NA
	#igraph::E(g_network)$bray <- df_network_bray$bray
	return(g_network)
}

gmgc10_env_gene_abundance_network_bray <- function (g_network, df_gmgc_summary) {
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

gmgc10_env_gene_relative_abundance <- function (df_gmgc_summary) {
	df_gmgc_summary_norm <- df_gmgc_summary %>%
    mutate(total_unigenes = rowSums(df_gmgc_summary %>% select(-OGs))) %>% 
    pivot_longer(cols = c(-OGs, -total_unigenes), names_to = "habitat", values_to = "n_unigenes") %>% 
    mutate(ratio_unigenes = n_unigenes / total_unigenes) %>% 
    select(-n_unigenes) %>% 
    pivot_wider(names_from = "habitat", values_from = "ratio_unigenes")
  return(df_gmgc_summary_norm)
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

gmgc10_env_gene_abundance_network_bray_vis <- function (g_network, df_gmgc_summary, ogs) {
	set.seed(1)
	df_g_network <- ggnetwork_fr(g_network)
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend, color = bray),
			size = 1,
			data = df_g_network %>% filter(v_e == "edge")
		) +
		scale_color_viridis_c(na.value = "grey70") +
		#scale_color_viridis_c(option = "inferno", begin = 0, end = 0.8, na.value = "grey70") +
		geom_scatterpie(
			aes(x = x, y = y, group = name, r = 0.009),
			cols = c("gut","human-associated","wastewater","marine","freshwater","soil","built-environment"),
			color = NA,
			data = left_join(
					df_g_network,
					df_gmgc_summary %>% dplyr::rename(name = "OGs") %>% pivot_wider(names_from = envs, values_from = sum_unigenes),
					by = "name"
				) %>% filter(v_e == "vertex")
		) +
		coord_equal() +
		geom_text_repel(
			aes(x = x, y = y, label = name_custom_2),
			min.segment.length = 0,
			segment.size = 0.1,
			segment.linetype = 1,
			segment.color = "black",
			size = 3,
			data = df_g_network %>% filter(v_e == "vertex") %>% mutate(name_custom_2 = sub("\n.*","",name_custom)) %>% filter(name %in% ogs)
		) +
		scale_fill_manual(values = unname(palette.colors(palette = "Tableau 10"))) +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0, size = 18, face = "bold"),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

gmgc10_env_gene_abundance_network_bray_vis_enlarge <- function (g_network, g_network_ori, df_gmgc_summary, ogs, r_val, seg_size) {
	set.seed(1)
	bray_min <- igraph::E(g_network_ori)$bray %>% min(na.rm = TRUE)
	bray_max <- igraph::E(g_network_ori)$bray %>% max(na.rm = TRUE)
	df_g_network <- ggnetwork_fr(g_network)
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend, color = bray),
			size = seg_size,
			data = df_g_network %>% filter(v_e == "edge")
		) +
		scale_color_viridis_c(limits = c(bray_min, bray_max), na.value = "grey70") +
		geom_scatterpie(
			aes(x = x, y = y, group = name, r = r_val),
			cols = c("gut","human-associated","wastewater","marine","freshwater","soil","built-environment"),
			color = NA,
			data = left_join(
					df_g_network,
					df_gmgc_summary %>% dplyr::rename(name = "OGs") %>% pivot_wider(names_from = envs, values_from = sum_unigenes),
					by = "name"
				) %>% filter(v_e == "vertex")
		) +
		guides(fill = "none") +
		coord_equal() +
		geom_text_repel(
			aes(x = x, y = y, label = name_custom_2),
			min.segment.length = 0,
			segment.size = 0.1,
			segment.linetype = 1,
			segment.color = "black",
			size = 3,
			data = df_g_network %>% filter(v_e == "vertex") %>% mutate(name_custom_2 = sub("\n.*","",name_custom)) %>% filter(name %in% ogs)
		) +
		scale_fill_manual(values = unname(palette.colors(palette = "Tableau 10"))) +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0, size = 18, face = "bold"),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

gmgc10_env_gene_abundance_network_bray_vis_enlarge2 <- function (g_network, g_network_ori, df_gmgc_summary, ogs, r_val, seg_size, bray_thre) {
	set.seed(1)
	bray_min <- igraph::E(g_network_ori)$bray %>% min(na.rm = TRUE)
	bray_max <- igraph::E(g_network_ori)$bray %>% max(na.rm = TRUE)
	df_g_network <- ggnetwork_fr(g_network)
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend, color = bray, linetype = bray_tf),
			size = seg_size,
			data = df_g_network %>% filter(v_e == "edge") %>% mutate(bray_tf = bray >= bray_thre)
		) +
		scale_color_viridis_c(limits = c(bray_min, bray_max), na.value = "grey70") +
		scale_linetype_manual(values = c("21","solid"), na.value = "21") +
		#scale_linetype_manual(values = c("longdash","solid"), na.value = "longdash") +
		geom_scatterpie(
			aes(x = x, y = y, group = name, r = r_val),
			cols = c("gut","human-associated","wastewater","marine","freshwater","soil","built-environment"),
			color = "black",
			size = 0.1,
			data = left_join(
					df_g_network,
					df_gmgc_summary %>% dplyr::rename(name = "OGs") %>% pivot_wider(names_from = envs, values_from = sum_unigenes),
					by = "name"
				) %>% filter(v_e == "vertex")
		) +
		guides(fill = "none") +
		coord_equal() +
		geom_text_repel(
			aes(x = x, y = y, label = name_custom_2),
			min.segment.length = 0,
			segment.size = 0.1,
			segment.linetype = 1,
			segment.color = "black",
			size = 3,
			data = df_g_network %>% filter(v_e == "vertex") %>% mutate(name_custom_2 = sub("\n.*","",name_custom)) %>% filter(name %in% ogs)
		) +
		scale_fill_manual(values = unname(palette.colors(palette = "Tableau 10"))) +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0, size = 18, face = "bold"),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

gmgc10_env_gene_abundance_network_bray_vis_enlarge2_no_habitat <- function (g_network, g_network_ori, df_gmgc_summary, ogs, r_val, seg_size, bray_thre) {
	set.seed(1)
	bray_min <- igraph::E(g_network_ori)$bray %>% min(na.rm = TRUE)
	bray_max <- igraph::E(g_network_ori)$bray %>% max(na.rm = TRUE)
	df_g_network <- ggnetwork_fr(g_network)
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend, color = bray, linetype = bray_tf),
			size = seg_size,
			data = df_g_network %>% filter(v_e == "edge") %>% mutate(bray_tf = bray >= bray_thre)
		) +
		scale_color_viridis_c(limits = c(bray_min, bray_max), na.value = "grey70") +
		scale_linetype_manual(values = c("21","solid"), na.value = "21") +
		#scale_linetype_manual(values = c("longdash","solid"), na.value = "longdash") +
		#geom_scatterpie(
		#	aes(x = x, y = y, group = name, r = r_val),
		#	cols = c("gut","human-associated","wastewater","marine","freshwater","soil","built-environment"),
		#	color = "black",
		#	size = 0.1,
		#	data = left_join(
		#			df_g_network,
		#			df_gmgc_summary %>% dplyr::rename(name = "OGs") %>% pivot_wider(names_from = envs, values_from = sum_unigenes),
		#			by = "name"
		#		) %>% filter(v_e == "vertex")
		#) +
		geom_point(
			aes(x = x, y = y),
			#cols = c("gut","human-associated","wastewater","marine","freshwater","soil","built-environment"),
			shape = 21,
			color = "grey50",
			fill = "grey80",
			size = r_val,
			data = df_g_network %>% filter(v_e == "vertex")
		) +
		guides(fill = "none") +
		coord_equal() +
		geom_text_repel(
			aes(x = x, y = y, label = name_custom_2),
			min.segment.length = 0,
			segment.size = 0.1,
			segment.linetype = 1,
			segment.color = "black",
			size = 3,
			data = df_g_network %>% filter(v_e == "vertex") %>% mutate(name_custom_2 = sub("\n.*","",name_custom)) %>% filter(name %in% ogs)
		) +
		scale_fill_manual(values = unname(palette.colors(palette = "Tableau 10"))) +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0, size = 18, face = "bold"),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

gmgc_bar <- function (df_gmgc_summary, ogs) {
	df_gmgc_bar <- df_gmgc_summary %>% filter(OGs %in% ogs)
	df_gmgc_bar$OGs <- factor(df_gmgc_bar$OGs, levels = rev(ogs))
	p_gmgc_bar <- ggplot(df_gmgc_bar) +
		geom_bar(
			aes(x = sum_unigenes, y = OGs, fill = envs),
			stat = "identity",
			position = "fill"
		) +
		scale_fill_manual(values = unname(palette.colors(palette = "Tableau 10")))
	return(p_gmgc_bar)
}

gmgc_bar_horizon <- function (df_gmgc_summary, df_ogs) {
	df_gmgc_bar <- df_gmgc_summary %>% right_join(df_ogs, by = "OGs")
	#df_gmgc_bar <- df_gmgc_summary %>% filter(OGs %in% ogs)
	#df_gmgc_bar$OGs <- factor(df_gmgc_bar$OGs, levels = rev(ogs))
	p_gmgc_bar <- ggplot(df_gmgc_bar) +
		geom_bar(
			aes(y = sum_unigenes, x = new_name, fill = envs),
			stat = "identity",
			position = "fill"
		) +
		scale_fill_manual(values = unname(palette.colors(palette = "Tableau 10"))) +
		labs(y = "Distribution", fill = "Habitat") +
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 13))
	return(p_gmgc_bar)
}

gmgc_bar_horizon2 <- function (df_gmgc_summary, df_ogs) {
  df_gmgc_bar <- df_gmgc_summary %>%
    pivot_longer(cols = c(-OGs, -total_unigenes), names_to = "envs", values_to = "unigene_ratio") %>% 
    right_join(df_ogs, by = "OGs")
  df_gmgc_bar$envs <- factor(df_gmgc_bar$envs, levels = c(
    c("marine"),
    c("human.oral","human.skin","human.nose","human.vagina"),
    c("soil"),
    c("human.gut","pig.gut","mouse.gut","cat.gut","dog.gut"),
    c("wastewater","freshwater","built.environment")
    ))
  p_gmgc_bar <- ggplot(df_gmgc_bar) +
    geom_bar(
      aes(y = unigene_ratio, x = new_name, fill = envs),
      stat = "identity",
      position = "fill"
    ) +
    scale_fill_manual(values = c(
      pals::stepped2()[1],
      pals::stepped2()[5:8],
      pals::stepped2()[9],
      pals::stepped2()[c(13:16,16)],
      pals::stepped2()[17:19]
      ), drop = FALSE) +
    labs(y = "Distribution", fill = "Habitat") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 13))
  return(p_gmgc_bar)
}

gmgc_bar_horizon3 <- function (df_gmgc_summary, ogs) {
  df_gmgc_bar <- df_gmgc_summary %>% filter(OGs %in% ogs) %>%
    pivot_longer(cols = c(-OGs, -total_unigenes), names_to = "envs", values_to = "unigene_ratio")
  df_gmgc_bar$envs <- factor(df_gmgc_bar$envs, levels = c(
    c("marine"),
    c("human.oral","human.skin","human.nose","human.vagina"),
    c("soil"),
    c("human.gut","pig.gut","mouse.gut","cat.gut","dog.gut"),
    c("wastewater","freshwater","built.environment")
    ))
  df_gmgc_bar$OGs <- factor(df_gmgc_bar$OGs, levels = ogs)
  p_gmgc_bar <- ggplot(df_gmgc_bar) +
    geom_bar(
      aes(y = unigene_ratio, x = OGs, fill = envs),
      stat = "identity",
      position = "fill"
    ) +
    scale_fill_manual(values = c(
      pals::stepped2()[1],
      pals::stepped2()[5:8],
      pals::stepped2()[9],
      pals::stepped2()[c(13:16,16)],
      pals::stepped2()[17:19]
      ), drop = FALSE) +
    labs(y = "Distribution", fill = "Habitat") + theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 13))
  return(p_gmgc_bar)
}

gmgc_binary_horizon <- function (df_ogs) {
  p_gmgc_binary <- ggplot(df_ogs) +
  geom_point(aes(x = new_name, y = "Habitat_group", fill = cluster_id, shape = habitat_group), color = "black", size = 4) +
	scale_fill_manual(values = cols25(15), na.value = "black", drop = FALSE) +
  scale_shape_manual(values = c(23,21,24,22), drop = FALSE) +
  theme_blank()
  return(p_gmgc_binary)
}

gmgc_binary_horizon2 <- function (df_ogs) {
  p_gmgc_binary <- ggplot(df_ogs) +
  geom_point(aes(x = new_name, y = "Habitat_group", fill = cluster_id, shape = habitat_group), color = "black", size = 3) +
	scale_fill_manual(values = cols25(15), na.value = "black", drop = FALSE) +
  scale_shape_manual(values = c(23,21,24,22,25), drop = FALSE) +
  theme_blank()
  return(p_gmgc_binary)
}

gmgc_bar_facet <- function (df_gmgc_summary, df_ogs_sort) {
	df_gmgc_bar <- inner_join(df_ogs_sort, df_gmgc_summary, by = "OGs")
	p_gmgc_bar <- ggplot(df_gmgc_bar) +
		geom_bar(
			aes(y = sum_unigenes, x = OGs, fill = envs),
			stat = "identity",
			position = "fill"
		) +
		scale_fill_manual(values = unname(palette.colors(palette = "Tableau 10"))) +
		facet_grid(. ~ cluster_id)
	return(p_gmgc_bar)
}

gmgc_bar_abundance <- function (df_gmgc_summary, ogs) {
	df_gmgc_bar <- df_gmgc_summary %>% filter(OGs %in% ogs)
	df_gmgc_bar$OGs <- factor(df_gmgc_bar$OGs, levels = rev(ogs))
	p_gmgc_bar <- ggplot(df_gmgc_bar) +
		geom_bar(
			aes(x = sum_unigenes, y = OGs, fill = envs),
			stat = "identity"
		) +
		scale_fill_manual(values = unname(palette.colors(palette = "Tableau 10")))
	return(p_gmgc_bar)
}

gmgc_pca <- function (df_gmgc_summary, ogs) {
	df_gmgc_ratio <- df_gmgc_summary %>% filter(OGs %in% ogs) %>%
		group_by(OGs) %>%
		nest() %>%
		ungroup() %>%
		mutate(total_unigenes = map(data, ~ .x$sum_unigenes %>% sum())) %>%
		unnest(c(data,total_unigenes)) %>%
		mutate(ratio_unigenes = sum_unigenes / total_unigenes) %>%
		select(-sum_unigenes,-total_unigenes) %>%
		pivot_wider(names_from = envs, values_from = ratio_unigenes)
	mx_gmgc_ratio <- as.matrix(df_gmgc_ratio %>% select(-OGs))
	rownames(mx_gmgc_ratio) <- df_gmgc_ratio$OGs
	pca <- prcomp(mx_gmgc_ratio, scale = F)
	return(pca)
}

gmgc_pca2 <- function (df_gmgc_summary, ogs) {
  df_gmgc_ratio <- df_gmgc_summary %>% filter(total_unigenes != 0)
  mx_gmgc_ratio <- as.matrix(df_gmgc_ratio %>% select(-OGs,-total_unigenes))
  rownames(mx_gmgc_ratio) <- df_gmgc_ratio$OGs
  pca <- prcomp(mx_gmgc_ratio, scale = F)
  return(pca)
}

gmgc10_env_gene_abundance_network_bray_vis_binary <- function (g_network, df_gmgc_summary, ogs, bray_thre) {
	#g_network <- g_cog_8_contra_bray
	#df_gmgc_summary <- df_gmgc_cog
	#ogs <- ""
	set.seed(1)
	df_g_network <- ggnetwork_fr(g_network)
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend, color = bray >= bray_thre),
			size = 1,
			data = df_g_network %>% filter(v_e == "edge")
		) +
		scale_color_manual(values = hcl.colors(2,"Zissou 1"), na.value = "grey70") +
		#scale_color_viridis_c(option = "inferno", begin = 0, end = 0.8, na.value = "grey70") +
		geom_scatterpie(
			aes(x = x, y = y, group = name, r = 0.009),
			cols = c("gut","human-associated","wastewater","marine","freshwater","soil","built-environment"),
			color = NA,
			data = left_join(
					df_g_network,
					df_gmgc_summary %>% dplyr::rename(name = "OGs") %>% pivot_wider(names_from = envs, values_from = sum_unigenes),
					by = "name"
				) %>% filter(v_e == "vertex")
		) +
		coord_equal() +
		geom_text_repel(
			aes(x = x, y = y, label = name_custom_2),
			min.segment.length = 0,
			segment.size = 0.1,
			segment.linetype = 1,
			segment.color = "black",
			size = 3,
			data = df_g_network %>% filter(v_e == "vertex") %>% mutate(name_custom_2 = sub("\n.*","",name_custom)) %>% filter(name %in% ogs)
		) +
		scale_fill_manual(values = unname(palette.colors(palette = "Tableau 10"))) +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0, size = 18, face = "bold"),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

gmgc10_env_gene_abundance_network_bray_vis_binary_no_habitat <- function (g_network, df_gmgc_summary, ogs, bray_thre) {
	#g_network <- g_cog_8_contra_bray
	#df_gmgc_summary <- df_gmgc_cog
	#ogs <- ""
	set.seed(1)
	df_g_network <- ggnetwork_fr(g_network)
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend, color = bray >= bray_thre),
			size = 1,
			data = df_g_network %>% filter(v_e == "edge")
		) +
		scale_color_manual(values = hcl.colors(2,"Zissou 1"), na.value = "grey70") +
		#scale_color_viridis_c(option = "inferno", begin = 0, end = 0.8, na.value = "grey70") +
		geom_point(
			aes(x = x, y = y),
			#cols = c("gut","human-associated","wastewater","marine","freshwater","soil","built-environment"),
			shape = 21,
			color = "grey50",
			fill = "grey80",
			size = 2,
			data = df_g_network %>% filter(v_e == "vertex")
		) +
		coord_equal() +
		geom_text_repel(
			aes(x = x, y = y, label = name_custom_2),
			min.segment.length = 0,
			segment.size = 0.1,
			segment.linetype = 1,
			segment.color = "black",
			size = 3,
			data = df_g_network %>% filter(v_e == "vertex") %>% mutate(name_custom_2 = sub("\n.*","",name_custom)) %>% filter(name %in% ogs)
		) +
		scale_fill_manual(values = unname(palette.colors(palette = "Tableau 10"))) +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0, size = 18, face = "bold"),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

gmgc10_env_gene_abundance_network_bray_vis_enlarge_binary <- function (g_network, df_gmgc_summary, ogs, r_val, seg_size) {
	set.seed(1)
	df_g_network <- ggnetwork_fr(g_network)
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend, color = bray >= 0.65),
			size = seg_size,
			data = df_g_network %>% filter(v_e == "edge")
		) +
		scale_color_manual(values = hcl.colors(2,"Zissou 1"), na.value = "grey70") +
		geom_scatterpie(
			aes(x = x, y = y, group = name, r = r_val),
			cols = c("gut","human-associated","wastewater","marine","freshwater","soil","built-environment"),
			color = NA,
			data = left_join(
					df_g_network,
					df_gmgc_summary %>% dplyr::rename(name = "OGs") %>% pivot_wider(names_from = envs, values_from = sum_unigenes),
					by = "name"
				) %>% filter(v_e == "vertex")
		) +
		guides(fill = "none") +
		coord_equal() +
		geom_text_repel(
			aes(x = x, y = y, label = name_custom_2),
			min.segment.length = 0,
			segment.size = 0.1,
			segment.linetype = 1,
			segment.color = "black",
			size = 3,
			data = df_g_network %>% filter(v_e == "vertex") %>% mutate(name_custom_2 = sub("\n.*","",name_custom)) %>% filter(name %in% ogs)
		) +
		scale_fill_manual(values = unname(palette.colors(palette = "Tableau 10"))) +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0, size = 18, face = "bold"),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

network_nodes_pick <- function (g_network, ogs) {
	v_ego <- ego(g_network, order = 1, nodes = ogs) %>%
		unlist() %>%
		names() %>%
		unique()
	sub_network <- subgraph(g_network, v_ego)
	return(sub_network)
}

gmgc_env_bray_plot <- function (g_network) {
	p_plot <- as_data_frame(g_network, what = "edges") %>%
		as_tibble() %>%
		ggplot(aes(y = "Edges", x = bray)) +
		geom_violin(fill = NA, draw_quantiles = c(0.25, 0.5, 0.75)) +
		geom_jitter(
			size = 1.5,
			alpha = 0.3,
			shape = 16
		) +
		labs(x = "Bray-Curtis") +
		theme(
			axis.title.y = element_blank(),
			axis.text.y = element_blank(),
			axis.ticks.y = element_blank()
		)
	return(p_plot)
}

gmgc_env_bray_plot_binary <- function (g_network, fun_profile) {
	p_plot <- as_data_frame(g_network, what = "edges") %>%
		as_tibble() %>%
		ggplot(aes(y = .data[[fun_profile]], x = bray)) +
		geom_jitter(
			size = 1.5,
			color = hcl.colors(2,"Zissou 1")[1],
			shape = 16
		) +
		geom_violin(fill = NA, draw_quantiles = c(0.25, 0.5, 0.75)) +
		labs(x = "Bray-Curtis")
	return(p_plot)
}

gmgc_env_bray_plot_binary2 <- function (g_network, fun_profile) {
	p_plot <- as_data_frame(g_network, what = "edges") %>%
		as_tibble() %>%
		mutate(bray_thre = bray >= 0.7) %>%
		ggplot(aes(y = .data[[fun_profile]], x = bray)) +
		geom_jitter(
			aes(color = bray_thre),
			size = 1.5,
			#color = hcl.colors(2,"Zissou 1")[1],
			shape = 16
		) +
		scale_color_manual(values = hcl.colors(2,"Zissou 1")) +
		geom_violin(fill = NA, draw_quantiles = c(0.25, 0.5, 0.75)) +
		labs(x = "Bray-Curtis")
	return(p_plot)
}

gmgc_env_bray_plot_binary3 <- function (g_network, fun_profile) {
	p_plot <- as_data_frame(g_network, what = "edges") %>%
		as_tibble() %>%
		mutate(bray_thre = bray >= 0.7) %>%
		ggplot(aes(y = .data[[fun_profile]], x = bray)) +
		geom_jitter(
			aes(color = bray_thre),
			size = 1.5,
			#color = hcl.colors(2,"Zissou 1")[1],
			shape = 16
		) +
		scale_color_manual(values = c("grey70", hcl.colors(2,"Zissou 1")[2])) +
		geom_violin(fill = NA, draw_quantiles = c(0.25, 0.5, 0.75)) +
		labs(x = "Bray-Curtis")
	return(p_plot)
}

gmgc_env_bray_plot_binary4 <- function (g_network, fun_profile, n_tax_thre, bray_thre) {
  set.seed(123)
	p_plot <- as_data_frame(g_network, what = "edges") %>%
		as_tibble() %>%
	  filter(n_tax >= n_tax_thre) %>% 
		mutate(bray_thre = bray >= bray_thre) %>%
		ggplot(aes(x = .data[[fun_profile]], y = bray)) +
		geom_jitter(
			aes(color = bray_thre),
			size = 1.5,
			#color = hcl.colors(2,"Zissou 1")[1],
			shape = 16
		) +
		scale_color_manual(values = c("grey70", hcl.colors(2,"Zissou 1")[2])) +
		geom_boxplot(fill = NA, outlier.color = NA) +
		labs(x = "Bray-Curtis")
	return(p_plot)
}

gmgc_env_bray_plot_binary5_df <- function (g_network, n_tax_thre, bray_thre) {
  set.seed(123)
  df_plot <- as_data_frame(g_network, what = "edges") %>%
    as_tibble() %>%
    filter(n_tax >= n_tax_thre) %>% 
    mutate(bray_thre = bray >= bray_thre, n_tax_threshold = paste0("n_tax >= ", as.character(n_tax_thre)))
  return(df_plot)
}

gmgc_env_bray_plot_binary5 <- function (df_network, fun_profile, bray_thre) {
  set.seed(123)
  p_plot <- df_network %>%
    ggplot(aes(x = .data[[fun_profile]], y = bray, color = .data[[fun_profile]], fill = .data[[fun_profile]])) +
    geom_jitter(color = "grey90", size = 1.5, shape = 1) +
    geom_hline(yintercept = bray_thre, linetype = "dashed") +
    geom_violin(alpha = 0.5) +
    geom_boxplot(color = "black", fill = "white", outlier.color = NA, width = 0.1) +
    scale_color_manual(values = hcl.colors(2,"TealRose")) +
    scale_fill_manual(values = hcl.colors(2,"TealRose")) +
    facet_grid(. ~ n_tax_threshold) +
    theme(legend.position = "none") +
    labs(y = "Bray-Curtis", x = "Sharing function")
  return(p_plot)
}

gmgc10_env_gene_abundance_network_bray_cytoscape <- function (g_graph, df_gmgc_summary, network_title) {
	df_gmgc_summary_norm <- df_gmgc_summary %>%
		group_by(OGs) %>%
		nest() %>%
		mutate(total_unigenes = map(data, ~ sum(.x$sum_unigenes))) %>%
		unnest(c(data,total_unigenes)) %>%
		mutate(norm_unigenes = sum_unigenes / total_unigenes) %>%
		select(-sum_unigenes,-total_unigenes) %>%
		pivot_wider(names_from = envs, values_from = norm_unigenes) %>%
		ungroup() %>%
		dplyr::rename(name = "OGs")
	df_node_def <- as_data_frame(g_graph, what = "vertices") %>%
		as_tibble() %>%
		left_join(df_gmgc_summary_norm, by = "name")
	g_network <- graph_from_data_frame(as_data_frame(g_graph, what = "edges"), directed = FALSE, vertices = df_node_def)
	createNetworkFromIgraph(g_network, network_title)
	setLayoutProperties('kamada-kawai', list(singlePartition=1, randomize=0))
	layoutNetwork('kamada-kawai')
	setNodeShapeDefault('ellipse')
	setNodeCustomPieChart(
		c("gut","human.associated","wastewater","marine","freshwater","soil","built.environment"),
		colors = unname(palette.colors(palette = "Tableau 10"))
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

network_genomic_feature <- function (g_network, OG, gfeat, bin_width) {
	if (OG == "cog") {
		df_gfeat <- fread("data/tsv/genomic_feature_cog.tsv.gz", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>% select(OG,.data[[gfeat]]) #%>% group_by(OG) %>% nest() %>% ungroup()
	} else if (OG == "kegg") {
		df_gfeat <- fread("data/tsv/genomic_feature_kegg.tsv.gz", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>% dplyr::rename(OG_old = "OG") %>% mutate(OG = sub("ko:", "", OG_old)) %>% select(OG,.data[[gfeat]]) #%>% group_by(OG) %>% nest() %>% ungroup()
	}
	df_gfeat_og <- df_gfeat %>%
		group_by(OG) %>%
		summarize(mean_val = mean(.data[[gfeat]]))
	mean_vals <- df_gfeat_og$mean_val
	names(mean_vals) <- df_gfeat_og$OG
	df_gfeat_count <- df_gfeat %>%
		mutate(bin = cut(.data[[gfeat]], breaks = bin_width)) %>%
		group_by(OG,bin) %>%
		summarize(count = dplyr::n()) %>%
		ungroup() %>%
		pivot_wider(names_from = bin, values_from = count)
	df_gfeat_count[is.na(df_gfeat_count)] <- 0
	mx_gfeat_prob <- df_gfeat_count %>% select(-OG) %>% as.matrix() %>% prop.table(margin = 1) #%>% t() %>% as.data.frame()
	rownames(mx_gfeat_prob) <- df_gfeat_count$OG
	paired_ogs_kl <- function (og_a, og_b) {
		kl <- KL(mx_gfeat_prob[c(og_a,og_b),])
		return(kl)
	}
	df_network_kl <- as_data_frame(g_network, what = "edges") %>%
		as_tibble() %>%
		select(from, to) %>%
		mutate(kl = map2(from, to, ~ paired_ogs_kl(.x, .y))) %>%
		unnest(c(kl))
	igraph::E(g_network)$kl <- df_network_kl$kl
	V(g_network)$mean_val <- mean_vals[V(g_network)$name]
	return(g_network)
}

network_genomic_feature_kl <- function (g_network, df_gfeat, OG, gfeat, bin_width) {
  df_gfeat_og <- df_gfeat %>%
    group_by(OG) %>%
    summarize(mean_val = mean(.data[[gfeat]]))
  mean_vals <- df_gfeat_og$mean_val
  names(mean_vals) <- df_gfeat_og$OG
  df_gfeat_count <- df_gfeat %>%
    mutate(bin = cut(.data[[gfeat]], breaks = bin_width)) %>%
    group_by(OG,bin) %>%
    summarize(count = dplyr::n()) %>%
    ungroup() %>%
    pivot_wider(names_from = bin, values_from = count)
  df_gfeat_count[is.na(df_gfeat_count)] <- 0
  mx_gfeat_prob <- df_gfeat_count %>% select(-OG) %>% as.matrix() %>% prop.table(margin = 1) #%>% t() %>% as.data.frame()
  rownames(mx_gfeat_prob) <- df_gfeat_count$OG
  paired_ogs_kl <- function (og_a, og_b) {
    kl <- KL(mx_gfeat_prob[c(og_a,og_b),])
    return(kl)
  }
  df_network_kl <- as_data_frame(g_network, what = "edges") %>%
    as_tibble() %>%
    select(from, to) %>%
    mutate(kl = map2(from, to, ~ paired_ogs_kl(.x, .y))) %>%
    unnest(c(kl))
  igraph::E(g_network)$kl <- df_network_kl$kl
  V(g_network)$mean_val <- mean_vals[V(g_network)$name]
  return(g_network)
}

network_genomic_feature_vis_kl <- function (g_network, kl_thre) {
	set.seed(1)
	ogs <- subgraph.edges(g_network, igraph::E(g_network)[igraph::E(g_network)$kl >= kl_thre]) %>% V() %>% .$name
	df_g_network <- ggnetwork_fr(g_network)
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend, color = kl, size = n_tax),
			data = df_g_network %>% filter(v_e == "edge")
		) +
		scale_size_continuous(range = c(0.8,2.5)) +
		guides(size = "none") +
		new_scale("size") +
		scale_color_viridis_c(na.value = "grey70") +
		geom_point(
			aes(x = x, y = y, fill = mean_val, size = degree),
			shape = 21,
			color = "white",
			data = df_g_network %>% filter(v_e == "vertex")
		) +
		guides(size = "none") +
		scale_size_continuous(range = c(1.3,2.6)) +
		scale_fill_viridis_c(option = "inferno", begin = 0, end = 0.8, na.value = "grey70") +
		geom_text_repel(
			aes(x = x, y = y, label = name_custom_2),
			min.segment.length = 0,
			segment.size = 0.1,
			segment.linetype = 1,
			segment.color = "black",
			size = 3,
			data = df_g_network %>% filter(v_e == "vertex") %>% mutate(name_custom_2 = sub("\n.*","",name_custom)) %>% filter(name %in% ogs)
		) +
		coord_equal() +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0, size = 18, face = "bold"),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

histo_genomic_feature <- function (g_network, OG, gfeat, ogs, scale10) {
	if (OG == "cog") {
		df_gfeat <- fread("data/tsv/genomic_feature_cog.tsv.gz", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>% select(OG,.data[[gfeat]]) #%>% group_by(OG) %>% nest() %>% ungroup()
	} else if (OG == "kegg") {
		df_gfeat <- fread("data/tsv/genomic_feature_kegg.tsv.gz", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>% select(OG,.data[[gfeat]]) #%>% group_by(OG) %>% nest() %>% ungroup()
	}
	og_names <- sub("\n.*", "", V(g_network)$name_custom)[V(g_network)$name %in% ogs]
	names(og_names) <- V(g_network)$name[V(g_network)$name %in% ogs]
	df_gfeat_plot <- df_gfeat %>%
		filter(OG %in% ogs) %>%
		left_join(
			tibble(
				OG = names(og_names),
				Name = factor(og_names, levels = og_names[ogs])
			),
			by = "OG"
		)
	p_gfeat <- df_gfeat_plot %>%
		ggplot(aes(x = .data[[gfeat]] / scale10, fill = Name)) +
		geom_histogram(alpha = 0.5, position = "identity") +
		scale_fill_manual(values = hcl.colors(2,"Blue-Red 2"))
	return(p_gfeat)
}

dens_genomic_feature <- function (g_network, OG, gfeat, ogs, scale10) {
	if (OG == "cog") {
		df_gfeat <- fread("data/tsv/genomic_feature_cog.tsv.gz", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>% select(OG,.data[[gfeat]]) #%>% group_by(OG) %>% nest() %>% ungroup()
	} else if (OG == "kegg") {
		df_gfeat <- fread("data/tsv/genomic_feature_kegg.tsv.gz", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
			as_tibble() %>% select(OG,.data[[gfeat]]) %>% dplyr::rename(OG_old = "OG") %>% mutate(OG = sub("^ko:", "", OG_old))
	}
	og_names <- sub("\n.*", "", V(g_network)$name_custom)[V(g_network)$name %in% ogs]
	names(og_names) <- V(g_network)$name[V(g_network)$name %in% ogs]
	df_gfeat_plot <- df_gfeat %>%
		filter(OG %in% ogs) %>%
		left_join(
			tibble(
				OG = names(og_names),
				Name = factor(og_names, levels = og_names[ogs])
			),
			by = "OG"
		)
	p_gfeat <- df_gfeat_plot %>%
		ggplot(aes(x = .data[[gfeat]] / scale10, fill = Name, color = Name)) +
		geom_density(alpha = 0.5) +
		scale_fill_manual(values = hcl.colors(2,"Blue-Red 2")) +
		scale_color_manual(values = hcl.colors(2,"Blue-Red 2"))
	return(p_gfeat)
}

mkdf_genomic_feature <- function (g_network, OG, gfeat, ogs) {
  if (OG == "cog") {
    df_gfeat <- fread("data/tsv/genomic_feature_cog.tsv.gz", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
      as_tibble() %>% select(accession, OG,.data[[gfeat]]) #%>% group_by(OG) %>% nest() %>% ungroup()
  } else if (OG == "kegg") {
    df_gfeat <- fread("data/tsv/genomic_feature_kegg.tsv.gz", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
      as_tibble() %>% select(accession, OG,.data[[gfeat]]) %>% dplyr::rename(OG_old = "OG") %>% mutate(OG = sub("^ko:", "", OG_old))
  }
  og_names <- sub("\n.*", "", V(g_network)$name_custom)[V(g_network)$name %in% ogs]
  names(og_names) <- V(g_network)$name[V(g_network)$name %in% ogs]
  df_gfeat_plot <- df_gfeat %>%
    filter(OG %in% ogs) %>%
    left_join(
      tibble(
        OG = names(og_names),
        Name = factor(og_names, levels = og_names[ogs])
      ),
      by = "OG"
    )
  return(df_gfeat_plot)
}

visualize_og_displacement_rect_gfeat <- function (OG, dom, og_a, og_b, gene_a, gene_b, tax_level) {
	if (OG == "cog" & dom == "archaea") {
		og_a_path <- paste0("data/pastml_cog2020_ar/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_cog2020_ar/", og_b, ".txt")
		tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Archaea"
		node_plot_size <- 1.4
	} else if (OG == "cog" & dom == "bacteria") {
		og_a_path <- paste0("data/pastml_cog2020_bac/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_cog2020_bac/", og_b, ".txt")
		tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Bacteria"
		node_plot_size <- 0.8
	} else if (OG == "kegg" & dom == "archaea") {
		og_a_path <- paste0("data/pastml_ko_ar/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_ko_ar/", og_b, ".txt")
		tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Archaea"
		node_plot_size <- 1.4
	} else if (OG == "kegg" & dom == "bacteria") {
		og_a_path <- paste0("data/pastml_ko_bac/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_ko_bac/", og_b, ".txt")
		tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Bacteria"
		node_plot_size <- 0.8
	}
	df_og_a <- read.csv(og_a_path, header = TRUE, sep = "\t", quote = "") %>% as_tibble()
	df_og_b <- read.csv(og_b_path, header = TRUE, sep = "\t", quote = "") %>% as_tibble()
	gtdb_node_labeled_tree <- read.tree(tree_path)
	df_genomes <- read_genomes_analyzed() %>% select(accession, Domain:Species, genome_size, gc_percentage) %>% rename(label = "accession")
	df_og_a_pivot <- df_og_a %>% mutate(presence = 1) %>% group_by(node) %>%
		pivot_wider(names_from = all_of(og_a), values_from = presence, values_fill = 0) %>% rename(label = "node", og_a_TRUE = "TRUE", og_a_FALSE = "FALSE")
	df_og_b_pivot <- df_og_b %>% mutate(presence = 1) %>% group_by(node) %>%
		pivot_wider(names_from = all_of(og_b), values_from = presence, values_fill = 0) %>% rename(label = "node", og_b_TRUE = "TRUE", og_b_FALSE = "FALSE")
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
		unnest(c(og_a_presence, og_b_presence))
	p_tree <- ggtree(gtdb_node_labeled_tree, branch.length = 'none', layout = 'circular', color = "gray", size = 0.5)
	df_tip_data <- df_tree_mapped[which(df_tree_mapped$label %in% gtdb_node_labeled_tree$tip.label), c("label", "og_a_presence", "og_b_presence")] %>%
		pivot_longer(-label, names_to = "ogs_tip", values_to = "pastml_tip") %>%
		mutate(
			og_presence = pastml_tip >= 1,
			OGs_tip = map(ogs_tip, ~ if (.x == "og_a_presence") {paste0(og_a, "\n", gene_a)} else if (.x == "og_b_presence") {paste0(og_b, "\n", gene_b)})
		) %>%
		unnest(c(OGs_tip)) %>%
		left_join(p_tree$data, by = "label") %>%
		rename(ID = "label")
	df_tip_data$OGs_tip <- factor(df_tip_data$OGs_tip, levels = c(paste0(og_a, "\n", gene_a), paste0(og_b, "\n", gene_b)))
	df_tax <- df_genomes[which(df_genomes$label %in% p_tree$data[p_tree$data$isTip,]$label), c("label", tax_level,"genome_size","gc_percentage")]
	df_tax_top <- df_tax %>% group_by(.data[[tax_level]]) %>% summarize(count = dplyr::n()) %>% arrange(desc(count)) %>%
		head(11) %>% select(-count) %>% mutate(tax_color = .data[[tax_level]])
	df_tax <- df_tax %>% left_join(df_tax_top, by = tax_level)
	df_tax$tax_color[is.na(df_tax$tax_color)] <- "Others"
	df_tax$tax_color <- factor(df_tax$tax_color, levels = c(df_tax_top$tax_color, "Others"))
	df_tax_data <- left_join(p_tree$data, df_tax, by = "label") %>% filter(isTip == TRUE)
	df_node_data <- df_tree_mapped[which(df_tree_mapped$label %in% gtdb_node_labeled_tree$node.label), c("label", "og_a_presence", "og_b_presence")] %>%
		pivot_longer(-label, names_to = "ogs_node", values_to = "pastml_node") %>%
		mutate(
			node_size = factor(pastml_node),
			OGs_node = map(ogs_node, ~ if (.x == "og_a_presence") {paste0(og_a, "\n", gene_a)} else if (.x == "og_b_presence") {paste0(og_b, "\n", gene_b)})
		) %>%
		unnest(c(OGs_node)) %>%
		left_join(p_tree$data, by = "label")
	df_node_data$node_size[df_node_data$pastml_node == 0] <- NA
	df_node_data$OGs_node <- factor(df_node_data$OGs_node, levels = c(paste0(og_a, "\n", gene_a), paste0(og_b, "\n", gene_b)))
	p_tree_plotted <- p_tree +
		geom_point(
			data = df_node_data,
			aes(x = x, y = y, color = OGs_node, size = node_size),
			shape = 16,
			alpha = 0.5
		) +
		guides(size = "none") +
		scale_size_manual(values = c(0, node_plot_size/2, node_plot_size), na.translate = FALSE, drop = FALSE) +
		scale_color_manual(values = hcl.colors(2,"Blue-Red 2"), guide = "none") +
		new_scale_fill() +
		geom_fruit(
			data = df_tip_data %>% filter(og_presence == TRUE) %>% select(ID, OGs_tip),
			geom = geom_tile,
			mapping = aes(y = ID, x = OGs_tip, fill = OGs_tip),
			color = NA,
			pwidth = 0.4,
			offset = 0.12
		) +
		labs(fill = "OGs") +
		scale_fill_manual(values = hcl.colors(2,"Blue-Red 2")) +
		new_scale_fill() +
		geom_fruit(
			data = df_tax_data %>% rename(ID = "label") %>% select(ID, tax_color, genome_size),
			geom = geom_bar,
			mapping = aes(y = ID, x = genome_size, fill = tax_color),
			orientation = "y", 
                    	stat = "identity",
			color = NA,
			pwidth = 0.3,
			offset = 0.27
		) +
		geom_fruit(
			data = df_tax_data %>% rename(ID = "label") %>% select(ID, tax_color, gc_percentage),
			geom = geom_bar,
			mapping = aes(y = ID, x = gc_percentage, fill = tax_color),
			orientation = "y", 
                    	stat = "identity",
			color = NA,
			pwidth = 0.3,
			offset = 0.2
		) +
		guides(fill = "none") +
		scale_fill_manual(values = cbp12, na.value = cbp12[12]) +
		labs(size = "Presence\nprobability", fill = tax_level) +
    		theme_void() +
  		layout_rectangular() +
		coord_flip() +
		scale_y_reverse() +
    		theme(
			plot.background = element_rect(fill = "white",color = NA),
    			text = element_text(family = "Arial", size = 12),
    			axis.title.x = element_text(size = 13),
    			axis.title.y = element_text(size = 13),
    			plot.title = element_text(face = "bold", hjust = 0, size = 18),
    			plot.subtitle = element_text(hjust = 0.5, size = 12),
    			plot.caption = element_text(hjust = 0.5, size = 10),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13),
			legend.position = "right",
			legend.title = element_blank()
		)
}

visualize_og_displacement_rect_gfeat2 <- function (OG, dom, og_a, og_b, gene_a, gene_b, tax_level) {
  if (OG == "cog" & dom == "archaea") {
    og_a_path <- paste0("data/pastml_cog2020_ar/", og_a, ".txt")
    og_b_path <- paste0("data/pastml_cog2020_ar/", og_b, ".txt")
    tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Archaea"
    node_plot_size <- 1.4
  } else if (OG == "cog" & dom == "bacteria") {
    og_a_path <- paste0("data/pastml_cog2020_bac/", og_a, ".txt")
    og_b_path <- paste0("data/pastml_cog2020_bac/", og_b, ".txt")
    tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Bacteria"
    node_plot_size <- 0.8
  } else if (OG == "kegg" & dom == "archaea") {
    og_a_path <- paste0("data/pastml_ko_ar/", og_a, ".txt")
    og_b_path <- paste0("data/pastml_ko_ar/", og_b, ".txt")
    tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Archaea"
    node_plot_size <- 1.4
  } else if (OG == "kegg" & dom == "bacteria") {
    og_a_path <- paste0("data/pastml_ko_bac/", og_a, ".txt")
    og_b_path <- paste0("data/pastml_ko_bac/", og_b, ".txt")
    tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Bacteria"
    node_plot_size <- 0.8
  }
  df_og_a <- read.csv(og_a_path, header = TRUE, sep = "\t", quote = "") %>% as_tibble()
  df_og_b <- read.csv(og_b_path, header = TRUE, sep = "\t", quote = "") %>% as_tibble()
  gtdb_node_labeled_tree <- read.tree(tree_path)
  df_genomes <- read_genomes_analyzed() %>% select(accession, Domain:Species, genome_size, gc_percentage) %>% rename(label = "accession")
  df_og_a_pivot <- df_og_a %>% mutate(presence = 1) %>% group_by(node) %>%
    pivot_wider(names_from = all_of(og_a), values_from = presence, values_fill = 0) %>% rename(label = "node", og_a_TRUE = "TRUE", og_a_FALSE = "FALSE")
  df_og_b_pivot <- df_og_b %>% mutate(presence = 1) %>% group_by(node) %>%
    pivot_wider(names_from = all_of(og_b), values_from = presence, values_fill = 0) %>% rename(label = "node", og_b_TRUE = "TRUE", og_b_FALSE = "FALSE")
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
    unnest(c(og_a_presence, og_b_presence))
  p_tree <- ggtree(gtdb_node_labeled_tree, branch.length = 'none', layout = 'circular', color = "gray", size = 0.5)
  df_tip_data <- df_tree_mapped[which(df_tree_mapped$label %in% gtdb_node_labeled_tree$tip.label), c("label", "og_a_presence", "og_b_presence")] %>%
    pivot_longer(-label, names_to = "ogs_tip", values_to = "pastml_tip") %>%
    mutate(
      og_presence = pastml_tip >= 1,
      OGs_tip = map(ogs_tip, ~ if (.x == "og_a_presence") {paste0(og_a, "\n", gene_a)} else if (.x == "og_b_presence") {paste0(og_b, "\n", gene_b)})
    ) %>%
    unnest(c(OGs_tip)) %>%
    left_join(p_tree$data, by = "label") %>%
    rename(ID = "label")
  df_tip_data$OGs_tip <- factor(df_tip_data$OGs_tip, levels = c(paste0(og_a, "\n", gene_a), paste0(og_b, "\n", gene_b)))
  df_tax <- df_genomes[which(df_genomes$label %in% p_tree$data[p_tree$data$isTip,]$label), c("label", tax_level,"genome_size","gc_percentage")]
  df_tax_top <- df_tax %>% group_by(.data[[tax_level]]) %>% summarize(count = dplyr::n()) %>% arrange(desc(count)) %>%
    head(11) %>% select(-count) %>% mutate(tax_color = .data[[tax_level]])
  df_tax <- df_tax %>% left_join(df_tax_top, by = tax_level)
  df_tax$tax_color[is.na(df_tax$tax_color)] <- "Others"
  df_tax$tax_color <- factor(df_tax$tax_color, levels = c(df_tax_top$tax_color, "Others"))
  df_tax_data <- left_join(p_tree$data, df_tax, by = "label") %>% filter(isTip == TRUE)
  df_node_data <- df_tree_mapped[which(df_tree_mapped$label %in% gtdb_node_labeled_tree$node.label), c("label", "og_a_presence", "og_b_presence")] %>%
    pivot_longer(-label, names_to = "ogs_node", values_to = "pastml_node") %>%
    mutate(
      node_size = factor(pastml_node),
      OGs_node = map(ogs_node, ~ if (.x == "og_a_presence") {paste0(og_a, "\n", gene_a)} else if (.x == "og_b_presence") {paste0(og_b, "\n", gene_b)})
    ) %>%
    unnest(c(OGs_node)) %>%
    left_join(p_tree$data, by = "label")
  df_node_data$node_size[df_node_data$pastml_node == 0] <- NA
  df_node_data$OGs_node <- factor(df_node_data$OGs_node, levels = c(paste0(og_a, "\n", gene_a), paste0(og_b, "\n", gene_b)))
  p_tree_plotted <- p_tree +
    geom_point(
      data = df_node_data,
      aes(x = x, y = y, color = OGs_node, size = node_size),
      shape = 16,
      alpha = 0.5
    ) +
    guides(size = "none") +
    scale_size_manual(values = c(0, node_plot_size/2, node_plot_size), na.translate = FALSE, drop = FALSE) +
    scale_color_manual(values = hcl.colors(2,"Blue-Red 2"), guide = "none") +
    new_scale_fill() +
    geom_fruit(
      data = df_tip_data %>% filter(og_presence == TRUE) %>% select(ID, OGs_tip),
      geom = geom_tile,
      mapping = aes(y = ID, x = OGs_tip, fill = OGs_tip),
      color = NA,
      pwidth = 0.4,
      offset = 0.13
    ) +
    labs(fill = "OGs") +
    scale_fill_manual(values = hcl.colors(2,"Blue-Red 2")) +
    new_scale_fill() +
    geom_fruit(
      data = df_tax_data %>% rename(ID = "label") %>% select(ID, tax_color, genome_size),
      geom = geom_bar,
      mapping = aes(y = ID, x = genome_size, fill = tax_color),
      orientation = "y", 
      stat = "identity",
      color = NA,
      pwidth = 0.55,
      offset = 0.25,
      axis.params=list(
        axis       = "x",
        text.size  = 1,
        hjust      = -0.3,
        vjust      = 0.5,
        nbreak     = 2,
        line.size = 0.05
      ),
      grid.params=list()
    ) +
    geom_fruit(
      data = df_tax_data %>% rename(ID = "label") %>% select(ID, tax_color, gc_percentage),
      geom = geom_bar,
      mapping = aes(y = ID, x = gc_percentage, fill = tax_color),
      orientation = "y", 
      stat = "identity",
      color = NA,
      pwidth = 0.55,
      offset = 0.06,
      axis.params=list(
        axis       = "x",
        text.size  = 1,
        hjust      = -0.3,
        vjust      = 0.5,
        nbreak     = 2,
        line.size = 0.05
      ),
      grid.params=list()
    ) +
    guides(fill = "none") +
    scale_fill_manual(values = cbp12, na.value = cbp12[12]) +
    labs(size = "Presence\nprobability", fill = tax_level) +
    theme_void() +
    layout_rectangular() +
    ylim(-2000,26841) +
    coord_flip() +
    scale_y_reverse() +
    theme(
      plot.background = element_rect(fill = "white",color = NA),
      text = element_text(family = "Arial", size = 12),
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text(size = 13),
      plot.title = element_text(face = "bold", hjust = 0, size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      plot.caption = element_text(hjust = 0.5, size = 10),
      plot.tag = element_text(face = "bold", size = 18),
      strip.text.x = element_text(size = 13),
      strip.text.y = element_text(size = 13),
      legend.position = "right",
      legend.title = element_blank()
    )
}

phylogenetic_anova <- function (OG, dom, og_a) {
  if (OG == "cog" & dom == "archaea") {
    og_a_path <- paste0("data/pastml_csv_cog2020_ar/", og_a, ".csv")
    tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Archaea"
  } else if (OG == "cog" & dom == "bacteria") {
    og_a_path <- paste0("data/pastml_csv_cog2020_bac/", og_a, ".csv")
    tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Bacteria"
  } else if (OG == "kegg" & dom == "archaea") {
    og_a_path <- paste0("data/pastml_csv_ko_ar/", og_a, ".csv")
    tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Archaea"
  } else if (OG == "kegg" & dom == "bacteria") {
    og_a_path <- paste0("data/pastml_csv_ko_bac/", og_a, ".csv")
    tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
    Dom <- "Bacteria"
  }
  df_og_a <- read.csv(og_a_path, header = TRUE, sep = ",", quote = "") %>% as_tibble()
  gtdb_node_labeled_tree <- read.tree(tree_path)
  df_genomes <- read_genomes_analyzed() %>% select(accession, Domain:Species, genome_size, gc_percentage) %>%
    filter(Domain == paste0("d__", Dom)) %>% rename(label = "accession")
  df_genomes_ogs <- df_genomes %>%
    left_join(df_og_a %>% dplyr::rename(label = "accession"), by = "label")
  og_presence <- df_genomes_ogs[,og_a] %>% pull()
  names(og_presence) <- df_genomes_ogs$label
  gc_percent <- df_genomes_ogs$gc_percentage
  names(gc_percent) <- df_genomes_ogs$label
  phylANOVA(gtdb_node_labeled_tree, og_presence, gc_percent)
}

visualize_og_displacement_rect_gfeat_pick <- function (OG, dom, og_a, og_b, gene_a, gene_b, tax_level, pick_node) {
	if (OG == "cog" & dom == "archaea") {
		og_a_path <- paste0("data/pastml_cog2020_ar/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_cog2020_ar/", og_b, ".txt")
		tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Archaea"
		node_plot_size <- 1.4
	} else if (OG == "cog" & dom == "bacteria") {
		og_a_path <- paste0("data/pastml_cog2020_bac/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_cog2020_bac/", og_b, ".txt")
		tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Bacteria"
		node_plot_size <- 0.8
	} else if (OG == "kegg" & dom == "archaea") {
		og_a_path <- paste0("data/pastml_ko_ar/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_ko_ar/", og_b, ".txt")
		tree_path <- "raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Archaea"
		node_plot_size <- 1.4
	} else if (OG == "kegg" & dom == "bacteria") {
		og_a_path <- paste0("data/pastml_ko_bac/", og_a, ".txt")
		og_b_path <- paste0("data/pastml_ko_bac/", og_b, ".txt")
		tree_path <- "raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree"
		Dom <- "Bacteria"
		node_plot_size <- 0.8
	}
	df_og_a <- read.csv(og_a_path, header = TRUE, sep = "\t", quote = "") %>% as_tibble()
	df_og_b <- read.csv(og_b_path, header = TRUE, sep = "\t", quote = "") %>% as_tibble()
	gtdb_node_labeled_tree <- read.tree(tree_path) %>% tree_subset(node = pick_node, levels_back = 0)
	df_genomes <- read_genomes_analyzed() %>% select(accession, Domain:Species, genome_size, gc_percentage) %>% rename(label = "accession")
	df_og_a_pivot <- df_og_a %>% mutate(presence = 1) %>% group_by(node) %>%
		pivot_wider(names_from = all_of(og_a), values_from = presence, values_fill = 0) %>% rename(label = "node", og_a_TRUE = "TRUE", og_a_FALSE = "FALSE")
	df_og_b_pivot <- df_og_b %>% mutate(presence = 1) %>% group_by(node) %>%
		pivot_wider(names_from = all_of(og_b), values_from = presence, values_fill = 0) %>% rename(label = "node", og_b_TRUE = "TRUE", og_b_FALSE = "FALSE")
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
		unnest(c(og_a_presence, og_b_presence))
	p_tree <- ggtree(gtdb_node_labeled_tree, branch.length = 'none', layout = 'circular', color = "gray", size = 0.5)
	df_tip_data <- df_tree_mapped[which(df_tree_mapped$label %in% gtdb_node_labeled_tree$tip.label), c("label", "og_a_presence", "og_b_presence")] %>%
		pivot_longer(-label, names_to = "ogs_tip", values_to = "pastml_tip") %>%
		mutate(
			og_presence = pastml_tip >= 1,
			OGs_tip = map(ogs_tip, ~ if (.x == "og_a_presence") {paste0(og_a, "\n", gene_a)} else if (.x == "og_b_presence") {paste0(og_b, "\n", gene_b)})
		) %>%
		unnest(c(OGs_tip)) %>%
		left_join(p_tree$data, by = "label") %>%
		rename(ID = "label")
	df_tip_data$OGs_tip <- factor(df_tip_data$OGs_tip, levels = c(paste0(og_a, "\n", gene_a), paste0(og_b, "\n", gene_b)))
	df_tax <- df_genomes[which(df_genomes$label %in% p_tree$data[p_tree$data$isTip,]$label), c("label", tax_level,"genome_size","gc_percentage")]
	df_tax_top <- df_tax %>% group_by(.data[[tax_level]]) %>% summarize(count = dplyr::n()) %>% arrange(desc(count)) %>%
		head(11) %>% select(-count) %>% mutate(tax_color = .data[[tax_level]])
	df_tax <- df_tax %>% left_join(df_tax_top, by = tax_level)
	df_tax$tax_color[is.na(df_tax$tax_color)] <- "Others"
	df_tax$tax_color <- factor(df_tax$tax_color, levels = c(df_tax_top$tax_color, "Others"))
	df_tax_data <- left_join(p_tree$data, df_tax, by = "label") %>% filter(isTip == TRUE)
	df_node_data <- df_tree_mapped[which(df_tree_mapped$label %in% gtdb_node_labeled_tree$node.label), c("label", "og_a_presence", "og_b_presence")] %>%
		pivot_longer(-label, names_to = "ogs_node", values_to = "pastml_node") %>%
		mutate(
			node_size = factor(pastml_node),
			OGs_node = map(ogs_node, ~ if (.x == "og_a_presence") {paste0(og_a, "\n", gene_a)} else if (.x == "og_b_presence") {paste0(og_b, "\n", gene_b)})
		) %>%
		unnest(c(OGs_node)) %>%
		left_join(p_tree$data, by = "label")
	df_node_data$node_size[df_node_data$pastml_node == 0] <- NA
	df_node_data$OGs_node <- factor(df_node_data$OGs_node, levels = c(paste0(og_a, "\n", gene_a), paste0(og_b, "\n", gene_b)))
	p_tree_plotted <- p_tree +
		geom_point(
			data = df_node_data,
			aes(x = x, y = y, color = OGs_node, size = node_size),
			shape = 16,
			alpha = 0.5
		) +
		guides(size = "none") +
		scale_size_manual(values = c(0, node_plot_size/2, node_plot_size), na.translate = FALSE, drop = FALSE) +
		scale_color_manual(values = hcl.colors(2,"Blue-Red 2"), guide = "none") +
		new_scale_fill() +
		geom_fruit(
			data = df_tip_data %>% filter(og_presence == TRUE) %>% select(ID, OGs_tip),
			geom = geom_tile,
			mapping = aes(y = ID, x = OGs_tip, fill = OGs_tip),
			color = NA,
			pwidth = 0.4,
			offset = 0.12
		) +
		labs(fill = "OGs") +
		scale_fill_manual(values = hcl.colors(2,"Blue-Red 2")) +
		new_scale_fill() +
		geom_fruit(
			data = df_tax_data %>% rename(ID = "label") %>% select(ID, tax_color, genome_size),
			geom = geom_bar,
			mapping = aes(y = ID, x = genome_size, fill = tax_color),
			orientation = "y", 
                    	stat = "identity",
			color = NA,
			pwidth = 0.3,
			offset = 0.27
		) +
		geom_fruit(
			data = df_tax_data %>% rename(ID = "label") %>% select(ID, tax_color, gc_percentage),
			geom = geom_bar,
			mapping = aes(y = ID, x = gc_percentage, fill = tax_color),
			orientation = "y", 
                    	stat = "identity",
			color = NA,
			pwidth = 0.3,
			offset = 0.2
		) +
		#guides(fill = "none") +
		scale_fill_manual(values = cbp12, na.value = cbp12[12]) +
		labs(size = "Presence\nprobability", fill = tax_level) +
    		theme_void() +
  		layout_rectangular() +
		#coord_flip() +
		scale_y_reverse() +
    		theme(
			plot.background = element_rect(fill = "white",color = NA),
    			text = element_text(family = "Arial", size = 12),
    			axis.title.x = element_text(size = 13),
    			axis.title.y = element_text(size = 13),
    			plot.title = element_text(face = "bold", hjust = 0, size = 18),
    			plot.subtitle = element_text(hjust = 0.5, size = 12),
    			plot.caption = element_text(hjust = 0.5, size = 10),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13),
			legend.position = "right",
			legend.title = element_blank()
		)
}

contra_network_large_component_habi <- function (g_network_prokatlas, g_network_large_compo) {
	df_network_prokatlas <- g_network_prokatlas %>% as_data_frame() %>% as_tibble()
	ogs_m_f_s <- df_network_prokatlas %>% filter(to == "marine" | to == "freshwater" | to == "soil") %>% select(from) %>% pull() %>% unique()
	ogs_a_c_f <- df_network_prokatlas %>% filter(to == "anaerobic_digester" | to == "compost" | to == "fossil") %>% select(from) %>% pull() %>% unique()
	ogs_guts <- df_network_prokatlas %>% filter(to == "human_gut" | to == "chicken_gut" | to == "pig_gut" | to == "mouse_gut" | to == "human_bile") %>% select(from) %>% pull() %>% unique()
	ogs_rock <- df_network_prokatlas %>% filter(to == "rock_porewater") %>% select(from) %>% pull() %>% unique()
	df_ogs_habitat <- rbind(
		tibble(
			name = c(ogs_m_f_s),
			habitats = "marine;freshwater;soil"
		),
		tibble(
			name = c(ogs_a_c_f),
			habitats = "anaerobic_digester;compost;fossil"
		),
		tibble(
			name = c(ogs_guts),
			habitats = "Guts;human_bile"
		),
		tibble(
			name = c(ogs_rock),
			habitats = "rock_porewater"
		)
	)
	df_ogs_habitat <- left_join(
		g_contra_large_component_network %>% as_data_frame(what = "vertices") %>% as_tibble() %>% select(name),
		df_ogs_habitat,
		by = "name"
	)
	df_ogs_habitat$habitats[is.na(df_ogs_habitat$habitats)] <- "Others"
	V(g_contra_large_component_network)[df_ogs_habitat$name]$habitat <- df_ogs_habitat$habitats
	set.seed(1)
	df_g_large_component <- g_contra_large_component_network %>%
		ggnetwork(layout = igraph::layout_with_fr(g_contra_large_component_network)) %>%
		as_tibble() %>%
		mutate(
			x_same = x == xend,
			y_same = y == yend,
			v_e = map2(x_same, y_same, ~ if (.x == TRUE & .y == TRUE) {"vertex"} else {"edge"})
		) %>%
		unnest(c(v_e))
	df_g_large_component$habitat <- factor(df_g_large_component$habitat, levels = c("Guts;human_bile", "marine;freshwater;soil", "anaerobic_digester;compost;fossil", "rock_porewater", "Others"))
	df_g_large_component$n_replace_prok[df_g_large_component$n_replace_prok >= 100] <- 100
	p_network <- ggplot() +
		geom_segment(
			aes(x = x, y = y, xend = xend, yend = yend),
			size = 1,
			color = "grey70",
			data = df_g_large_component %>% filter(v_e == "edge")
		) +
		geom_point(
			aes(x = x, y = y, fill = habitat, size = betweenness),
			color = "white",
			shape = 21,
			data = df_g_large_component %>% filter(v_e == "vertex")
		) +
		scale_size_continuous(range = c(2,5)) +
		scale_fill_manual(values = c(cbp8[1:4], "grey60"), na.value = "black") +
		labs(fill = "Habitat group") +
		guides(size = "none") +
		#geom_text_repel(
		#	aes(x = x, y = y, label = plot_text),
		#	data = df_g_large_component %>% filter(v_e == "vertex", !is.na(vertex_color))
		#) +
		theme_blank() +
		theme(
			text = element_text(family = "Arial", size = 12),
			plot.title = element_text(hjust = 0.5, size = 13),
			plot.background = element_rect(fill = "white",color = NA),
			plot.tag = element_text(face = "bold", size = 18),
			strip.text.x = element_text(size = 13),
			strip.text.y = element_text(size = 13)
		)
	return(p_network)
}

distribution_function_degree <- function (g_graph, tax_cutoff, edgetype) {
	df_distribution_degree <- tibble(degree = igraph::degree(g_graph)) %>%
	  group_by(degree) %>%
	  summarize(
	    count = dplyr::n(),
	    total_nodes = length(V(g_graph))
	  ) %>%
	  mutate(
	    Density = count / total_nodes,
	    `Taxa cutoff` = tax_cutoff,
	    edge_type = edgetype
	  )
	return(df_distribution_degree)
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

cytoscape_anticorrelate_taxon <- function (g_contra_displace, node_color_fun, network_title) {
	createNetworkFromIgraph(g_contra_displace, network_title)
	setLayoutProperties('kamada-kawai', list(singlePartition=1, randomize=0))
	layoutNetwork('kamada-kawai')
	setNodeShapeDefault('ellipse')
	setNodeBorderWidthDefault(4)
	setNodeBorderColorDefault("#FFFAFA")
	setEdgeLineWidthMapping('n_tax', widths=c(10,40))
	cividis_pseudo <- function (n, alpha = 1, begin = 0, end = 0.9, direction = 1, option = "cividis") {
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
	    table.column = 'n_tax',
	    colors = cividis_pseudo,
	    mapping.type = "c"
	  )
	setNodeSizeMapping('degree', sizes=c(50,100))
	if (node_color_fun == "cog_fun_cat") {
		df_cog_node_color <- tibble(cog_fun_cat = V(g_contra_displace)$cog_fun_cat) %>%
		  distinct() %>% arrange(cog_fun_cat) %>%
		  mutate(
		    node_color = gsub("FF$","",viridis(length(unique((V(g_contra_displace)$cog_fun_cat)))))
		  )
		setNodeColorMapping(
		    table.column = 'cog_fun_cat',
		    table.column.values = unique(V(g_contra_displace)$cog_fun_cat),
		    colors = tibble(cog_fun_cat = unique(V(g_contra_displace)$cog_fun_cat)) %>%
		        left_join(df_cog_node_color, by = "cog_fun_cat") %>%
		        select(node_color) %>% pull(),
		    mapping.type='d'
		  )
	} else if (node_color_fun == "ko_path_num") {
		ko_path_nums <- V(g_contra_displace)$ko_path_num
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
}

read_network_data_each_taxon <- function (tax_x, edgetype, OG, g_whole) {
	#tax_x <- "p__Patescibacteria"
	#edgetype <- "Anticorrelating"
	#OG <- "cog"
	#g_whole <- g_cog_0_contra
	if (edgetype == "Anticorrelating" & OG == "cog") {
		df_network <- fread(paste0("data/contrapositive_pairs_cog2020/contrapositive_pairs_cog2020_", tax_x, ".tsv"), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("", "NULL")) %>%
		  as_tibble() %>%
		  mutate(
		    edge_type = "Anticorrelating"
		  )
	} else if (edgetype == "Correlating" & OG == "cog") {
		df_network <- fread(paste0("data/cooccurring_pairs_cog2020/cooccurring_pairs_cog2020_", tax_x, ".tsv"), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("", "NULL")) %>%
		  as_tibble() %>%
		  mutate(
		    edge_type = "Anticorrelating"
		  )
	} else if (edgetype == "Anticorrelating" & OG == "kegg") {
		df_network <- fread(paste0("data/contrapositive_pairs_ko/contrapositive_pairs_ko_", tax_x, ".tsv"), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("", "NULL")) %>%
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
		df_network <- fread(paste0("data/cooccurring_pairs_ko/cooccurring_pairs_ko_", tax_x, ".tsv"), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings=c("", "NULL")) %>%
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
	igraph::E(g_network)$n_tax[igraph::E(g_network)$n_tax >= 10] <- 10
	return(g_network)
}

subgraph_by_n_tax <- function (g_network, i) {
	g_subgraph <- subgraph.edges(g_network, as_ids(igraph::E(g_network)[igraph::E(g_network)$n_tax >= i]), delete.vertices = TRUE)
	return(g_subgraph)
}

cytoscape_to_ndex_anticorrelate <- function (tax_x, edgetype, OG, g_whole, node_color_fun, usr, password) {
	#tax_x <- "d__Bacteria"
	#edgetype <- "Anticorrelating"
	#OG <- "cog"
	#g_whole <- g_cog_0_contra
	#node_color_fun <- "cog_fun_cat"

	# n_tax_thre will be max of =< 1,000 edges
	g_og_contra_tax_x <- read_network_data_each_taxon(tax_x, edgetype, OG, g_whole)
	n_tax_thre <- 1
	while (length(igraph::E(g_og_contra_tax_x)[n_tax >= n_tax_thre]) > 1000) {
		n_tax_thre <- n_tax_thre + 1
	}
	net_title <- paste0("anticorrelating-", OG, "-", n_tax_thre, "-", tax_x)
	gsub_og_contra_tax_x <- subgraph_by_n_tax(g_og_contra_tax_x, n_tax_thre)
	V(gsub_og_contra_tax_x)$degree <- igraph::degree(gsub_og_contra_tax_x)
	cytoscape_anticorrelate_taxon(gsub_og_contra_tax_x, node_color_fun, net_title)
	ndex_uuid <- exportNetworkToNDEx(usr, password, FALSE)
	deleteAllNetworks()
	df_og_tax <- tibble(
		OGs = V(gsub_og_contra_tax_x)$name,
		taxon = tax_x
	)
	df_ndex_net <- tibble(
		Taxon = tax_x,
		OG = OG,
		`Network title (links)` = net_title,
		`n_tax threshold` = n_tax_thre,
		`# edges` = igraph::E(gsub_og_contra_tax_x) %>% length(),
		`# nodes` = igraph::V(gsub_og_contra_tax_x) %>% length(),
		`NDEx UUID` = ndex_uuid
	)
	return(list(df_og_tax, df_ndex_net))
}

ndex_network_degree <- function (tax_x, edgetype, OG, g_whole, node_color_fun) {
	#tax_x <- "d__Bacteria"
	#edgetype <- "Anticorrelating"
	#OG <- "cog"
	#g_whole <- g_cog_0_contra
	#node_color_fun <- "cog_fun_cat"

	# n_tax_thre will be max of =< 1,000 edges
	g_og_contra_tax_x <- read_network_data_each_taxon(tax_x, edgetype, OG, g_whole)
	n_tax_thre <- 1
	while (length(igraph::E(g_og_contra_tax_x)[n_tax >= n_tax_thre]) > 1000) {
		n_tax_thre <- n_tax_thre + 1
	}
	net_title <- paste0("anticorrelating-", OG, "-", n_tax_thre, "-", tax_x)
	gsub_og_contra_tax_x <- subgraph_by_n_tax(g_og_contra_tax_x, n_tax_thre)
	df_degree <- tibble(
		OGs = V(gsub_og_contra_tax_x)$name,
		taxon = tax_x,
		Degree = igraph::degree(gsub_og_contra_tax_x)
	)
	return(df_degree)
}
