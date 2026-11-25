source('scripts/r/functions.R')

args <- commandArgs(trailingOnly = T)
OG <- args[1]
#OG <- "cog"

df_taxa_summarized_50 <- read.csv("data/tsv/taxa_summarized_50_stats_pd.tsv", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>% as_tibble()


if (OG == "cog") {
	df_contra_pairs <- fread("data/tsv/contrapositive_pairs_cog2020.tsv", stringsAsFactors=FALSE,header = TRUE, sep = "\t", quote = "") %>%
		select(og_id_A, og_id_B, tax) %>%
		mutate(edge_type = "contrapositive")
	df_cooccu_pairs <- fread("data/tsv/cooccurring_pairs_cog2020.tsv", stringsAsFactors=FALSE,header = TRUE, sep = "\t", quote = "") %>%
		select(og_id_A, og_id_B, tax) %>%
		mutate(edge_type = "cooccurring")
} else if (OG == "kegg") {
	df_contra_pairs <- fread("data/tsv/contrapositive_pairs_ko.tsv", stringsAsFactors=FALSE,header = TRUE, sep = "\t", quote = "") %>%
		select(og_id_A, og_id_B, tax) %>%
		mutate(edge_type = "contrapositive")
	df_cooccu_pairs <- fread("data/tsv/cooccurring_pairs_ko.tsv", stringsAsFactors=FALSE,header = TRUE, sep = "\t", quote = "") %>%
		select(og_id_A, og_id_B, tax) %>%
		mutate(edge_type = "cooccurring")
}

df_graphs <- rbind(
		as_tibble(df_contra_pairs),
		as_tibble(df_cooccu_pairs)
	) %>%
	nest(data = c(-tax, -edge_type)) %>%
	mutate(
		graphs = map(data, ~ graph_from_data_frame(., directed = FALSE) %>% simplify(edge.attr.comb = "concat"))
	)

# Graph stats

df_graph_structure <- df_graphs %>% mutate(
		edge_density = map(graphs, ~ edge_density(.)),
		n_node = map(graphs, ~ vcount(.)),
		n_edge = map(graphs, ~ ecount(.)),
		assortativity_degree = map(graphs, ~ assortativity_degree(., directed = FALSE)),
		mean_distance = map(graphs, ~ mean_distance(., directed = FALSE)),
		n_components = map(graphs, ~ count_components(.)),
		max_component_size = map(graphs, ~ max(components(.)$csize)),
		n_triangles = map(graphs, ~ sum(count_triangles(.)))
	) %>%
	select(-data, -graphs) %>%
	unnest(c(
		edge_density,
		n_node,
		n_edge,
		assortativity_degree,
		mean_distance,
		n_components,
		max_component_size,
		n_triangles
	)) %>%
	left_join(df_taxa_summarized_50, by = "tax")

df_contra_graph <- df_contra_pairs %>%
	as_tibble() %>%
	left_join(df_taxa_summarized_50, by = "tax") %>%
	select(og_id_A, og_id_B, tax, faith_pd) %>%
	graph_from_data_frame(directed = FALSE) %>%
	simplify(edge.attr.comb = "concat") %>%
	as_data_frame() %>%
	as_tibble() %>%
	mutate(
		n_tax = map(tax, ~ length(.)),
		faith_pd_sum = map(faith_pd, ~ sum(.))
	) %>%
	unnest(c(n_tax, faith_pd_sum)) %>%
	select(-tax, -faith_pd) %>%
	mutate(edge_type = "contrapositive")

df_cooccu_graph <- df_cooccu_pairs %>%
	as_tibble() %>%
	left_join(df_taxa_summarized_50, by = "tax") %>%
	select(og_id_A, og_id_B, tax, faith_pd) %>%
	graph_from_data_frame(directed = FALSE) %>%
	simplify(edge.attr.comb = "concat") %>%
	as_data_frame() %>%
	as_tibble() %>%
	mutate(
		n_tax = map(tax, ~ length(.)),
		faith_pd_sum = map(faith_pd, ~ sum(.))
	) %>%
	unnest(c(n_tax, faith_pd_sum)) %>%
	select(-tax, -faith_pd) %>%
	mutate(edge_type = "cooccurring")

df_graph_structure_prok <- tibble(
		tax = c("k__Prokaryotes", "k__Prokaryotes"),
		edge_type = c("contrapositive", "cooccurring"),
		data = list(df_contra_graph, df_cooccu_graph)
	) %>%
	mutate(
		graphs =  map(data, ~ graph_from_data_frame(., directed = FALSE)),
		edge_density = map(graphs, ~ edge_density(.)),
		n_node = map(graphs, ~ vcount(.)),
		n_edge = map(graphs, ~ ecount(.)),
		assortativity_degree = map(graphs, ~ assortativity_degree(., directed = FALSE)),
		mean_distance = map(graphs, ~ mean_distance(., directed = FALSE)),
		n_components = map(graphs, ~ count_components(.)),
		max_component_size = map(graphs, ~ max(components(.)$csize)),
		n_triangles = map(graphs, ~ sum(count_triangles(.)))
	) %>%
	select(-data, -graphs) %>%
	unnest(c(edge_density, n_node, n_edge, assortativity_degree, mean_distance, n_components, max_component_size, n_triangles)) %>%
	left_join(df_taxa_summarized_50, by = "tax")

write.table(
		dplyr::bind_rows(df_graph_structure_prok, df_graph_structure),
		paste0("data/tsv/network_structure_", OG,".tsv"),
		append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE, na = ""
	)

# Edge stats

df_whole_graph <- dplyr::bind_rows(df_contra_graph, df_cooccu_graph) %>%
	graph_from_data_frame(directed = FALSE) %>%
	simplify(edge.attr.comb = "concat") %>%
	as_data_frame() %>%
	as_tibble() %>%
	mutate(
		n_edge_type = map(edge_type, ~ length(.))
	) %>%
	unnest(c(-from, -to)) %>%
	arrange(desc(n_tax))

if (OG == "cog") {
	#df_ogs_def <- read_cog2020_def_na() %>% dplyr::rename(og_id = "cog_id")
	df_ogs_def <- read_cog2020_def_fixed() %>% dplyr::rename(og_id = "cog_id")
} else if (OG == "kegg") {
	df_ogs_def <- read_ko_def_na() %>% dplyr::rename(og_id = "ko_id")
}

colnames(df_ogs_def) <- gsub("$", "_A", colnames(df_ogs_def))

df_whole_graph_def <- df_whole_graph %>%
	dplyr::rename(og_id_A = "from", og_id_B = "to") %>%
	left_join(df_ogs_def, by = "og_id_A")

colnames(df_ogs_def) <- gsub("_A$", "_B", colnames(df_ogs_def))

df_whole_graph_def <- df_whole_graph_def %>%
	left_join(df_ogs_def, by = "og_id_B")

if (OG == "cog") {

	df_whole_graph_def_tfcol <- df_whole_graph_def %>% mutate(
			cog_fun_cats_A = strsplit(cog_fun_cat_A, ""),
			n_cog_fun_cats_A = map(cog_fun_cats_A, ~ length(.)),
			cog_fun_cats_B = strsplit(cog_fun_cat_B, ""),
			n_cog_fun_cats_B = map(cog_fun_cats_B, ~ length(.)),
			share_cog_fun_cat = map2(cog_fun_cats_A, cog_fun_cats_B, ~ length(intersect(.x, .y)) >= 1),
			share_cog_fun_path = map2(cog_fun_path_A, cog_fun_path_B, ~ if (is.na(.x[[1]]) | is.na(.y[[1]])) {FALSE} else {length(intersect(.x, .y)) >= 1})
		) %>%
		select(-cog_fun_cats_A, -cog_fun_cats_B) %>%
		unnest(c(n_cog_fun_cats_A, n_cog_fun_cats_B, share_cog_fun_cat, share_cog_fun_path))

} else if (OG == "kegg") {

	df_whole_graph_def_tfcol <- df_whole_graph_def %>% mutate(
			ko_cogs_A = strsplit(ko_cog_A, ";"),
			ko_cogs_B = strsplit(ko_cog_B, ";"),
			n_ko_cogs_A = map(ko_cogs_A, ~ length(.)),
			n_ko_cogs_B = map(ko_cogs_B, ~ length(.)),
			share_ko_cog = map2(ko_cogs_A, ko_cogs_B, ~ if (is.na(.x[[1]]) | is.na(.y[[1]])) {FALSE} else {length(intersect(.x, .y)) >= 1}),
			ko_ecs_A = strsplit(ko_ec_A, ";"),
			ko_ecs_B = strsplit(ko_ec_B, ";"),
			n_ko_ecs_A = map(ko_ecs_A, ~ length(.)),
			n_ko_ecs_B = map(ko_ecs_B, ~ length(.)),
			share_ko_ec = map2(ko_ecs_A, ko_ecs_B, ~ if (is.na(.x[[1]]) | is.na(.y[[1]])) {FALSE} else {length(intersect(.x, .y)) >= 1}),
			ko_ec3s_A = map(ko_ecs_A, ~ sub("\\.[0-9]*$", "", .x)),
			ko_ec3s_B = map(ko_ecs_B, ~ sub("\\.[0-9]*$", "", .x)),
			share_ko_ec3 = map2(ko_ec3s_A, ko_ec3s_B, ~ if (is.na(.x[[1]]) | is.na(.y[[1]])) {FALSE} else {length(intersect(.x, .y)) >= 1}),
			ko_ec2s_A = map(ko_ec3s_A, ~ sub("\\.[0-9]*$", "", .x)),
			ko_ec2s_B = map(ko_ec3s_B, ~ sub("\\.[0-9]*$", "", .x)),
			share_ko_ec2 = map2(ko_ec2s_A, ko_ec2s_B, ~ if (is.na(.x[[1]]) | is.na(.y[[1]])) {FALSE} else {length(intersect(.x, .y)) >= 1}),
			ko_ec1s_A = map(ko_ec2s_A, ~ sub("\\.[0-9]*$", "", .x)),
			ko_ec1s_B = map(ko_ec2s_B, ~ sub("\\.[0-9]*$", "", .x)),
			share_ko_ec1 = map2(ko_ec1s_A, ko_ec1s_B, ~ if (is.na(.x[[1]]) | is.na(.y[[1]])) {FALSE} else {length(intersect(.x, .y)) >= 1}),
			ko_rns_A = strsplit(ko_rn_A, ";"),
			ko_rns_B = strsplit(ko_rn_B, ";"),
			n_ko_rns_A = map(ko_rns_A, ~ length(.)),
			n_ko_rns_B = map(ko_rns_B, ~ length(.)),
			share_ko_rn = map2(ko_rns_A, ko_rns_B, ~ if (is.na(.x[[1]]) | is.na(.y[[1]])) {FALSE} else {length(intersect(.x, .y)) >= 1}),
			ko_paths_num_specific_A = strsplit(ko_path_num_specific_A, ";"),
			ko_paths_num_specific_B = strsplit(ko_path_num_specific_B, ";"),
			n_ko_paths_num_specific_A = map(ko_paths_num_specific_A, ~ length(.)),
			n_ko_paths_num_specific_B = map(ko_paths_num_specific_B, ~ length(.)),
			share_ko_path_num_specific = map2(ko_paths_num_specific_A, ko_paths_num_specific_B, ~ if (is.na(.x[[1]]) | is.na(.y[[1]])) {FALSE} else {length(intersect(.x, .y)) >= 1}),
			ko_brs_specific_A = strsplit(ko_br_specific_A, ";"),
			ko_brs_specific_B = strsplit(ko_br_specific_B, ";"),
			n_ko_brs_specific_A = map(ko_brs_specific_A, ~ length(.)),
			n_ko_brs_specific_B = map(ko_brs_specific_B, ~ length(.)),
			share_ko_br_specific = map2(ko_brs_specific_A, ko_brs_specific_B, ~ if (is.na(.x[[1]]) | is.na(.y[[1]])) {FALSE} else {length(intersect(.x, .y)) >= 1}),
		) %>%
		select(-ko_cogs_A, -ko_cogs_B, -ko_ecs_A, -ko_ecs_B, -ko_ec3s_A, -ko_ec3s_B, -ko_ec2s_A, -ko_ec2s_B, -ko_ec1s_A, -ko_ec1s_B, -ko_rns_A, -ko_rns_B, -ko_paths_num_specific_A, -ko_paths_num_specific_B, -ko_brs_specific_A, -ko_brs_specific_B) %>%
		unnest(c(n_ko_cogs_A:share_ko_br_specific))

}

write.table(
		df_whole_graph_def_tfcol,
		paste0("data/tsv/network_edge_stats_", OG, ".tsv"),
		append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE, na = ""
	)

# Degree of nodes

colnames(df_ogs_def) <- gsub("_B$", "", colnames(df_ogs_def))

df_degree <- df_graphs %>%
	mutate(
		degree = map(graphs, ~ degree(.)),
		og_id = map(degree, ~ names(.))
	) %>%
	select(-data, -graphs) %>%
	unnest(c(degree, og_id))

prok_degree_contra <- df_contra_graph %>%
	graph_from_data_frame(directed = FALSE) %>%
	degree()

prok_degree_cooccu <- df_cooccu_graph %>%
	graph_from_data_frame(directed = FALSE) %>%
	degree()

df_degree_whole <- rbind(
	tibble(
		tax = "k__Prokaryotes",
		edge_type = "contrapositive",
		degree = prok_degree_contra,
		og_id = names(prok_degree_contra)
	),
	tibble(
		tax = "k__Prokaryotes",
		edge_type = "cooccurring",
		degree = prok_degree_cooccu,
		og_id = names(prok_degree_cooccu)
	),
	df_degree
)

write.table(
		df_degree_whole %>%
			left_join(df_ogs_def, by = "og_id") %>%
			left_join(df_taxa_summarized_50, by = "tax"),
		paste0("data/tsv/network_degree_", OG, ".tsv"),
		append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE, na = ""
	)

# Distance

df_distance <- df_graphs %>%
	mutate(
		distances = map(graphs, ~ distance_table(., directed = FALSE)),
		n_edge = map(distances, ~ c(.$res, .$unconnected)),
		distance = map(distances, ~ c(1:length(.$res), Inf))
	) %>%
	select(-data, -graphs, -distances) %>%
	unnest(c(distance, n_edge))

prok_distance_contra <- df_contra_graph %>%
	graph_from_data_frame(directed = FALSE) %>%
	distance_table()

prok_distance_cooccu <- df_cooccu_graph %>%
	graph_from_data_frame(directed = FALSE) %>%
	distance_table()

df_distance_whole <- rbind(
	tibble(
		tax = "k__Prokaryotes",
		edge_type = "contrapositive",
		n_edge = c(prok_distance_contra$res, prok_distance_contra$unconnected),
		distance = c(1:length(prok_distance_contra$res), Inf)
	),
	tibble(
		tax = "k__Prokaryotes",
		edge_type = "cooccurring",
		n_edge = c(prok_distance_cooccu$res, prok_distance_cooccu$unconnected),
		distance = c(1:length(prok_distance_cooccu$res), Inf)
	),
	df_distance
)

write.table(
		left_join(df_distance_whole, df_taxa_summarized_50, by = "tax"),
		paste0("data/tsv/network_distance_", OG, ".tsv"),
		append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE, na = ""
	)

# Betweenness

df_betweenness <- df_graphs %>%
	mutate(
		betweenness = map(graphs, ~ betweenness(., directed = FALSE)),
		og_id = map(betweenness, ~ names(.))
	) %>%
	select(-data, -graphs) %>%
	unnest(c(betweenness, og_id))

prok_betweenness_contra <- df_contra_graph %>%
	graph_from_data_frame(directed = FALSE) %>%
	betweenness()

prok_betweenness_cooccu <- df_cooccu_graph %>%
	graph_from_data_frame(directed = FALSE) %>%
	betweenness()

df_betweenness_whole <- rbind(
	tibble(
		tax = "k__Prokaryotes",
		edge_type = "contrapositive",
		betweenness = prok_betweenness_contra,
		og_id = names(prok_betweenness_contra)
	),
	tibble(
		tax = "k__Prokaryotes",
		edge_type = "cooccurring",
		betweenness = prok_betweenness_cooccu,
		og_id = names(prok_betweenness_cooccu)
	),
	df_betweenness
)

write.table(
		df_betweenness_whole %>%
			left_join(df_ogs_def, by = "og_id") %>%
			left_join(df_taxa_summarized_50, by = "tax"),
		paste0("data/tsv/network_betweenness_", OG, ".tsv"),
		append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE, na = ""
	)

# Clustering coefficient

df_clustering_coef <- df_graphs %>%
	mutate(
		clustering_coef = map(graphs, ~ transitivity(., type = "local")),
		og_id = map(graphs, ~ V(.)$name)
	) %>%
	select(-data, -graphs) %>%
	unnest(c(clustering_coef, og_id))

prok_clustering_coef_contra <- df_contra_graph %>%
	graph_from_data_frame(directed = FALSE) %>%
	transitivity(type = "local")

prok_clustering_coef_cooccu <- df_cooccu_graph %>%
	graph_from_data_frame(directed = FALSE) %>%
	transitivity(type = "local")

df_clustering_coef_whole <- rbind(
	tibble(
		tax = "k__Prokaryotes",
		edge_type = "contrapositive",
		clustering_coef = prok_clustering_coef_contra,
		og_id = V(df_contra_graph %>% graph_from_data_frame(directed = FALSE))$name
	),
	tibble(
		tax = "k__Prokaryotes",
		edge_type = "cooccurring",
		clustering_coef = prok_clustering_coef_cooccu,
		og_id = V(df_cooccu_graph %>% graph_from_data_frame(directed = FALSE))$name
	),
	df_clustering_coef
)

write.table(
		df_clustering_coef_whole %>%
			left_join(df_ogs_def, by = "og_id") %>%
			left_join(df_taxa_summarized_50, by = "tax"),
		paste0("data/tsv/network_clustering_coef_", OG, ".tsv"),
		append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE, na = ""
	)

# Connected components

df_components <- df_graphs %>%
	mutate(
		membership = map(graphs, ~ components(.)$membership),
		og_id = map(membership, ~ names(.))
	) %>%
	select(-data, -graphs) %>%
	unnest(c(membership, og_id))

prok_components_contra <- df_contra_graph %>%
	graph_from_data_frame(directed = FALSE) %>%
	components()

prok_components_cooccu <- df_cooccu_graph %>%
	graph_from_data_frame(directed = FALSE) %>%
	components()

df_components_whole <- rbind(
	tibble(
		tax = "k__Prokaryotes",
		edge_type = "contrapositive",
		membership = prok_components_contra$membership,
		og_id = names(prok_components_contra$membership)
	),
	tibble(
		tax = "k__Prokaryotes",
		edge_type = "cooccurring",
		membership = prok_components_cooccu$membership,
		og_id = names(prok_components_cooccu$membership)
	),
	df_components
)

write.table(
		df_components_whole %>%
			left_join(df_ogs_def, by = "og_id") %>%
			left_join(df_taxa_summarized_50, by = "tax"),
		paste0("data/tsv/network_components_", OG, ".tsv"),
		append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE, na = ""
	)

# Motif (3 nodes)

#df_motif_3 <- df_graphs %>%
#	mutate(
#		motif_3 = map(graphs, ~ motifs(., size = 3)),
#		isomorphism_class = map(motif_3, ~ 1:length(.))
#	) %>%
#	select(-data, -graphs) %>%
#	unnest(c(motif_3, isomorphism_class))
#
#prok_motif_3_contra <- df_contra_graph %>%
#	graph_from_data_frame(directed = FALSE) %>%
#	motifs(size = 3)
#
#prok_motif_3_cooccu <- df_cooccu_graph %>%
#	graph_from_data_frame(directed = FALSE) %>%
#	motifs(size = 3)
#
#df_motif_3_whole <- rbind(
#	tibble(
#		tax = "k__Prokaryotes",
#		edge_type = "contrapositive",
#		motif_3 = prok_motif_3_contra,
#		isomorphism_class = 1:length(prok_motif_3_contra)
#	),
#	tibble(
#		tax = "k__Prokaryotes",
#		edge_type = "cooccurring",
#		motif_3 = prok_motif_3_cooccu,
#		isomorphism_class = 1:length(prok_motif_3_cooccu)
#	),
#	df_motif_3
#)
#
#write.table(
#		df_motif_3_whole %>%
#			left_join(df_taxa_summarized_50, by = "tax"),
#		paste0("data/tsv/network_motif_3_", OG, ".tsv"),
#		append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE, na = ""
#	)

# Motif (4 nodes)

#df_motif_4 <- df_graphs %>%
#	mutate(
#		motif_4 = map(graphs, ~ motifs(., size = 4)),
#		isomorphism_class = map(motif_4, ~ 1:length(.))
#	) %>%
#	select(-data, -graphs) %>%
#	unnest(c(motif_4, isomorphism_class))
#
#prok_motif_4_contra <- df_contra_graph %>%
#	graph_from_data_frame(directed = FALSE) %>%
#	motifs(size = 4)
#
#prok_motif_4_cooccu <- df_cooccu_graph %>%
#	graph_from_data_frame(directed = FALSE) %>%
#	motifs(size = 4)
#
#df_motif_4_whole <- rbind(
#	tibble(
#		tax = "k__Prokaryotes",
#		edge_type = "contrapositive",
#		motif_4 = prok_motif_4_contra,
#		isomorphism_class = 1:length(prok_motif_4_contra)
#	),
#	tibble(
#		tax = "k__Prokaryotes",
#		edge_type = "cooccurring",
#		motif_4 = prok_motif_4_cooccu,
#		isomorphism_class = 1:length(prok_motif_4_cooccu)
#	),
#	df_motif_4
#)
#
#write.table(
#		df_motif_4_whole %>%
#			left_join(df_taxa_summarized_50, by = "tax"),
#		paste0("data/tsv/network_motif_4_", OG, ".tsv"),
#		append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE, na = ""
#	)
