source('scripts/r/functions.R')

args <- commandArgs(trailingOnly = TRUE)
OG <- args[1]
dom <- args[2]
n_cpu <- as.numeric(args[3])
#OG <- "cog"
#dom <- "archaea"
#n_cpu <- 8
#check COG0207 & COG1351 (x <- 62067)
#check COG0436 (x <- 61547)

df_net_edge_stats <- fread(paste0("data/tsv/network_edge_stats_", OG, ".tsv"), header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
	as_tibble() %>%
	mutate(edge_id = paste(og_id_A, og_id_B, sep = "_"))

#df_net_edge_stats_cooccu <- edge_stats4og_replacement_count_cooccu(df_net_edge_stats)
#df_net_edge_stats_contra <- edge_stats4og_replacement_count_contra(df_net_edge_stats)

#df_net_edge_stats_select <- rbind(df_net_edge_stats_cooccu, df_net_edge_stats_contra) %>% select(og_id_A, og_id_B) %>% distinct()

df_net_edge_stats_select <- df_net_edge_stats %>% filter(edge_type == "contrapositive")

df_og_replacement <- do.call(dplyr::bind_rows, mclapply(1:dim(df_net_edge_stats_select)[1], count_og_replacement, mc.cores = n_cpu))

df_og_replacement_arrange <- df_og_replacement %>% select(
	og_id_A,
	og_id_B,
	tree_dom,
	n_replace,
	n_A_N,
	n_A_A,
	n_A_B,
	n_A_AB,
	n_A_U,
	n_B_N,
	n_B_A,
	n_B_B,
	n_B_AB,
	n_B_U,
	n_AB_N,
	n_AB_A,
	n_AB_B,
	n_AB_AB,
	n_AB_U,
	n_N_N,
	n_N_A,
	n_N_B,
	n_N_AB,
	n_N_U,
	tree_edges_A_B,
	tree_edges_B_A
)

write.table(df_og_replacement_arrange, paste0("data/tsv/og_displacement_", dom, "_", OG, ".tsv"), append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE, na = "")
