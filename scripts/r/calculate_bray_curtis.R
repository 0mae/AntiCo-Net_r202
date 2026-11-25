source('scripts/r/functions_ms_fisg_tables.R')

args <- commandArgs(trailingOnly = T)
OG <- args[1]
graph_tsv <- args[2]
out_tsv <- args[3]
#OG <- "cog"
#graph_tsv <- "data/tsv/anticorrelogome_network_cog_ntax8.tsv.gz"
#out_tsv <- "data/tsv/anticorrelogome_network_cog_ntax8_bray.tsv"

# Load network data (extract first two columns)
df_graph <- read.csv(graph_tsv, header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
  as_tibble()
df_edge <- df_graph[,c(1,2)]
colnames(df_edge) <- c("from","to")
# Load normalize GMGC data
if (OG == "cog") {
  df_gmgc_summary <- read.csv("data/tsv/gmgc10_cog_norm.tsv.gz", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
    as_tibble()
} else if (OG == "kegg") {
  df_gmgc_summary <- read.csv("data/tsv/gmgc10_kegg_norm.tsv.gz", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL")) %>%
    as_tibble()
}
# Make list (hash) for normalize GMGC data
l_gmgc_summary <- pmap(df_gmgc_summary %>% select(-OGs,-total_unigenes), list)
names(l_gmgc_summary) <- df_gmgc_summary$OGs
# Make function for calculating Bray-Curtis dissimilarity
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
# A few minutes step
df_network_bray <- df_edge %>%
  mutate(
    bray = map2(from, to, ~ 
                  # Skip OGs that are not found in eggnog (e.g. "COG4985")
                  if (is.null(l_gmgc_summary[[.x]]) || is.null(l_gmgc_summary[[.y]])) {
                    NA
                  # Skip OGs that have 0 abundance (e.g. "COG1851")
                  } else if (sum(is.nan(unlist(l_gmgc_summary[[.x]]))) > 0 || sum(is.nan(unlist(l_gmgc_summary[[.y]]))) > 0) {
                    NA
                  } else {
                    paired_ogs_bray(.x, .y)
                  })
  ) %>% unnest(c(bray))
# Write tsv
write.table(df_network_bray, out_tsv, append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE, na = "")
