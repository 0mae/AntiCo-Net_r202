library(tidyverse)
library(ape)

if (exists("params")) {
  dom <- params$dom
  og_a <- params$og_a
  og_b <- params$og_b
} else {
  args <- commandArgs(trailingOnly = TRUE)
  dom <- args[1]
  og_a <- args[2]
  og_b <- args[3]
  #dom <- "bacteria"
  #og_a <- "COG0207"
  #og_b <- "COG1351"
}

if (dom == "archaea") {
  Dom <- "ar"
  tree_path <- "ogs_distribution/species_ar.tree.gz"
} else if (dom == "bacteria") {
  Dom <- "bac"
  tree_path <- "ogs_distribution/species_bac.tree.gz"
}

untar_csv <- function(Dom, og_x) {
  OG_prefix <- gsub("[0-9]", "", og_x)
  if (OG_prefix == "COG"){
    og_tag <- "cog2020"
  } else if (OG_prefix == "K"){
    og_tag <- "ko"
  }
  tar_header <- paste0("pastml_csv_", og_tag, "_", Dom)
  untar(
      paste0("ogs_distribution/", tar_header, ".tar.gz"),
      paste0(tar_header, "/", og_x, ".csv"),
      exdir = "ogs_distribution"
      )
  df_og_x <- read.csv(paste0("ogs_distribution/", tar_header, "/", og_x, ".csv"), header = TRUE, quote = "", na.strings=c("", "NULL")) %>%
      as_tibble()
  system(paste0("rm -rf ogs_distribution/", tar_header))
  return(df_og_x)
}

read_genomes_analyzed <- function() {
	df_genomes_analyzed <- rbind(
		read.csv("ogs_distribution/genomes_analyzed_ar.tsv", header = TRUE, sep = "\t", quote = "") %>% as_tibble(),
		read.csv("ogs_distribution/genomes_analyzed_bac.tsv", header = TRUE, sep = "\t", quote = "") %>% as_tibble()
	)
	return(df_genomes_analyzed)
}

df_og_a <- untar_csv(Dom, og_a)
colnames(df_og_a) <- c("accession", "A")
df_og_b <- untar_csv(Dom, og_b)
colnames(df_og_b) <- c("accession", "B")
gtdb_node_labeled_tree <- read.tree(tree_path)
df_genomes <- read_genomes_analyzed() %>% select(accession,Domain:Species,gc_percentage,genome_size) %>% rename(node = "accession")

df_ogs <- left_join(df_og_a, df_og_b, by = "accession") %>% 
  mutate(OG_presence = map2(
    A, B, 
    ~ if (.x == TRUE & .y == TRUE) {
      "Both"
    } else if (.x == TRUE & .y == FALSE) {
      og_a
    } else if (.x == FALSE & .y == TRUE) {
      og_b
    } else {
      NA
    })
  ) %>% 
  unnest(OG_presence) %>% 
  select(accession, OG_presence, A, B)
colnames(df_ogs) <- c("accession", "OG presence", og_a, og_b)

df_taxonium_meta <- right_join(df_ogs, df_genomes %>% dplyr::rename(accession = "node"), by = "accession")

write.table(df_taxonium_meta, paste0("ogs_distribution/species_", Dom, "_", og_a, "_", og_b, ".tsv"), append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE)

system(paste0("gzip ogs_distribution/species_", Dom, "_", og_a, "_", og_b, ".tsv"))
