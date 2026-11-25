source('scripts/r/functions.R')

args <- commandArgs(trailingOnly = TRUE)
in_file <- args[1]
out_path <- args[2]
n_cpu <- as.numeric(args[3])
#in_file <- "data/proteins_raw_gtdbtax_cog2020_tidy/proteins_raw_gtdbtax_cog2020_tidy_d__Archaea.tsv"
#out_path <- "data/pastml_csv_cog2020_ar/"

if (length(grep("cog", in_file)) >= 1) {
	dt_genomes_ogs_tidy <- read.csv(in_file, header = TRUE, sep = "\t", quote = "") %>%
		dplyr::rename(og_id = "cog_id") %>%
		mutate(presence = count > 0) %>%
		data.table() %>%
		lazy_dt()
} else if (length(grep("ko", in_file)) >= 1) {
	dt_genomes_ogs_tidy <- read.csv(in_file, header = TRUE, sep = "\t", quote = "") %>%
		mutate(
			og_id = map(ko_id, ~ sub("ko:", "", .x)),
			presence = count > 0
		) %>%
		select(-ko_id) %>%
		unnest(c(og_id)) %>%
		data.table() %>%
		lazy_dt()
}

df_pivot <- dt_genomes_ogs_tidy %>% select(-tax,-count) %>%
	pivot_wider(names_from = "og_id", values_from = presence, values_fill = FALSE) %>%
	as_tibble()

write_pastml_csv <- function (each_og) {
	df_out <- df_pivot %>% select(accession, .data[[each_og]])
	write.table(
		df_out,
		paste0(out_path, each_og, ".csv"),
		append=FALSE, quote=FALSE, sep = ",", row.names=FALSE
	)
}

mclapply(colnames(df_pivot)[2:dim(df_pivot)[2]], write_pastml_csv, mc.cores = n_cpu)
