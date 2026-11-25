source('scripts/r/functions.R')

args <- commandArgs(trailingOnly = T)
search_result_file <- args[1]
soft <- args[2]
out_path <- args[3]
#search_result_file <- "data/mmseqs_tsv_cog2020/proteins_raw_cog2020.0000.tsv"
#soft <- "mmseqs"
#out_path <- "data/mmseqs_tsv_filtered_cog2020/"

if (soft == "mmseqs") {
	out_file <- sub(".tsv", ".filtered.tsv", sub(".*/", "", search_result_file))
	#df_search_result <- mmseqs_search_result(search_result_file) %>% filter(tcov >= 0.8, qcov >= 0.8, bits >= 50)
	df_search_result <- mmseqs_search_result(search_result_file) %>% filter(tcov >= 0.8, bits >= 50)
} else if (soft == "diamond") {
	out_file <- sub(".tsv", ".filtered.tsv", sub(".*/", "", search_result_file))
	#df_search_result <- diamond_search_result(search_result_file) %>% filter(tcov >= 0.8, qcov >= 0.8, bits >= 50)
	df_search_result <- diamond_search_result(search_result_file) %>% filter(tcov >= 0.8, bits >= 50)
}

# Filter results by 0.8 * max(bits) in each range
df_search_result_filtered <- filter_mmseqs()

write.table(df_search_result_filtered, paste0(out_path, out_file), append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE)
