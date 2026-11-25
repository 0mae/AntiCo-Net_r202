source('scripts/r/functions.R')

df_ko_genes_list <- read_ko_genes_list()
#df_ko_genes_list <- df_ko_genes_list %>% head(100000)

df_geneid_koid <- df_ko_genes_list %>%
	group_by(gene_id) %>%
	summarize(ko_id = list(ko_id)) %>%
	rowwise() %>%
	mutate(
		is_multi_domain = length(ko_id) > 1,
		ko_ids = paste(unique(sort(ko_id)), collapse = ";")
	) %>%
	ungroup() %>%
	unnest(c(ko_id))

write.table(df_geneid_koid, "db/kegg_fixed/kegg_geneid_koid.tsv", append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE)
