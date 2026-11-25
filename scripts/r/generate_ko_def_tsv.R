source('scripts/r/functions.R')

args <- commandArgs(trailingOnly = T)
out_path <- args[1]

df_koid_kosymbol_koname <- read.csv("db/kegg_fixed/kegg_koid_kosymbol_koname.tsv", header = TRUE, sep = "\t", quote = "") %>% as_tibble()

df_koid_cog <- read_ko_lists("db/kegg/genes/ko/ko_cog.list", c("ko_id", "ko_cog")) %>%
    group_by(ko_id) %>% summarize(ko_cog = paste(ko_cog, collapse = ";"))
df_koid_ec <- read_ko_lists("db/kegg/genes/ko/ko_enzyme.list", c("ko_id", "ko_ec")) %>%
    group_by(ko_id) %>% summarize(ko_ec = paste(ko_ec, collapse = ";"))
df_koid_path <- read_ko_lists("db/kegg/genes/ko/ko_pathway.list", c("ko_id", "ko_path")) %>%
    group_by(ko_id) %>% summarize(ko_path = paste(ko_path, collapse = ";"))
df_koid_rn <- read_ko_lists("db/kegg/genes/ko/ko_reaction.list", c("ko_id", "ko_rn")) %>%
    group_by(ko_id) %>% summarize(ko_rn = paste(ko_rn, collapse = ";"))
df_koid_br <- read_ko_lists("db/kegg/genes/ko/ko_brite.list", c("ko_id", "ko_br")) %>%
    group_by(ko_id) %>% summarize(ko_br = paste(ko_br, collapse = ";"))

df_ko_def <- df_koid_kosymbol_koname %>%
    left_join(df_koid_cog, by = "ko_id") %>%
    left_join(df_koid_ec, by = "ko_id") %>%
    left_join(df_koid_path, by = "ko_id") %>%
    left_join(df_koid_rn, by = "ko_id") %>%
    left_join(df_koid_br, by = "ko_id") %>%
	mutate(
		ko_paths_num = map(ko_path, ~ unique(sort(sub("path:ko", "", sub("path:map", "", unlist(strsplit(.x, ";"))))))),
		ko_paths_num_specific = map(ko_paths_num, ~ setdiff(.x, c("01100", "01110", "01120", "01200", "01210", "01212", "01230", "01240", "01220"))),
		ko_path_num = map(ko_paths_num, ~ if (length(.x) == 0) {NA} else {paste0(.x, collapse = ";")}),
		ko_path_num_specific = map(ko_paths_num_specific, ~ if (length(.x) == 0) {NA} else {paste0(.x, collapse = ";")}),
        ko_brs = map(ko_br, ~ unique(sort(unlist(strsplit(.x, ";"))))),
        ko_brs_specific = map(ko_brs, ~ setdiff(.x, c("br:ko00001"))),
        ko_br_specific = map(ko_brs_specific, ~ if (length(.x) == 0) {NA} else {paste0(.x, collapse = ";")})
	) %>%
	select(-ko_paths_num, -ko_paths_num_specific, -ko_brs, -ko_brs_specific) %>%
	unnest(c(ko_path_num, ko_path_num_specific, ko_br_specific))


write.table(df_ko_def, out_path, append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE, na = "")
