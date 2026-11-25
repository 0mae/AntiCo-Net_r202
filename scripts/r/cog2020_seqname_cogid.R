source('scripts/r/functions.R')

df_cog2020 <- read_cog2020_tsv()

df_cog2020_seqname_cogid <- df_cog2020 %>% distinct(seq_name, cog_id, cog_membership) %>% group_by(seq_name, cog_id) %>% summarize(cog_memberships = paste(sort(list(cog_membership)[[1]]), collapse = ";"), cog_membership_rep = sort(list(cog_membership)[[1]])[1])

df_cog2020_n_cogs <- df_cog2020 %>% distinct(seq_name, cog_id) %>% group_by(seq_name) %>% summarize(n_cogs = dplyr::n())

df_cog2020_multi_domain <- df_cog2020_n_cogs %>% filter(n_cogs > 1) %>% select(seq_name) %>% mutate(is_multi_domain = TRUE)

df_cog2020_seqname_cogid_multi <- left_join(df_cog2020_seqname_cogid, df_cog2020_multi_domain, by = "seq_name")

df_cog2020_seqname_cogid_multi$is_multi_domain[is.na(df_cog2020_seqname_cogid_multi$is_multi_domain)] <- FALSE

df_cog2020_cogids <- df_cog2020_seqname_cogid_multi %>% group_by(seq_name, is_multi_domain) %>% summarize(cog_ids = paste(sort(list(cog_id)[[1]]), collapse = ";"))

#filter(df_cog2020_0_1_cogids, seq_name == "WP_128698900_1")$cog_ids

df_cog2020_seqname_cogid_multi_cogids <- inner_join(df_cog2020_seqname_cogid_multi, select(df_cog2020_cogids, -is_multi_domain), by = "seq_name") %>% arrange(seq_name)

write.table(df_cog2020_seqname_cogid_multi_cogids, "db/cog2020_fixed/cog-20.seqname_cogid.tsv", append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE)
