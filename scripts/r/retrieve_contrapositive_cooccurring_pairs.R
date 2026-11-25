source('scripts/r/functions.R')

args <- commandArgs(trailingOnly = T)
OG <- args[1]
in_file <- args[2]
out_path <- args[3]
og_pair_type <- args[4]

if (og_pair_type == "contrapositive") {
    df_contingency <- read.csv(in_file, header = TRUE, sep = "\t", quote = "") %>%
        as_tibble() %>%
        filter(
            contrapositive_ratio >= 0.8,
            ratio_present_A >= 0.1,
            ratio_present_A <= 0.9,
            ratio_present_B >= 0.1,
            ratio_present_B <= 0.9,
            TF / (TF + FT) >= 0.1,
            TF / (TF + FT) <= 0.9
        )
} else if (og_pair_type == "cooccurring") {
    df_contingency <- read.csv(in_file, header = TRUE, sep = "\t", quote = "") %>%
        as_tibble() %>%
        filter(
            cooccurring_ratio >= 0.8,
            ratio_present_A >= 0.1,
            ratio_present_A <= 0.9,
            ratio_present_B >= 0.1,
            ratio_present_B <= 0.9,
            TT / (TT + FF) >= 0.1,
            TT / (TT + FF) <= 0.9
        )
}

if (OG == "cog") {
    df_deff <- read_cog2020_def_fixed() %>% dplyr::rename(og_id = "cog_id")
} else if (OG == "kegg") {
    df_deff <- read_ko_def() %>% dplyr::rename(og_id = "ko_id")
}

colnames(df_deff) <- gsub("$", "_A", colnames(df_deff))

df_contingency <- df_contingency %>% left_join(df_deff, by = "og_id_A")

colnames(df_deff) <- gsub("_A$", "_B", colnames(df_deff))

df_contingency <- df_contingency %>% left_join(df_deff, by = "og_id_B")

write.table(df_contingency, out_path, append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE, na = "")
