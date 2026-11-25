ls data/contrapositive_pairs_cog2020/* |\
  xargs -I {} sh -c "tail -n +2 {}" |\
  cut -f 1-14 |\
  cat <(head -n1 data/contrapositive_pairs_cog2020/contrapositive_pairs_cog2020_c__ABY1.tsv | cut -f 1-14) - \
  > data/tsv/contrapositive_pairs_cog2020.tsv
ls data/contrapositive_pairs_ko/* |\
  xargs -I {} sh -c "tail -n +2 {}" |\
  cut -f 1-14 |\
  cat <(head -n1 data/contrapositive_pairs_ko/contrapositive_pairs_ko_c__ABY1.tsv | cut -f 1-14) - \
  > data/tsv/contrapositive_pairs_ko.tsv
