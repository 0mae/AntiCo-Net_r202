cat db/kegg/genes/ko/ko | \
  perl scripts/perl/parse_ko.pl | \
  sort -k 1,1 \
  > db/kegg_fixed/kegg_koid_kosymbol_koname.tsv
