seqkit fx2tab db/kegg/genes/fasta/prokaryotes.pep | \
  perl -ne 's/ .+?\t/\t/;print' | \
  sort -k 1,1 | \
  join -t "$(printf '\011')" -1 1 -2 1 -a 1 <(bash scripts/bash/kegg_genes_prokaryotes_with_ko.sh) - | \
  seqkit tab2fx \
  > db/kegg_fixed/prokaryotes.pep.with_ko.faa
