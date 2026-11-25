source /etc/profile.d/modules.sh
cd $PBS_O_WORKDIR

grep '^>' db/kegg/genes/fasta/prokaryotes.pep | wc -l
zcat db/kegg/genes/fasta/prokaryotes.dat.gz | cut -f 1 | sort | uniq | wc -l
bash scripts/bash/kegg_genes_prokaryotes_with_ko.sh | wc -l
bash scripts/bash/kegg_prokaryotes_ko.sh | wc -l
cut -f 2 db/kegg/genes/ko/ko_genes.list | sort | uniq | wc -l
cut -f 1 db/kegg/genes/ko/ko_genes.list | sort | uniq | wc -l
zcat db/kegg/genes/links/genes_ko.list.gz | cut -f 1 | sort | uniq | wc -l
zcat db/kegg/genes/links/genes_ko.list.gz | cut -f 2 | sort | uniq | wc -l
