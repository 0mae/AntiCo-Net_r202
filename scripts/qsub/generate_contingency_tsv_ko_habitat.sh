source /etc/profile.d/modules.sh
module load R/4.1.2
cd $PBS_O_WORKDIR

Rscript --vanilla --slave scripts/r/generate_contingency_tsv_habitat.R kegg k__Prokaryotes data/tsv/contingency_tsv_ko_habitat_hp5.tsv 5
Rscript --vanilla --slave scripts/r/generate_contingency_tsv_habitat.R kegg k__Prokaryotes data/tsv/contingency_tsv_ko_habitat_hp0.tsv 0
Rscript --vanilla --slave scripts/r/generate_contingency_tsv_habitat.R kegg k__Prokaryotes data/tsv/contingency_tsv_ko_habitat_hp10.tsv 10
