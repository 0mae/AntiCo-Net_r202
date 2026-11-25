source /etc/profile.d/modules.sh
module load R/4.1.2
cd $PBS_O_WORKDIR

Rscript --vanilla --slave scripts/r/generate_ko_def_tsv.R db/kegg_fixed/kegg_ko_def.tsv
