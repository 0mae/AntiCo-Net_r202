source /etc/profile.d/modules.sh
module load R/4.1.2
cd $PBS_O_WORKDIR

Rscript --vanilla --slave scripts/r/generate_pastml_csv.R \
  data/proteins_raw_gtdbtax_ko_tidy/proteins_raw_gtdbtax_ko_tidy_d__Archaea.tsv \
  data/pastml_csv_ko_ar/ \
  12
