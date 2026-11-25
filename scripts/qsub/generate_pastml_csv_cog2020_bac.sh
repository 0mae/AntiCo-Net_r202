source /etc/profile.d/modules.sh
module load R/4.1.2
cd $PBS_O_WORKDIR

Rscript --vanilla --slave scripts/r/generate_pastml_csv.R \
  data/proteins_raw_gtdbtax_cog2020_tidy/proteins_raw_gtdbtax_cog2020_tidy_d__Bacteria.tsv \
  data/pastml_csv_cog2020_bac/ \
  12
