source /etc/profile.d/modules.sh
cd $PBS_O_WORKDIR

bash scripts/bash/generate_genomes_ogs_tidy_tsv.sh data/tsv/proteins_raw_gtdbtax_ko.tsv data/tsv/proteins_raw_gtdbtax_ko_tidy.tsv
