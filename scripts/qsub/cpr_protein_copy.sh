source /etc/profile.d/modules.sh
cd $PBS_O_WORKDIR

source scripts/bash/functions.sh
cat data/gtdb/cpr_metadata.tsv | gtdbacc2path | xargs -I {} -P 8 -n 1 sh -c "cp raw/gtdb/gtdb_genomes_reps_r202/{}*_genomic.fna.gz data/gtdb/cpr_genomes_reps/"
