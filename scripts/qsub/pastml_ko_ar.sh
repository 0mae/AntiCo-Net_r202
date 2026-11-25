source /etc/profile.d/modules.sh
module load pastml/1.9.34
cd $PBS_O_WORKDIR

ls data/pastml_csv_ko_ar/ | sed 's/.csv$//' | xargs -I {} -P 12 -n 1 sh -c "bash scripts/bash/run_pastml.sh kegg archaea {}"
