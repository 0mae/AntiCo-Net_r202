source /etc/profile.d/modules.sh
module load pastml/1.9.34
cd $PBS_O_WORKDIR

ls data/pastml_csv_cog2020_bac/ | sed 's/.csv$//' | xargs -I {} -P 40 -n 1 sh -c "bash scripts/bash/run_pastml.sh cog bacteria {}"
