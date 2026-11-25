source /etc/profile.d/modules.sh
cd $PBS_O_WORKDIR

bash scripts/bash/generate_koid_kosymbol_koname_tsv.sh
