source /etc/profile.d/modules.sh
module load CheckM/1.1.3
cd $PBS_O_WORKDIR

checkm qa --tab_table raw/gtdb/cpr_43_markers.hmm data/gtdb/cpr_checkm > data/gtdb/checkm_result_cpr.tsv
