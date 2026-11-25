source /etc/profile.d/modules.sh
module load seqkit/0.16.1
cd $PBS_O_WORKDIR

seqkit stats db/kegg/genes/fasta/prokaryotes.pep
seqkit stats db/kegg_fixed/prokaryotes.pep.with_ko.faa
