source /etc/profile.d/modules.sh
module load diamond/2.0.11
cd $PBS_O_WORKDIR

diamond makedb --in db/diamond_kegg/prokaryotes.pep.with_ko.faa -d db/diamond_kegg/prokaryotes.pep.with_ko -p 12
