if [ $1 = "cog" ] && [ $2 = "archaea" ]; then

	pastml \
		--prediction_method MPPA \
		--tree raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree \
		--data data/pastml_csv_cog2020_ar/$3.csv \
		-o data/pastml_cog2020_ar/$3.txt \
		--work_dir data/pastml_cog2020_ar/$3_pastml \
		--data_sep , \
		--tip_size_threshold 10000000000 \
		--threads 1 \
		2>qsub_out/pastml_cog2020_ar/$3.err && \
		targz -p 1 --rm data/pastml_cog2020_ar/$3_pastml

elif [ $1 = "cog" ] && [ $2 = "bacteria" ]; then

	pastml \
		--prediction_method MPPA \
		--tree raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree \
		--data data/pastml_csv_cog2020_bac/$3.csv \
		-o data/pastml_cog2020_bac/$3.txt \
		--work_dir data/pastml_cog2020_bac/$3_pastml \
		--data_sep , \
		--tip_size_threshold 10000000000 \
		--threads 1 \
		2>qsub_out/pastml_cog2020_bac/$3.err && \
		targz -p 1 --rm data/pastml_cog2020_bac/$3_pastml

elif [ $1 = "kegg" ] && [ $2 = "archaea" ]; then

	pastml \
		--prediction_method MPPA \
		--tree raw/gtdb/ar122_r202_ete3fmt5_node_labeled.tree \
		--data data/pastml_csv_ko_ar/$3.csv \
		-o data/pastml_ko_ar/$3.txt \
		--work_dir data/pastml_ko_ar/$3_pastml \
		--data_sep , \
		--tip_size_threshold 10000000000 \
		--threads 1 \
		2>qsub_out/pastml_ko_ar/$3.err && \
		targz -p 1 --rm data/pastml_ko_ar/$3_pastml

elif [ $1 = "kegg" ] && [ $2 = "bacteria" ]; then

	pastml \
		--prediction_method MPPA \
		--tree raw/gtdb/bac120_r202_ete3fmt5_node_labeled.tree \
		--data data/pastml_csv_ko_bac/$3.csv \
		-o data/pastml_ko_bac/$3.txt \
		--work_dir data/pastml_ko_bac/$3_pastml \
		--data_sep , \
		--tip_size_threshold 10000000000 \
		--threads 1 \
		2>qsub_out/pastml_ko_bac/$3.err && \
		targz -p 1 --rm data/pastml_ko_bac/$3_pastml

fi
