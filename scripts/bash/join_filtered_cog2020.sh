if [[ $3 == "cog" ]]; then
	cat $2 |\
		tail -n +2 |\
		sort -k 8,8 |\
		join -t "$(printf '\011')" -1 8 -2 1 -a 1 - <(tail -n +2 $1 | sort -k 1,1) |\
		cat <(echo -ne "target\tquery\tranges.x.start\tranges.x.end\tranges.x.width\tranges.y.start\tranges.y.end\tranges.y.width\tfident\talnlen\tmismatch\tgapopen\ttstart\ttend\tevalue\tbits\tqlen\tqcov\ttlen\ttcov\tis_multi_domain_2\tcog_id\tcog_memberships\tcog_membership_rep\tis_multi_domain\tcog_ids\n") -
elif [[ $3 == "kegg" ]]; then
	cat $2 |\
		tail -n +2 |\
		sort -k 8,8 |\
		join -t "$(printf '\011')" -1 8 -2 1 -a 1 - <(tail -n +2 $1 | sort -k 1,1) |\
		cat <(echo -ne "target\tquery\tranges.x.start\tranges.x.end\tranges.x.width\tranges.y.start\tranges.y.end\tranges.y.width\tfident\talnlen\tmismatch\tgapopen\ttstart\ttend\tevalue\tbits\tqlen\tqcov\ttlen\ttcov\tis_multi_domain_2\tko_id\tis_multi_domain\tko_ids\n") -
fi
