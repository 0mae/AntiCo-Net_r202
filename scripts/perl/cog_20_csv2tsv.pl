#!/usr/bin/perl

use strict;

## USAGE: perl scripts/perl/cog_20_csv2tsv.pl db/cog2020_fixed/cog-20.wrong_names.seq_csv.tsv db/cog2020/cog-20.cog.csv db/cog2020_fixed/cog-20.wrong_columns.geneid_ncol.tsv > db/cog2020_fixed/cog-20.cog.tsv

my %hash;

open IN, $ARGV[0] or die "Cannot open $ARGV[0]: $!";
while(<IN>){
	chomp;
	my @array = split(/\t/, $_, -1);
	$hash{$array[1]} = $array[0];
}
close IN;

my %wrong_column;

open IN, $ARGV[2] or die "Cannot open $ARGV[2]: $!";
while(<IN>){
	chomp;
	my @array = split(/\t/, $_, -1);
	$wrong_column{$array[0]} = $array[1];
}
close IN;

print "gene_id\tncbi_assembly_id\tprotein_id\tprotein_length\tcog_footprint_on_protein\tcog_footprint_length\tcog_id\treserved\tcog_membership\tbitscore\te_val\tcog_profile_length\tprotein_footprint_on_cog\tseq_name\n";

open IN, $ARGV[1] or die "Cannot open $ARGV[1]: $!";
while(<IN>){
	chomp;
	my @array = split(/,/, $_, -1);
	if ($wrong_column{$array[0]} == 14) {
		pop @array;
		my $str = join("\t", @array);
		if (length($hash{$array[2]}) > 0) {
			print "$str\t$hash{$array[2]}\n";
		} else {
			print "$str\t$array[2]\n";
		}
	} else {
		my $str = join("\t", @array);
		if (length($hash{$array[2]}) > 0) {
			print "$str\t$hash{$array[2]}\n";
		} else {
			print "$str\t$array[2]\n";
		}
	}
}
close IN;



