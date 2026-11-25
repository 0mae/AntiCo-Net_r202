#!/usr/bin/perl

use strict;

## USAGE: perl scripts/perl/transform_kegg_pathway_list.pl db/kegg/pathway/pathway.list > db/kegg_fixed/pathway_mod.list

open IN, $ARGV[0] or die "Cannot open $ARGV[0]";
my @in = <IN>;
close IN;

my $col_first;
my $col_second;

foreach(@in) {
	chomp;
	if (m/^#([^#].+)/) {
		$col_first = $1;
	} elsif (m/^##([^#].+)/) {
		$col_second = $1;
	} else {
		print "$_\t$col_first\t$col_second\n";
	}
}
