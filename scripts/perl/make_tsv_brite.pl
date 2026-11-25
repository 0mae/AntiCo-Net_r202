#!/usr/bin/perl

use strict;

## USAGE: perl scripts/perl/make_tsv_brite.pl db/kegg_fixed/kegg_brite_ko.tsv > db/kegg_fixed/ko_brite_text.tsv

my $A;
my $B;
my $C;
my $D;

print "A\tB\tC\tko_id\tD\n";

open IN, $ARGV[0] or die "Cannot open $ARGV[0]: $!";
while(<IN>){
	chomp;
	my @array = split(/\t/, $_);
	if ($array[0] eq "A") {
		$A = $array[1];
	} elsif ($array[0] eq "B") {
		$B = $array[1];
	} elsif ($array[0] eq "C") {
		$C = $array[1];
	} elsif ($array[0] eq "D") {
		$D = $array[1];
		my $ko = $array[1];
		$ko =~ s/(K\d+) .*/$1/;
		print "$A\t$B\t$C\t$ko\t$D\n";
	}
}
close IN;
