#!/usr/bin/perl

use strict;

my @in = <STDIN>;

my @array = split /\/\/\//, (join('', @in));

pop(@array);

print "ko_id\tko_symbol\tko_name\n";

foreach(@array){
	my @str = split /\n/;
	my $ko;
	my $sym;
	my $name;
	foreach(@str){
		if(m/^ENTRY/){
			s/^ENTRY +//;
			s/ .*//;
			$ko = "ko:$_";
		}elsif(m/^SYMBOL/){
			s/^SYMBOL +//;
			$sym = $_;
		}elsif(m/^NAME/){
			s/^NAME +//;
			$name = $_;
		}
	}
	print "$ko\t$sym\t$name\n";
}
