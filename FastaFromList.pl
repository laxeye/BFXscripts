#!/usr/bin/perl
use strict;
use warnings;

open FILE,$ARGV[0] or die $!; # List of proteins/genes/contigs, each in a new line
	my @query = <FILE>;
close FILE;

my %hash;

foreach(@query){
	chomp;
	$hash{$_} = 1;
}

open FILE,$ARGV[1] or die $!; # Multifasta file
	my @orig = <FILE>;
close FILE;
my $outfile = $ARGV[0] . "\.fna";
open FILE,">$outfile" or die $!; # Output

my $printnext=0;

foreach (@orig){
	if (/\>/) {
		chomp;
		$_ =~ s/\>//;
		if(exists($hash{$_})){
			print FILE ">$_\n";
			$printnext=1;
			next;
		}else{
			$printnext=0;
		}
	}
	next if $printnext==0;
	print FILE if $printnext == 1;
}

close FILE;
