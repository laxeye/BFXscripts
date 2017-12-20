#!/usr/bin/perl
use strict;
use warnings;

local $/ = ">";

usage() if $#ARGV != 0;

open F,"$ARGV[0]" || die;
	my @data = <F>;
close F;

shift @data;

my ($title, $seq, $length);
my @nucl = ("A", "C", "G", "T");
my @tetramers;
my @result;

foreach my $i (@nucl){
	foreach my $j (@nucl){
		foreach my $k (@nucl){
			push @tetramers,$i.$j.$k.$_ foreach(@nucl)
		}
	}
}

print "#SequenceName";

foreach(@tetramers){
	print "\t$_";
}
print "\n";

my $count;

foreach(@data) {
	@result = ();
	($title, $seq) = split /\n/,$_,2;
	$seq =~ s/\n//;
	$length = length($seq);
	foreach(@tetramers){
		$count = () = ($seq =~ /$_/gi);
		push @result,sprintf "%.4f", 100*$count/$length;
	}
	print "$title";
	foreach(@result){
		print "\t$_";
	}
	print "\n";	
}

sub usage{
	print "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n";
	print "\#Tetranucleotide statistics for (meta)genomic DNA sequences \#\n";
	print "\#\#\#\#\#\#\# Usage: tetramers.pl <input.fasta> > outut.csv \#\#\#\#\#\#\#\n";
	print "\#\#\#\#\#\# The output is tab-separated and includes header \#\#\#\#\#\#\n";
	print "\#\#\#\#\#\#\#\#\#\#\# The values are presented as percents \#\#\#\#\#\#\#\#\#\#\#\#\n";
	print "\#\#\#\#\#\#\#\#\#\#\# Licensed under The MIT License (MIT) \#\#\#\#\#\#\#\#\#\#\#\#\n";
	print "\#\#\#\#\#\# Copyright (c) Aleksei A. Korzhenkov, 2016-2017 \#\#\#\#\#\#\#\n";
	print "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n";
	exit;
}