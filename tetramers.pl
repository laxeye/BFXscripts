#!/usr/bin/perl
use strict;
use warnings;

local $/ = ">";

usage() if not defined $ARGV[0];

open F,"$ARGV[0]" || die;
	my @data = <F>;
close F;

shift @data;

my $title;
my $seq;
my $length;
my @nucl = ("A", "C", "G", "T");
my @bimers;
my @result;
my @names;
my %hashtable;

foreach my $oldest (@nucl){
foreach my $older (@nucl){
	foreach my $old (@nucl){
		foreach(@nucl){
			push @bimers,$oldest.$older.$old.$_;
		}
	}
}
}

foreach(@bimers){
	print "\t$_";
}
print "\n";

foreach(@data) {
	@result = "";
	($title, $seq) = split /\n/,$_,2;
	push @names,$title;
	$seq =~ s/\n//;
	$length = length($seq) + 1;
	foreach(@bimers){
		push @result,sprintf "%.4f", 100*($seq =~ s/$_//gi)/$length;
	}
#	$C = sprintf "%.3f", $seq =~ tr/C///$length;
#	$G = sprintf "%.3f", $seq =~ tr/G///$length;
#	$GC = $G + $C;
	print "$title";
	foreach(@result){
		print "\t$_";
	}
	print "\n";
	
}

sub usage{
	print "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n";
	print "\#\#\#\#\#\#\#Tetranucleotide statistics for (meta)genomic DNA\#\#\#\#\#\#\n";
	print "\#\#\#\#\#\#\#Usage: tetramers.pl <input.fasta> > outut.csv\#\#\#\#\#\#\#\#\#\n";
	print "\#\#\#The output is tab-separated and the 2nd column is empty\#\#\#\n";
	print "\#\#\#\#\#\#\#\#\#\#\#\#Licensed under The MIT License (MIT)\#\#\#\#\#\#\#\#\#\#\#\#\#\n";
	print "\#\#\#\#\#\#\#\#\#Copyright (c) Aleksei A. Korzhenkov, 2016\#\#\#\#\#\#\#\#\#\#\#\n";
	print "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n";
	exit;
}