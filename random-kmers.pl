#!/usr/bin/perl
use warnings;
use strict;

my $string;
my @list;
my @alpha = ("A","C","G","T");
my ($i,$j);

unless (defined $ARGV[0] && defined $ARGV[1] && defined $ARGV[2] ){
	print "Plese run as random-kmers.pl <kmers count> <kmer length> <AT count in k-mer>\n";
	exit;
}

my $limit = $ARGV[0];
my $klength = $ARGV[1];
my $ATcount = $ARGV[2];

if ($klength < $ATcount){
	print "Error: k-mer length less than AT-count!\n";
	exit;
}

for($j = 1; $j <= $limit;){
	for($i = 1; $i <= $klength; $i++){
		$string .= $alpha[int(rand 4)]
	}
	if (($string =~ tr/A// + $string =~ tr/T//) eq $ATcount){
		push @list,$string;
		$j++;
	}
	$string ="";
}

$j=1;

foreach(@list){
	print ">$j\n";
	print $_ . "\n";
	$j++;
}