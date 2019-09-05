#!/usr/bin/perl -w
# Written by Aleksei Korzhenkov (2019)
# Korzhenkov_AA@nrcki.ru
# using "Building Customized Data Pipelines Using the Entrez Programming Utilities (eUtils)" 
# by Eric Sayers and David Wheeler 
# (https://www.ncbi.nlm.nih.gov/books/NBK1058/).
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.

use warnings;
use LWP::Simple;

if (not defined $ARGV[0] or $ARGV[0] =~ /-h/){
	print "$0 - substitute NCBI accession identifiers to taxa in phylogenetic trees.\n";
	print "Usage $0 -aacc|-anam|-prot|-nucl <input.nwk> > output.nwk\n";
	print "\t-aacc Assembly accession\n";
	print "\t-anam Assembly name\n";
	print "\t-prot Protein accession\n";
	print "\t-nucl Nucleotide accesion.\n";
	exit;
}

# Parsing args
my $ftype;
$ftype = "AACC" if $ARGV[0] =~ /-aacc/;
$ftype = "ANAM" if $ARGV[0] =~ /-anam/;
$ftype = "PACC" if $ARGV[0] =~ /-prot/;
$ftype = "NACC" if $ARGV[0] =~ /-nucl/;

# Line separator undef to parse multiline trees
undef $/;

open FILE,"$ARGV[1]" or die $!;

my $file = <FILE>;
my @queries;

if ($ftype =~ "AACC") {
	@queries = ( $file =~ /(GC\w_\d+\.\d)/g );
} elsif($ftype =~ "ANAM") {
	@queries = ( $file =~ /(ASM\d+v\d)/g );
} else {
	# Protein and nucleotide accessions
	@queries = ( $file =~ /(\w+\d+\.\d)/g );
}

close FILE;

my ($ftpath,$dir,$file_to_get);
my @files;

# Setting DB
my $db = 'assembly';
if ($ftype =~ /PACC/) {
	$db = 'protein';
} elsif($ftype =~ /NACC/){
	$db = 'nucleotide';
}

# EUtils URL
my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';

foreach my $query (@queries){
	$file_to_get ="0";

	# Assemble and get the esearch URL
	$url = $base . "esearch.fcgi?db=$db&term=$query";
	$output = get($url);

	# Parse ID
	next unless ($output =~ /<Id>(\d+)<\/Id>/);
	$id = $1;

	# Parse Taxonomy ID
	$url = $base . "esummary.fcgi?db=$db&id=$id";
	$output = get($url);
	
	if ($db =~ /assembly/) {
		$org = $1 if ($output =~ /Taxid>(\d+)<\/Taxid>/);
	} else {
		$org = $1 if ($output =~ /TaxId.+>(\d+)</);
	}

	$url = $base . "esummary.fcgi?db=taxonomy&id=$org";
	$output = get($url);

	if ($output =~ /<Item Name="ScientificName".+>(.+)<\/Item>/){
		$org = "$1";
		$org =~ s/ /_/g;
		$org =~ s/[\[\]]//g;
	}

	#Print substitutions to STDERR
	print STDERR "Leaf: $query\tTaxon: $org\n";
	
	my $subst = $org."_".$query;
	$file =~ s/[\w\d\.]*$query[\w\d\.]*/$subst/;
}

print $file;
