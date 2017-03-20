#!/usr/bin/perl
# Written by Aleksei Korzhenkov (2017) 
# using "Building Customized Data Pipelines Using the Entrez Programming Utilities (eUtils)" 
# by Eric Sayers and David Wheeler 
# (https://www.ncbi.nlm.nih.gov/books/NBK1058/).
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.

use warnings;
use LWP::Simple;
use Net::FTP;

my $host = "ftp.ncbi.nlm.nih.gov";
my $user = "anonymous";
my $password = 'anonymous@mail.com';

if (not defined $ARGV[0] or $ARGV[0] =~ /-h/){
	print "Usage $0 -p|-n|-gff|-gbk <list of accesion numbers>.\n";
	print "List should have one accesion number in each line.\n";
	exit;
}
my $ftype;
$ftype = "p" if $ARGV[0] =~ /-p/;
$ftype = "n" if $ARGV[0] =~ /-n/;
$ftype = "gff" if $ARGV[0] =~ /-gff/;
$ftype = "gbk" if $ARGV[0] =~ /-gbk/;


open FILE,"$ARGV[1]" or die $!;
my @queries = <FILE>;	# Reading file
close FILE;
my $N = $#queries+1;
my $i=1;

my $f = Net::FTP->new($host,Timeout => 30, Passive => 1) or die "Can't open $host\n"; # Connecting ftp-server
$f->login($user, $password) or die "Can't log $user in\n"; # Logging in
$f->binary(); # Binary mode
print "Succesfully connected to NCBI!\n";
print "$N files will be downloaded.\n";



my ($ftpath,$dir,$file_to_get);
my @files;
my $db = 'assembly'; # What DB we are looking in
my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'; # EUtils URL

foreach my $query (@queries){
	$file_to_get ="0";

	#assemble the esearch URL
	$url = $base . "esearch.fcgi?db=$db&term=$query&usehistory=y";

	#post the esearch URL
	$output = get($url);

	#parse WebEnv and QueryKey
	$web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
	$key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);

	### include this code for ESearch-ESummary
	#assemble the esummary URL
	$url = $base . "esummary.fcgi?db=$db&query_key=$key&WebEnv=$web";

	#post the esummary URL
	$docsums = get($url);
	
	$ftpath = $1 if ($docsums =~ /<FtpPath_GenBank>(\S+)<\/FtpPath_GenBank>/);

	$ftpath =~ s!ftp://ftp.ncbi.nlm.nih.gov!!; # Cut off the server name

	$f->cwd($ftpath) or die "Can't cwd to $ftpath\n";

	local $, = "\n";

	@files = $f->ls() or die "Can't ls in $ftpath!\n";
	foreach(@files){
		if ($ftype =~ /p/ and /protein.faa.gz/){ 
			# Fasta AA CDS
			$file_to_get = $_;
			last;
		}
		if ($ftype =~ /n/ and /genomic/ && !/_rna_from_genomic/ && !/_cds_from_genomic/){ 
			# Fasta AA CDS
			$file_to_get = $_;
			last;
		}
		if ($ftype =~ /gff/ and /genomic.gff/){ 
			# Fasta AA CDS
			$file_to_get = $_;
			last;
		}
		if ($ftype =~ /gbk/ and /genomic.gb/){ 
			# Fasta AA CDS
			$file_to_get = $_;
			last;
		}
	}

	if ($file_to_get eq 0){
		print "Failed to download genome sequence for $query!\n";
		$i ++;
		next;
	}

	$f->get($file_to_get) or die "Can't get $file_to_get from $ftpath\n";
	print "Downloading $file_to_get. $i of $N files.\n";

	$i++;
}
print "Download complete!\n";