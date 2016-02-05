use strict;
use warnings;

open FILE,$ARGV[0] or die $!; # List of proteins, each in a new line
	my @query = <FILE>;
close FILE;

open FILE,$ARGV[1] or die $!; # Multifasta file
	my @orig = <FILE>;
close FILE;

open FILE,'>output.fa' or die $!; # Output file

foreach my $prot (@query){
	chomp $prot;
	my $printnext=0;
	foreach (@orig){
		if(/$prot$/){
			print FILE;
			$printnext=1;
			next;
		}
		last if /\>/ and $printnext == 1;
		print FILE if $printnext == 1;
	}
}

close FILE;
