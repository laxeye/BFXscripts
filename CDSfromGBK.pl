use warnings;
use strict;

open FILE, $ARGV[0] or die $!;
	my @file = <FILE>;
close FILE;

my $doprint = 0;

foreach(@file){
	if(/protein_id/){
		print(">" . (split/"/)[1] , "\n");
		$doprint = 1;
		next;
	}
	elsif($doprint eq 1){
		s/translation=//;
		if(/\"$/){
			s/\W//g;
			print "$_\n";
			$doprint = 0;
		}else{
		s/\W//g;
		print "$_\n";
		}
	}
}
