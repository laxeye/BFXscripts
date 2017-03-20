#!/usr/bin/env perl

#Usage: samtools depth -a [any options for samtools] <mapping.file.sorted.bam> | check_mg_coverage.pl
#Or: check_mg_coverage.pl < <file made with samtools>

use warnings;

my ($n,$p,$d);
my $oldname ="";
my @node=(0);
my $cont = 0;
my $m = 0;
my $depth = 0;

my $fold = 4;	#Coverage difference to report
my $window = 50; #Window size
my $slip = 10;	#Step size

while(<>){
	chomp;
	($n,$p,$d) = split/\t/;
	if ($n !~ /$oldname/){
		my @sorted = sort {$a <=> $b} @node;
		my $median = $sorted[int($#sorted/2)];
		my $i = 1;
		while($i < $#sorted - $window){
			foreach $j (0..$window-1){
				$depth += $node[$i+$j];
			}
			$depth = $depth/$window;
			if ($depth > $fold * $median || $depth < $median / $fold){
				if ($cont == 1){
					#Exception continues
					$m = $node[$i] if $m < $node[$i];
				}else{
					#New exception
					$start = $i;
					print "$oldname\t$median\t$node[$i]\t$i\n";
					$cont = 1;
				}
			}else{
				if ($cont == 1){
					$prv = $i - 1;
					if ($prv != $start){
						print "$oldname\t$median\t$m\tmax\n";
						print ("$oldname\t$median\t$node[$prv]\t$prv\n");
					}
					$cont = 0;
					$m = 0;
				}
			}
			$i+=$slip;
		}
		if ($cont == 1){
			print "$oldname\t$median\t$m\tmax\n";
			print "$oldname\t$median\t$node[$i]\t$i\n";
		}
		$m = 0;
		$cont = 0;
		@node=(0);
	}
	push @node,$d;
	$oldname = $n;
}

