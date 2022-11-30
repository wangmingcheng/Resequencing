#!/usr/bin/perl -w
use strict;

die "This version has bug, please use v1.1\n";

if (@ARGV != 2){
	print "Usage: perl $0 infile(gatk anno) outfile\n" ;
	exit(1);
}


my ($infile, $outfile) = @ARGV ;
open (IN, $infile) || die $! ;
my %hstat = ();
while (<IN>){
	chomp ;
	next if (m/^\#/ || m/^\s*$/);
	next if (m/INTERGENIC|INTRON|INTRAGENIC|DOWNSTREAM|UPSTREAM|UTR_3_PRIME|UTR_5_PRIME|SPLICE_SITE_ACCEPTOR|SPLICE_SITE_DONOR|START_GAINED|SYNONYMOUS_CODING|SYNONYMOUS_STOP|SYNONYMOUS_START/) ;
	if (m/SNPEFF_GENE_NAME=(.*?)\;/){
		$hstat{$1} ++ ;
	}
}
close(IN) ;

open (OUT, ">$outfile") || die $! ;
my $sum = keys %hstat ;
print OUT "#number\t$sum\n" ;
for my $key (keys %hstat){
	print OUT $key, "\n" ;
}
close(OUT) ;

