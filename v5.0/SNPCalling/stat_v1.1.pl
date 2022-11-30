#!/usr/bin/perl -w
use strict;

if (@ARGV != 2){
	print "\nVersion: v1.1\tfixed bug of statistic\n"; 
	print "Usage: perl $0 infile(gatk anno) outfile\n" ;
	exit(1);
}

my ($infile, $outfile) = @ARGV ;
open (IN, $infile) || die $! ;
open (OUT, ">$outfile.nonsy") || die "Can't creat $outfile.nonsy" ;
my %hstat = ();
while (<IN>){
	chomp ;
	next if (m/^\s*$/);
	if (m/^\#/){
		print OUT $_,"\n" ;
	}
	next if (m/INTERGENIC|INTRON|INTRAGENIC|DOWNSTREAM|UPSTREAM|UTR_3_PRIME|UTR_5_PRIME|SNPEFF_EFFECT=SPLICE_SITE_ACCEPTOR|SNPEFF_EFFECT=SPLICE_SITE_DONOR|SNPEFF_EFFECT=START_GAINED|SNPEFF_EFFECT=SYNONYMOUS_CODING|SNPEFF_EFFECT=SYNONYMOUS_STOP|SNPEFF_EFFECT=SYNONYMOUS_START/) ;
	if (m/SNPEFF_GENE_NAME=(.*?)\;/){
		 (my $name = $1) =~ s/Transcript_(.*?).cds\d+$/$1/e ;
		$hstat{$name} ++ ;
		print OUT $_,"\n" ;
	}
}
close(IN) ;
close(OUT) ;

open (OUT, ">$outfile") || die $! ;
my $sum = keys %hstat ;
print OUT "#number\t$sum\n" ;
for my $key (keys %hstat){
	print OUT $key, "\n" ;
}
close(OUT) ;

