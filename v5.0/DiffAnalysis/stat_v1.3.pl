#!/usr/bin/perl -w
use strict;

if (@ARGV != 3){
	print "\nVersion: v1.3\tselect mRNA or gene\n"; 
	print "Usage: perl $0 infile(gatk anno) outfile mRNA/gene\n" ;
	exit(1);
}

my ($infile, $outfile,$type) = @ARGV ;
my %hstat = ();
if($type eq "mRNA"){
	open IN,$infile or die $!;
	open OUT,">$outfile.nonsy" or die $!;
	while(<IN>){
		chomp;
		next if /^$/;
		print OUT "$_\n" if /^\#/;
		next if /INTERGENIC|INTRON|INTRAGENIC|DOWNSTREAM|UPSTREAM|UTR_3_PRIME|UTR_5_PRIME|SNPEFF_EFFECT=SPLICE_SITE_ACCEPTOR|SNPEFF_EFFECT=SPLICE_SITE_DONOR|SNPEFF_EFFECT=START_GAINED|SNPEFF_EFFECT=SYNONYMOUS_CODING|SNPEFF_EFFECT=SYNONYMOUS_STOP|SNPEFF_EFFECT=SYNONYMOUS_START/;
		if(/SNPEFF_TRANSCRIPT_ID=(.*?)\;/ || /SNPEFF_TRANSCRIPT_ID=(.*?)\s/){
			$hstat{$1}++;
			print OUT "$_\n";
		}
	}
	close IN;
	close OUT;
}

if($type eq "gene"){
	open IN,$infile or die $!;
	open OUT,">$outfile.nonsy" or die $!;
	while(<IN>){
		chomp;
		next if /^$/;
		print OUT "$_\n" if /^\#/;
		next if /INTERGENIC|INTRON|INTRAGENIC|DOWNSTREAM|UPSTREAM|UTR_3_PRIME|UTR_5_PRIME|SNPEFF_EFFECT=SPLICE_SITE_ACCEPTOR|SNPEFF_EFFECT=SPLICE_SITE_DONOR|SNPEFF_EFFECT=START_GAINED|SNPEFF_EFFECT=SYNONYMOUS_CODING|SNPEFF_EFFECT=SYNONYMOUS_STOP|SNPEFF_EFFECT=SYNONYMOUS_START/;
		if(/SNPEFF_GENE_NAME=(.*?)\;/ || /SNPEFF_GENE_NAME=(.*?)\s/){
			$hstat{$1}++;
			print OUT "$_\n";
		}
	}
	close IN;
	close OUT;
}

open OUT, ">$outfile" or die $!;
my $sum = keys %hstat ;
print OUT "#number\t$sum\n" ;
for my $key (keys %hstat){
	print OUT $key, "\n" ;
}
close(OUT) ;

