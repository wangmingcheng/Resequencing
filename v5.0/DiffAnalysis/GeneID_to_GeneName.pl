#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);

if (@ARGV != 2){
	print "Usage: perl $0 infile(anno.gatk.*.list) gff\n" ;
	exit(1);
}

my ($infile,$gff) = @ARGV;
$infile = abs_path($infile);
$gff = abs_path($gff);

open GFF,"$gff" or die $!;
my %hash;
while(<GFF>){
	chomp;
	next if /^$/;
	next if /^\#/;
	my ($gene,$info) = (split /\t/)[2,8];
	next if $gene ne 'gene';
	if($info =~ m/ID=(.*?);.*?Name=([0-9_\.A-Za-z]+);?/){
		$hash{$2} = $1;
	}
}
close GFF;

open LIST,"$infile" or die $!;
open OUT , ">${infile}2" or die $!;

my %gene;
while(<LIST>){
	chomp;
	next if /^$/;
	if (/^\#/){
		print OUT "$_\n"; 
		next;
	}else{
		print OUT $hash{$_},"\n" if(exists $hash{$_});
	}
}
close LIST;
close OUT;


