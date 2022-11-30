#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

# ------------------------------------------------------------------
my ($fIn,$fInn,$fInnn,$fOut,$type);

GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"region:s"=>\$fIn,
				"gff:s"=>\$fInn,
				"anno:s"=>\$fInnn,
				"type:s" => \$type,
				) or &USAGE;
&USAGE unless ($fIn and $fInn and $fInnn and $fOut);
# ------------------------------------------------------------------
open (IN,$fIn) or die $!;      # 打开.dedup.realn.bam_CNVs文件
open (INN,$fInn) or die $!;    # 打开参考基因组genes.gff文件
open (INNN,$fInnn) or die $!;  # 打开参考基因组注释文件Integrated_Function.annotation.xls
open (OUT,">$fOut") or die $!;

my %hash1;
my %hash2;
my %hash3;
my %hash4;
my @all;
$type ||= "gene";
while (<IN>) {   #获取染色体、起始位置、终止位置
	chomp;
	my @line=split/\t/,$_;
	$hash1{$_}=$line[0];
	$hash2{$_}=$line[1];
	$hash3{$_}=$line[2];
	push @all,$_;
}

while(<INN>){
	chomp;
	next if /\#/;
	my @line = split /\t+/,$_;
	if($type eq "gene"){
		next if $line[2] ne "gene";
		foreach my $all(@all){
			if ((($hash1{$all} eq $line[0]) || ("chr".$hash1{$all} eq $line[0])) and (($line[3] >=$hash2{$all} and $line[3] <=$hash3{$all} )or ($line[4] >=$hash2{$all} and $line[4] <= $hash3{$all}))){
				if(/ID=(.*?)\;/){
					$hash4{$1} = $all;
					$hash4{$1} = $all if /Name=(.*?)\;/;
				}
			}else{
				next;
			}
		}
	}elsif($type eq "mRNA"){
		next if $line[2] ne "mRNA";
		foreach my $all(@all){
			if ((($hash1{$all} eq $line[0]) || ("chr".$hash1{$all} eq $line[0])) and (($line[3] >=$hash2{$all} and $line[3] <=$hash3{$all} )or ($line[4] >=$hash2{$all} and $line[4] <= $hash3{$all}))){
				if(/ID=(.*?)\;/){
					$hash4{$1} = $all;
					$hash4{$1} = $all if /Name = (.*?)\;/;
				}else{
					next;
				}
			}
		}	
	}else{
		next;
	}
}

while (<INNN>) {
	chomp;
	my @line=split/\t+/,$_;
	my $gene_ID;
	if ($line[0] eq "#GeneID") { 
		$_=~s/^#//;
		print OUT  "Chr\tStart\tEnd\tPredicted_copy_number\tType_of_alteration","\t$_\n";
	}else{
		if (defined $hash4{$line[0]}) { 
			print OUT $hash4{$line[0]},"\t$_\n";
		}
		else{
			next;
		}
	}
}
close IN;
close INN;
close INNN;
close OUT;



################################################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub USAGE {
	my $usage=<<"USAGE";

      Usage:
		Options:
		-region <file>	input file,xxx format,forced

		-gff <file>	input file,xxx format,forced

		-anno <file>	Integrated_Function.annotation.xls

		-o <file>	output file,optional
		
		-type  <str>    gene/mRNA   [gene]

		-h		help

USAGE
	print $usage;
	exit;
}
