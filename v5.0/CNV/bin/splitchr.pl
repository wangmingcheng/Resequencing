#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;
#######################################################################################

# ==============================================================
# Get Options
# ==============================================================
my ($fIn,$fKey,$dOut,$log, $chr_id);

GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"od:s"=>\$dOut,
				"chr:s"=>\$chr_id,
				) or &USAGE;
&USAGE unless ($fIn);

#===============================================================
# Default optional value 
#===============================================================
$dOut||="./ref";
mkdir $dOut unless (-d $dOut) ;
mkdir $dOut."/chr";
#===============================================================
# Global value
#===============================================================
open(IN, $chr_id);
my %chr_hash;
while(<IN>){
	chomp;
	next if (/^$/);
	next if (/^\#/);
	my ($chr, $other) = split(/\s+/,$_, 2);
	$chr_hash{$chr} = $other;
}
close(IN);
print Dumper %chr_hash, "\n";
open (IN,$fIn) or die $!;
undef $/;
$/=">";
open (STAT,">$dOut/res.stat") or die $!;
while (<IN>) {
	next unless (my ($chr,$seq) = /(.*?)\n(.*)/s);
	my @line=split(/\s+/,$chr);
	next if(!exists $chr_hash{$line[0]});
	$seq=~s/>//g;
	open (OUT,">$dOut/chr/$line[0].fa") or die $!;
	print OUT ">$line[0]\n$seq\n";
	close (OUT);
	$seq =~ s/[\d\s>]//g;
	print STAT "$line[0]\t".length($seq)."\n";
}
$/="\n";
close (IN);
close(STAT);
#===============================================================
# Process
#===============================================================



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ==============================================================
# sub function
# ==============================================================
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Shi Tongwei <shitw\@biomarker.com.cn> 
Discription: split the reference fasta into chromosomes

Usage:
  Options:
  -i	<file>	reference genome, fasta format, forced
  -od	<str>	Directory where output file produced,optional,default [./ref]
  -chr  <file>  chromosome id file
  -h		Help

USAGE
	print $usage;
	exit;
}
