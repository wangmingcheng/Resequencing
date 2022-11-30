#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);

my ($fIn,$ploidy,$win,$ref,$reflen,$GCrange,$conffile,$outdir,$chr_id,$samtools);

GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"chr:s"=>\$chr_id,
	"ref:s"=>\$ref,
	"reflen:s"=>\$reflen,
	"ploidy:s"=>\$ploidy,
	"window:s"=>\$win,
	"GCrange:s"=>\$GCrange,
	"od:s"=>\$outdir,
	"samtools:s"=>\$samtools,
) or &USAGE;
&USAGE unless ($fIn and $ref and $ploidy);

####################################################################################################
my $_freec_="$Bin/bin/freec";
my $_splitchr_="$Bin/bin/splitchr.pl";
my $_svg2xxx_="$Bin/svg2xxx";

$conffile||=$outdir."/config.txt";
$outdir||="./FREEC_res";
mkdir $outdir if(!-e $outdir);
$outdir=abs_path($outdir);
$conffile=abs_path($conffile);
$fIn=abs_path($fIn);
$win||=50000;
$GCrange||="0.35-0.55";
my @gc=split(/-/,$GCrange);

system "perl $_splitchr_ -i $ref -od $outdir -chr $chr_id";
my $filename="$outdir/".basename($fIn); `ln -s $fIn $filename`;

open (OUT,">$conffile") or die $!;
my $line=<<"MYFORMAT";
[general]
chrLenFile = $outdir/res.stat
ploidy = $ploidy
breakPointThreshold = 0.8
window = $win
chrFiles = $outdir/chr
minExpectedGC = $gc[0]
maxExpectedGC = $gc[1]
numberOfProcesses = 8
outputDir = $outdir
samtools = $samtools
[sample]
mateFile = $filename
inputFormat = BAM
mateOrientation = FR
MYFORMAT

print OUT $line;
close (OUT) ;

system "$_freec_ -conf $conffile";
my $lenname="$outdir/".basename($ref).".len";
open(IN, $chr_id);
my %chr_hash;
while(<IN>){
	chomp;
	next if (/^$/);
	next if (/\#/);
	my ($chr_use, $others) = split(/\s+/, $_, 2);
	$chr_hash{$chr_use} = $others;
}
close(IN);
open(IN, $reflen);
open(OUT, ">$lenname") || die $!;
while(<IN>){
	chomp;
	next if(/^$/);
	if(/^\#/){
		print OUT "$_\n";
	}else{
		my $line_temp = $_;
		my ($chr_use, $others) = split(/\t/, $_, 2);
		if(exists $chr_hash{$chr_use}){
			print OUT "$line_temp\n";
		}
	}
}
close(IN);
close(OUT);

sub USAGE {#
	my $usage=<<"USAGE";

Usage:

  -i		    <file>	*.sam/*.bam, required
  -ref		    <str>	reference genome file ,fasta format,required
  -chr              <file>      chromosome ID file
  -reflen           <file>      chromosome/scaffold length file
  -ploidy	    <num>	ploidy number,required
  -window	    <num>	explicit window size , optional , default [50000]
  -GCrange	    <str>	exptected value of the GC-content, optional, default [0.35-0.55]
  -od	            <str>	Directory where FREEC output file produced,optional,default [./FREEC_res]

  -h		            Help
USAGE
	print $usage;
	exit;
}
