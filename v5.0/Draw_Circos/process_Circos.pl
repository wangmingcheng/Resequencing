#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($cfg);
GetOptions(
	"help|?" =>\&USAGE,
	"cfg:s"=>\$cfg,
) or &USAGE;
&USAGE unless ($cfg);

my %cfg_hash;
my %proinfo_hash;
my %info_hash;

open (IN,$cfg) or die $!;
while (<IN>) {
	chomp;
	next if /^$/ || /^\#/;
	my @data=split /\s+/,$_;
	$cfg_hash{$data[0]}=$data[1];
}
close IN;
my $assemble_level = 'chr';
my $chr_id = $cfg_hash{chr_id};
my $key = $cfg_hash{Project};
my $queue = $cfg_hash{queue};
my $maxproc = $cfg_hash{maxcpu};
my $indir = $cfg_hash{analysis_dir};
mkdir("$indir/Analysis/Circos",0755) unless -d "$indir/Analysis/Circos";
mkdir("$indir/Analysis/Circos/work_sh",0755) unless -d "$indir/Analysis/Circos/work_sh";
my %chr_id_hash;
open(IN, $chr_id);
while(<IN>){
	chomp;
	next if(/^$/);
	next if(/^\#/);
	my ($chr, $others) = split(/\s+/, $_, 2);
	$chr_id_hash{$chr} = $others;
}
close(IN);

my $ref = "$indir/Analysis/Ref/sequences.fa";
my $chr_id_sort_R = "$indir/Analysis/Circos/chr_id_sort.R";
my $chr_id_sorted = "$indir/Analysis/Circos/chr_id_sorted.txt";
open(OUT, ">$chr_id_sort_R") or die $!;
print OUT <<"	END";
id_tab <- read.table("$chr_id", header = T, stringsAsFactors = F, check.names = F, comment.char = \"\")
library(gtools)
library(ggplot2)
library(scales)
multi.mixedorder <- function(..., na.last = TRUE, decreasing = FALSE){
  do.call(order, c(
    lapply(list(...), function(l){
      if(is.character(l)){
        factor(l, levels=mixedsort(unique(l)))
      } else {
        l
      }
    }),
    list(na.last = na.last, decreasing = decreasing)
  ))
}
id_tab_order <- id_tab[multi.mixedorder(id_tab[,1]),]
names(id_tab_order)[1] <- "#Chrom"
write.table(id_tab_order, "$chr_id_sorted", sep = \"\\t\", quote = F, row.names = F, col.names = T)
	END
close(OUT);
my $cmd = "R --slave < $chr_id_sort_R";
&run_or_die($cmd);

open (IN,$chr_id_sorted) or die $!;
my $chr_num=0;
my %chr_hash;
while (<IN>){
	$chr_num++;
	chomp;
	next if(/^$/);
	next if(/^\#/);
	my @chr=split/\t/,$_;
	$chr_hash{chrnum}.=$chr[0].",";
	if ($chr_num>30){
		if($assemble_level ne "chr"){
			last;
		}
	}
}

my $Circle_type="scatter".","."line".","."histogram";
my @sample=glob "$indir/Analysis/Annotation/snp_anno/R*";
open (CIRCOS,">$indir/Analysis/Circos/work_sh/draw_circos.sh") or die $!;
my @svfiles = glob("$indir/Analysis/SV/$key*.max");
my@cnvFiles;
if (-d "$indir/Analysis/CNV/" ){
	@cnvFiles=glob("$indir/Analysis/CNV/*/$key*CNVs");
}

for (my $n=0;$n<@sample ;$n++){
	my $sample=basename($sample[$n]);
	next if ($sample !~ /^R\d+/);
	unless (-d "$indir/Analysis/Circos/$sample"){
		`mkdir $indir/Analysis/Circos/$sample`;
	}
	open (OUT,">$indir/Analysis/Circos/$sample/circos.config") or die $!;
	my $snpvcf=(glob("$sample[$n]/*.vcf"))[0];
	my $indelvcf=(glob("$sample[$n]/../../indel_anno/$sample/*.vcf"))[0];
	my( $svmax ,$cnv)=("","");
	for (my $i=0; $i<@svfiles; $i++){
		if ($svfiles[$i] =~ /$sample/){
			$svmax = $svfiles[$i] ;
		}
	}
		
	if (-d "$indir/Analysis/CNV/" ){
		for (my $i=0; $i<@cnvFiles; $i++){
			if ($cnvFiles[$i] =~ /$sample/){
				$cnv = $cnvFiles[$i] ;
			}
		}
	}
	print OUT "Ref\t$ref\n";
	print OUT "Chr\t$chr_hash{chrnum}\n";
	print OUT "Circle_type\t$Circle_type\n";
	print OUT "SNPvcf\t$snpvcf\n";
	print OUT "indelvcf\t$indelvcf\n";
	if (-f "$svmax" ){
		print OUT "SVinfo\t$svmax\n" ;
	}
	if (-d "$indir/Analysis/CNV/" ){
		print OUT "CNVinfo\t$cnv\n" ;
	}
	close(OUT);
	print CIRCOS "perl $Bin/getcircos.config.pl -cfg $indir/Analysis/Circos/$sample/circos.config -key $sample -o $indir/Analysis/Circos/$sample/ -detail $cfg\n";
}
close(CIRCOS);
&qsub("$indir/Analysis/Circos/work_sh/draw_circos.sh");


close (IN) ;
close (OUT) ;

#######################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &sub_format_datetime(localtime($time));
	print "$Time:\t$txt\n" ;
	return ($time) ;
}

sub run_or_die()
{
	my ($cmd) = @_ ;
	&show_log($cmd);
	my $flag = system($cmd) ;
	if ($flag != 0){
		&show_log("Error: command fail: $cmd");
		exit(1);
	}
	&show_log("done.");
	return ;
}
######################################################################
sub qsub()
{
	my ($shfile) = @_ ;
	my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc $maxproc --queue $queue --resource vf=15G --reqsub $shfile --independent" ;
	&run_or_die($cmd);
}
##########################################################################
sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";

Usage:

  -cfg <file>  input file,fasta format,forced 
  -h         Help

USAGE
	print $usage;
	exit;
}
