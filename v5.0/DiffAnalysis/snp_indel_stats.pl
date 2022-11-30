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
my %project_hash;
open(IN,$cfg) or die $!;
while (<IN>){
	chomp;
	next if ($_=~/\#/);
	next if ($_=~/^\s*$/);
	my @data=split/\s+/,$_;
	$cfg_hash{$data[0]}=$data[1];
}
my $type = $cfg_hash{type};
my $queue = $cfg_hash{queue};
my $maxcpu = $cfg_hash{maxcpu};

my $indir="$cfg_hash{analysis_dir}/";
if (-d "$indir/Analysis/Annotation/snp_anno/"){
	unless (-d "$indir/Analysis/Annotation/snp_anno/sh_dir"){
		`mkdir $indir/Analysis/Annotation/snp_anno/sh_dir`
	}
	my @sam=glob "$indir/Analysis/Annotation/snp_anno/R*";
	open (SH,">$indir/Analysis/Annotation/snp_anno/sh_dir/snp_indel.stat.sh");
	for (my $i=0;$i<@sam;$i++){
		my @vcf=glob "$sam[$i]/*.vcf";
		for (my $j=0;$j<@vcf ;$j++){
			my $bvcf=basename($vcf[$j]);
			(my $list=$bvcf)=~s/\.vcf$/\.list/;
			print SH "perl $Bin/stat_v1.3.pl $vcf[$j] $sam[$i]/$list gene\n" if $type eq "gene";
			print SH "perl $Bin/stat_v1.3.pl $vcf[$j] $sam[$i]/$list mRNA\n" if $type eq "mRNA";
		}
	}
	close(SH);
	&qsub("$indir/Analysis/Annotation/snp_anno/sh_dir/snp_indel.stat.sh");
}
if (-d "$indir/Analysis/Annotation/indel_anno/"){
	unless (-d "$indir/Analysis/Annotation/indel_anno/sh_dir"){
		`mkdir $indir/Analysis/Annotation/indel_anno/sh_dir`
	}
	my @sam=glob "$indir/Analysis/Annotation/indel_anno/R*";
	open (SH,">$indir/Analysis/Annotation/indel_anno/sh_dir/snp_indel.stat.sh");
	for (my $i=0;$i<@sam;$i++){
		my @vcf=glob "$sam[$i]/*.vcf";
		for (my $j=0;$j<@vcf ;$j++){
			my $bvcf=basename($vcf[$j]);
			(my $list=$bvcf)=~s/\.vcf$/\.list/;
			print SH "perl $Bin/stat_v1.3.pl $vcf[$j] $sam[$i]/$list gene\n" if $type eq "gene";
			print SH "perl $Bin/stat_v1.3.pl $vcf[$j] $sam[$i]/$list mRNA\n" if $type eq "mRNA";
		}
	}
	close(SH);
	&qsub("$indir/Analysis/Annotation/indel_anno/sh_dir/snp_indel.stat.sh");
}


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

sub qsub()
{
	my ($shfile) = @_ ;
	my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc $maxcpu --queue $queue --resource vf=15G --reqsub $shfile --independent" ;
	&run_or_die($cmd);
}


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
