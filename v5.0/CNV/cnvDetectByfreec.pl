#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use Cwd qw(abs_path getcwd);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use threads;

my (%opts,$od,$id,$GCrange);
GetOptions(\%opts,"c=s","od=s","id=s","GCrange=s","h");

if(!defined($opts{c}) || !defined($opts{od}) || !defined($opts{id}) ||!defined($opts{GCrange})  ||defined($opts{h}))
{
	print <<"	Usage End.";

		-od          outdir                            <dir>        must be given
		-id          indir of bam file                 <dir>        must be given
		-GCrange     GC range                                       must be given
		-c           config file 

	Usage End.

	exit(1);
}
####################################################################################################
my %hash;
my $config = $opts{c};
$GCrange = $opts{GCrange};
open CFG,$config or die $!;
while(<CFG>){
        chomp;
        next if /^$/ || /^\#/;
        my ($key,$value) = split /\s+/;
        $hash{$key} = $value;
}
close CFG;
my $key = $hash{Species};
my $chr_num = $hash{ChrNum};
my $queue = $hash{queue};
my $maxcpu = $hash{maxcpu};
my $type = $hash{type};
my $chr_id = $hash{chr_id};
my $ploidy = $hash{Ploidy};
my $ref = $hash{Ref1};
my $gff = $hash{Gff1};
my $samtools = $hash{samtools};
mkdir $opts{od} unless(-d $opts{od});
$od=abs_path($opts{od});
$id=abs_path($opts{id});

my $wk_dir="$od/work_sh";
mkdir($wk_dir,0755) unless -d $wk_dir;

my @t=();
my $process = 5;
my $child_num = 0;

&call_CNV();

####################################################################################################
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub call_CNV(){
	
	my $sh1 = "$wk_dir/step5_1_CNV.detect.sh";
	my $sh2 = "$wk_dir/step5_2_CNV.ann.sh";
	
	open (SH1, ">$sh1") || die "Can't open $sh1, $!\n" ;
	open (SH2, ">$sh2") || die "Can't open $sh2, $!\n" ;

	my @bams = glob("$id/*bam");
	for my$bamfile(@bams){
		my ($sample)=(basename($bamfile)=~/(R\d+)/);
		my$bname=basename($bamfile);
		print SH1 "perl $Bin/run_freec.pl -i $bamfile -ref $ref -ploidy $ploidy -od $od/$sample -GCrange $GCrange -chr $chr_id -reflen $ref\.len -samtools $samtools\n";
		my $ann=dirname($ref)."/Function_anno/Result/Integrated_Function.annotation.xls";
		print SH2 "perl $Bin/region_gff_anno_v1.1.pl -region $od/$sample/${bname}_CNVs -gff $gff -anno $ann -o  $od/$sample/$sample.CNV.gene.anno.xls -type $type\n";
	}
	close(SH1);
	close(SH2);
	&qsub($sh1);
	&qsub($sh2);	
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
	&show_log("$cmd  -> done.");
	return ;
}

## qsub
sub qsub()
{
	my ($shfile) = @_ ;
	my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh  --maxproc $maxcpu --queue $queue --resource vf=15G --reqsub $shfile --independent" ;
	&run_or_die($cmd);
	return ;
}



