#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($Config,$fOut,$Q);
GetOptions(
	"help|?" =>\&USAGE,
	"o:s"=>\$fOut,
	"dataconfig:s"=>\$Config,
	"Q:s"=>\$Q,
) or &USAGE;
&USAGE unless ($Config and $fOut);
$Q||=33;
my %seq_hash;
my $cfg=abs_path($Config);
my $outdir=abs_path($fOut);
unless (-d "$outdir/work_sh")
{
	`mkdir -p $outdir/work_sh`;
}
open (IN,$cfg) or die $!;
my $sample;
while (<IN>){
	chomp;
	next if /^$/ || /^\#/;
	my @data=split /\s+/,$_;
	if ($data[0]=~/Sample/i){
		$sample=$data[1];
	}
	if ($data[0]=~/^fq1/){
		$seq_hash{$sample}{1}=$data[1];
	}
	if ($data[0]=~/^fq2/){
		$seq_hash{$sample}{2}=$data[1];
	}
}
open (FASTP,">$outdir/work_sh/fastp.sh") or die $!;
foreach my $sam (sort keys %seq_hash){	
	chomp $sam;
	unless (-d "$outdir/$sam"){
		`mkdir $outdir/$sam`;
	}
	print  FASTP "cd $outdir/$sam && $qc_stat -a $seq_hash{$sam}{1} -b $seq_hash{$sam}{2} -f $sam -q 45 -Q $Q\n";
}
close(FASTP);
&qsub("$outdir/work_sh/gc_dis.sh");
my @Sample=glob("$outdir/R*/");
my @Paraf;
foreach my $sample (@Sample) 
{
	my @STAT=glob("$sample/*.stat");
	foreach my $stat (@STAT) 
	{
		push @Paraf,$stat unless $stat=~/.+cycleQ\.stat/;
	}
}
open (OUT,">$outdir/AllSample_GC_Q.stat");
print OUT "SampleID\tReadSum\tBaseSum\tGC(%)\tN(%)\tQ20(%)\tCycleQ20(%)\tQ30(%)\n";
foreach my $paraf (@Paraf) 
{
	my $now;
	my $file=basename($paraf);
	$file=~s/\.stat$//;
	open(IN,"$paraf")||die"can't open $paraf\n";
	<IN>;
	while (<IN>)
	{
		chomp;
		my @A=split/\s+/,$_,2;
		$now="$file\t$A[1]\n";
	}
	close(IN);
	print OUT $now;
}
close(OUT);

#######################################################################################
sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &sub_format_datetime(localtime($time));
	print "$Time:\t$txt\n" ;
	return ($time) ;
}
#############################################################
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
#######################################################################################

sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Usage:
  -dataconfig <file>  input config file
	
  -Q <file>  quality shift for different plotform [optional] [default 33]
  
  -o <file>  outdir 
 
  -h         Help

USAGE
	print $usage;
	exit;
}
