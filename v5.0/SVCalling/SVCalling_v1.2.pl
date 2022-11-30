#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);

my %opts;
GetOptions(\%opts,"id=s","c=s","od=s","h" );

if(!defined($opts{id}) || !defined($opts{c}) || !defined($opts{od}) || defined($opts{h}))
{
	print <<"	Usage End.";
	
	Usage:

		-id          bam files indir                           must be given

		-od          outdir                                    must be given
	
		-c           config                                    must be given

		-h           Help document

	Usage End.

	exit;
}
###############################################################################################
my $config = $opts{c};
my %hash;
open CFG,$config or die $!;
while(<CFG>){
        chomp;
        next if /^$/ || /^\#/;
        my ($key,$value) = split /\s+/;
        $hash{$key} = $value;
}
close CFG;

my $chr_id = $hash{chr_id};
my $gfffile = $hash{Gff2};
my $outprefix = $hash{Species};
my $maxproc = $hash{maxcpu};
my $queue = $hash{queue};
my $type = $hash{type};
my $breakdancer_dir = "$Bin/breakdancer/" ;
my $ref = $hash{Ref2};
my $score = 20;

my $indir = $opts{id} ;
$indir = abs_path($indir);
my $outdir = $opts{od} ;
mkdir($outdir,0755) unless -d $outdir;
$outdir = abs_path($outdir);

###############################################################################################
my @bamfiles = ();
&get_bamfiles(\@bamfiles);

# use breakdancer calling SV
&sv_calling(\@bamfiles) ;

# statistic result
&sv_statistic();

# sv annotation
&sv_anno();


###############################################################################################
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

## get bam files
sub get_bamfiles()
{
	my ($abamfiles) = @_ ;
	@{$abamfiles} = glob("$indir/*.bam");
	if (scalar (@{$abamfiles}) == 0){
		print "No bam file found, please check\n" ;
	}

	return ;
}

sub get_sample_bam_files()
{
	my ($abamfiles, $ahsamples_bam) = @_ ;
	for my $bamfile (@{$abamfiles}){
		my $basename = basename($bamfile);
		my $sample_id = $basename ;
		if ($basename =~ /\.(R\d+)\.sort.bam$/ || $basename =~ /\.(R\d+)L\d+\.sort.bam$/){
			$sample_id = $1 ;
		}
		push @{$ahsamples_bam->{$sample_id}}, $bamfile ;
	}
	return ;
}

sub sv_calling()
{
	my ($abamfiles) = @_ ;
	my $shfile = "$outdir/$outprefix.SV.sh" ;
	my $filter_shfile = "$outdir/$outprefix.filterSV.sh" ;
	open (SH, ">$shfile") || die "Can't creat $shfile, $!\n" ;
	open (FS, ">$filter_shfile") || die "Can't creat $filter_shfile, $!\n" ;
	my %hsamples_bam = ();
	&get_sample_bam_files($abamfiles, \%hsamples_bam);
	for my $sample (sort keys %hsamples_bam){
		my $basename = "$outprefix.$sample" ;
		my $out_config = "$outdir/$basename.config" ;
		my $out_SV = "$outdir/$basename.max.tmp" ;
		print SH "perl $breakdancer_dir/bam2cfg.pl" ;
		for my $bamfile (@{$hsamples_bam{$sample}}){
			print SH " $bamfile" ;
		}
		print SH " > $out_config && " ;
		print SH "perl $breakdancer_dir/BreakDancerMax.pl $out_config > $out_SV\n";
		my $filter_SV = "$outdir/$basename.max" ;
		print FS "perl $Bin/sv_filter.pl -i $out_SV -o $filter_SV -s $score -c $chr_id\n" ;
	}
	close(SH);
	close(FS);
	&qsub($shfile,$maxproc);
	&qsub($filter_shfile,$maxproc);
	return ;
}

sub sv_statistic()
{
	my $sv_dir = $outdir ;
	my $outfile = "$outdir/$outprefix.sv.stat.xls" ;
	my $cmd = "perl $Bin/sv_stat.pl -id $sv_dir -o $outfile" ;
	&run_or_die($cmd);

	return ;
}

sub sv_anno()
{
	my $sv_dir = $outdir ;
	my $outfile = "$outdir/$outprefix.sv.anno" ;
	my $cmd = "perl $Bin/sv_anno_v1.3.pl -id $sv_dir -g $gfffile -o $outfile -annotype $type" ;
	&run_or_die($cmd);

	return ;
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
	my $start_time = &show_log($cmd);
	my $flag = system($cmd) ;
	if ($flag != 0){
		my $end_time = &show_log("Error: command fail: $cmd");
		&show_log("Elaseped time: ".($end_time - $start_time)."s\n");
		exit(1);
	}
	my $end_time = &show_log("done.");
	&show_log("Elaseped time: ".($end_time - $start_time)."s\n");
	return ;
}

sub qsub()
{
	my ($shfile, $maxproc) = @_ ;
	my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc $maxproc --queue $queue --reqsub $shfile --independent" ;
	&run_or_die($cmd);
	return ;
}

