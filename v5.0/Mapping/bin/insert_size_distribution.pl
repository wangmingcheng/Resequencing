#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

my %opts;
GetOptions(\%opts,"i=s","o=s","m=s","x=s","c=s","h" );

if(!defined($opts{i}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";
	
    Usage:

		-i           alignment file (bam format)                    <infile>     must be given
		-o           result file prefix                             <outfile>    must be given
		-c	     config
		-m           minimum insert size for draw, default 100      [int]        optional
		-x           maximum insert size for draw, default 500      [int]        optional
		-h           Help document

	Usage End.

	exit;
}
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
my $java = $hash{java};
my $picard_dir = $hash{picard_dir};

# get parameters
my $bamfile = $opts{i} ;
my $outfile = $opts{o} ;
my $min_insert = defined $opts{m} ? $opts{m} : 300 ;
my $max_insert = defined $opts{x} ? $opts{x} : 800 ;

# get insert size dis
my $disfile = &stat_insert_size_dis($bamfile, $outfile);

# draw insert size dis
&draw_insert_size_dis($disfile, $outfile);


###############Subs
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &sub_format_datetime(localtime($time));
	print "$Time:\t$txt\n" ;
	return ($time) ;
}

#&run_or_die($cmd);
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


sub stat_insert_size_dis()
{
	my ($bamfile, $outfile) = @_ ;
	# stat insert size
	my $statfile = "$outfile.insert" ;
	my $histogram_file = "$outfile.insert.pdf" ;
	my $cmd = "$java -jar $picard_dir/CollectInsertSizeMetrics.jar VALIDATION_STRINGENCY=LENIENT HISTOGRAM_FILE=$histogram_file" ;
	$cmd .= " INPUT=$bamfile OUTPUT=$statfile STOP_AFTER=10000000 DEVIATIONS=20" ;
	&run_or_die($cmd) ;

	# get dis file
	my $disfile = "$outfile.insert.dis" ;
	$cmd = "perl $Bin/get_insert_dis_from_picard_result.pl -i $statfile -o $disfile" ;
	&run_or_die($cmd) ;

	return($disfile);
}

sub draw_insert_size_dis()
{
	my ($disfile, $outfile) = @_ ;
	my $cmd = "perl $Bin/draw_insert_size.pl -i $disfile -o $outfile -m $min_insert -x $max_insert" ;

	&run_or_die($cmd);

	return ;
}


