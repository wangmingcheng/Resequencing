#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

my %opts;
GetOptions(\%opts,"i=s","o=s","h" );
if(!defined($opts{i}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";

	Usage:

		-i           picard insert size result file               <infile>     must be given
		-o           insert size distribution result file         <outfile>    must be given
		-h           Help document

	Usage End.

	exit;
}

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
################
# get parameters
my $infile  = $opts{i} ;
my $outfile = $opts{o} ;

# reading infile and output
&reading_infile_and_output($infile, $outfile);

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";

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

sub reading_infile_and_output()
{
	my ($infile, $outfile) = @_ ;
	open (IN, $infile) || die "Can't open $infile, $!\n" ;
	open (OUT, ">$outfile") || die "Can't creat $outfile, $!\n" ;
	my $flag = 0 ;
	while (<IN>){
		chomp ;
		next if (m/^\#/ || m/^\s*$/) ;
		if (m/^insert_size/){
			$flag = 1 ;
			print OUT "#$_\n" ;
			next ;
		}
		if ($flag == 1){
			print OUT $_, "\n" ;
		}
	}

	close(IN);
	close(OUT);

	return ;
}


