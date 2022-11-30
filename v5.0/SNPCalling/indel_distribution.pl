#!/usr/bin/perl -w
# 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

####################

my %opts;
GetOptions(\%opts,"i=s","o=s","m=s","x=s","h" );

if(!defined($opts{i}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";

	Usage:

		-i           dis infile                                must be given

		-o           out png file                              must be given

		-m           minimum indel length, default -10         optional

		-x           maximum indel length, default 10          optional

		-h           Help document

	Usage End.

	exit;
}

# get parameters
my $infile  = $opts{i} ;
my $outfile = $opts{o} ;
my $min_length = defined $opts{m} ? $opts{m} : -10 ;
my $max_length = defined $opts{x} ? $opts{x} : 10 ;

my $disfile = &refine_file($infile);

my $cmd = "Rscript $Bin/bar_split.r --infile $disfile --outfile $outfile --main.col 1 --sub.col 2 --value.col 3 --main.lab \"Length of Indel\" --sub.lab \"Region\" --value.lab \"Number of Indel\" --title.lab \"Indel length distribution\" --legend.xpos 0.85 --legend.ypos 0.85" ;
&run_or_die($cmd);


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
	&show_log($cmd);
	my $flag = system($cmd) ;
	if ($flag != 0){
		&show_log("Error: command fail: $cmd");
		exit(1);
	}
	&show_log("done.");
	return ;
}


#my $disfile = &refine_file($infile);
sub refine_file()
{
	my ($infile) = @_ ;
	(my $disfile = $infile) =~ s/$/.cut/ ;
	open (IN, $infile) || die "Can't open $infile, $!\n" ;
	open (OUT, ">$disfile") || die "Can't creat $disfile, $!\n" ;
	while (<IN>){
		chomp ;
		next if (m/^\s*$/) ;
		if (m/^-?\d+/){
			my $length = (split)[0] ;
			if ($length >= $min_length && $length <= $max_length){
				print OUT $_,"\n" ;
			}
		}
		else{
			print OUT $_,"\n" ;
		}
	}
	close(IN);
	close(OUT);

	return ($disfile);
}


