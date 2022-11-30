#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

my %opts;
GetOptions(\%opts,"i=s","o=s","m=s","x=s","h" );

if(!defined($opts{i}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";

	Usage:

		-i           insert size distribution file                    <infile>     must be given
		-o           result file prefix                               <outfile>    must be given
		-m           minimum insert size for draw, default 100        [int]        optional
		-x           maximum insert size for draw, default 500        [int]        optional
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
my $min_insert = defined $opts{m} ? $opts{m} : 100 ;
my $max_insert = defined $opts{x} ? $opts{x} : 500 ;

# reading dis file and cut
my $cut_file = &generate_dis_file_for_draw($infile, $outfile, $min_insert, $max_insert);

# draw insert size dis png
&draw_insert_distribution($cut_file, $outfile);

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


# &show_log("txt")
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


sub generate_dis_file_for_draw()
{
	my ($infile, $outfile, $min_insert, $max_insert) = @_ ;
	my $basename = basename($infile);
	my $dirname = dirname($outfile);
	my $cut_file = "$dirname/$basename.cut" ;
	open (IN, $infile) || die "Can't open $infile, $!\n" ;
	open (OUT, ">$cut_file") || die "Can't creat $cut_file, $!\n" ;
	while (<IN>){
		chomp ;
		next if (m/^\s*$/);
		if (m/^\#/){
			print OUT $_, "\n" ;
			next ;
		}
		my ($insert, $num) = split ;
		if ($insert >= $min_insert && $insert <= $max_insert){
			print OUT $_, "\n" ;
		}
	}
	close(IN);
	close(OUT);

	return($cut_file);
}

#&draw_insert_distribution($cut_file, $outfile);
sub draw_insert_distribution()
{
	my ($cut_file, $outfile) = @_ ;
	my $cmd = "Rscript $Bin/singleArea.r --infile $cut_file --outfile $outfile.cut.png --x.col 1 --y.col 2 --x.lab \"Insert size\"" ;
	$cmd .= "  --y.lab \"Reads number\" --skip 1 --no.grid --color \"#FF66EE\" --title.lab \"Insert size distribution\"" ;

	&run_or_die($cmd);

	return ;
}


