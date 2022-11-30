#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

my %opts;
GetOptions(\%opts,"id=s","o=s","h" );

if(!defined($opts{id}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";

	Usage:

		-id          indir of depth stat files             must be given

		-o           outfile    must be given

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
my $indir = $opts{id} ;
my $outfile = $opts{o} ;

# get stat files
my @statfiles = glob("$indir/*depth.stat") ;

# merge stat files
&merge_stat_files(\@statfiles, $outfile);



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
	&show_log($cmd);
	my $flag = system($cmd) ;
	if ($flag != 0){
		&show_log("Error: command fail: $cmd");
		exit(1);
	}
	&show_log("done.");
	return ;
}

sub merge_stat_files()
{
	my ($astatfiles, $outfile) = @_ ;
	my %hstat = ();
	for my $statfile (@{$astatfiles}){
		&reading_stat_file($statfile, \%hstat);
	}

	# output result
	&output_stat_result($outfile, \%hstat);

	return ;
}

#&reading_stat_file($statfile, \%hstat);
sub reading_stat_file()
{
	my ($statfile, $ahstat) = @_ ;
	open (IN, $statfile) || die "Can't open $statfile, $!\n" ;
	my ($sample) = ((basename($statfile))=~/(R\d+L\d+|R\d+)/) ;
	while (<IN>){
		chomp ;
		if (m/^ave_depth\:\s+(\d+)/){
			$ahstat->{$sample}->[0] = $1 ;
		}
		elsif (m/^cov_ratio_1X\:\s+(\d+\.?\d+)/ || m/^cov_ratio_1X\:\s+(\d+)$/){
			$ahstat->{$sample}->[1] = int($1*10000)/100 ;
		}
		elsif (m/^cov_ratio_5X\:\s+(\d+\.?\d+)/ ||m/^cov_ratio_5X\:\s+(\d+)$/){
			$ahstat->{$sample}->[2] = int($1*10000)/100 ;
		}
		elsif (m/^cov_ratio_10X\:\s+(\d+\.?\d+)/ || m/^cov_ratio_10X\:\s+(\d+)$/){
			$ahstat->{$sample}->[3] = int($1*10000)/100 ;
		}
	}
	close(IN);

	return ;
}

#&output_stat_result($outfile, \%hstat);
sub output_stat_result()
{
	my ($outfile, $ahstat) = @_ ;
	open (OUT, ">$outfile") || die "Can't creat $outfile, $!\n" ;
	print OUT "BMK ID\tAve_depth\tCov_ratio_1X(\%)\tCov_ratio_5X(\%)\tCov_ratio_10X(\%)\n" ;
	for my $sample (sort keys %{$ahstat}){
		print OUT "$sample\t", $ahstat->{$sample}->[0],"\t", $ahstat->{$sample}->[1] ;
		print OUT "\t", $ahstat->{$sample}->[2],"\t", $ahstat->{$sample}->[3],"\n" ;
	}
	close(OUT);

	return ;
}

