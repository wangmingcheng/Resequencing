#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

my %opts;
GetOptions(\%opts,"id=s","o=s","h");

if(!defined($opts{id}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";

	Usage:

		-id          indir of map stat files             must be given

		-o           outfile                             must be given

		-h           Help document

	Usage End.

	exit;
}

# get parameters
my $indir = $opts{id} ;
my $outfile = $opts{o} ;

# get stat files
my @statfiles = glob("$indir/*.R*.stat");

# reading stat file 
&stat_mapping_ratio(\@statfiles, $outfile);


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


sub stat_mapping_ratio()
{
	my ($astatfiles, $outfile) = @_ ;
	my %hstat = ();
	for my $statfile (@{$astatfiles}){
		&reading_stat_file($statfile, \%hstat);
	}
	
	# output
	&output_stat_result($outfile, \%hstat);

	return ;
}

sub reading_stat_file()
{
	my ($statfile, $ahstat) = @_ ;
	my ($sample) = ((basename($statfile))=~/(R\d+L\d+|R\d+)/) ;
	open (IN, $statfile) || die "Can't open $statfile, $!\n" ;
	while (<IN>){
		chomp ;
		if (m/(\d+)\s+\+\s+\d+\s+in\s+total/){
			$ahstat->{$sample}->[0] = $1 ;
		}
		elsif (m/(\d+)\s+\+\s+\d+\s+duplicates/){
			$ahstat->{$sample}->[1] = $1 ;
		}
		elsif (m/(\d+)\s+\+\s+\d+\s+secondary/){
			$ahstat->{$sample}->[2] = $1 ;
		}
		elsif (m/(\d+)\s+\+\s+\d+\s+mapped\s+\((\d+.\d+)\%/){
			$ahstat->{$sample}->[3] = $1 ;
		}
		elsif (m/(\d+)\s+\+\s+\d+\s+paired\s+in\s+sequencing/){
			$ahstat->{$sample}->[4] = $1 ;
		}
		elsif (m/\d+\s+\+\s+\d+\s+properly\s+paired\s+\((\d+.\d+)\%/){
			$ahstat->{$sample}->[5] = $1 ;
		}
	}
	close(IN);

	return ;
}

sub output_stat_result()
{
	my ($outfile, $ahstat) = @_ ;
	open (OUT, ">$outfile") || die "Can't creat $outfile, $!\n" ;
	print OUT "BMK ID\tTotal_reads\tMapped(\%)\tProperly_mapped(\%)\n" ;
	for my $sample (sort keys %{$ahstat}){
		my $map_ratio = int(10000*($ahstat->{$sample}->[3] - $ahstat->{$sample}->[2])/$ahstat->{$sample}->[4])/100 ;
		print OUT "$sample\t", $ahstat->{$sample}->[4];
		print OUT "\t", $map_ratio, "\t", $ahstat->{$sample}->[5], "\n" ;
	}
	close(OUT);

	return ;
}


