#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
#===============================================================

my %opts;
GetOptions(\%opts,"i=s","o=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{i}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";

	Usage:

		-i           snp or indel vcf file                <infile>     must be given
		-o           selected snp result file             <outfile>    must be given

		-h           Help document

	Usage End.

	exit;
}

# get parameters
my $infile  = $opts{i} ;
my $outfile = $opts{o} ;

# reading vcf file and filters
&reading_vcf_file_and_filter($infile, $outfile);



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

sub reading_vcf_file_and_filter()
{
	my ($infile, $outfile) = @_ ;
	open (IN, $infile) || die "Can't open $infile, $!\n" ;
	open (OUT, ">$outfile") || die "Can't creat $outfile, $!\n" ;
	while(<IN>){
		chomp ;
		next if (m/^\s*$/);
		if (m/^\#/){
			print OUT $_,"\n" ;
			next ;
		}
		my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split ;
		my $flag = &judge_whether_is_different(\@samples);
		if ($flag == 1){
			print OUT $_,"\n" ;
		}
	}
	close(IN);
	close(OUT);

	return ;
}

#my $flag = &judge_whether_is_different(\@samples);
sub judge_whether_is_different()
{
	my ($asamples) = @_ ;
	my $flag = 0 ;
	my %hstat = ();
	for my $info (@{$asamples}){
		my ($type) = split /\:/, $info ;
		next if ($type eq './.');
		$hstat{$type} ++ ;
	}
	$flag = 1 if (scalar(keys %hstat) > 1);

	return($flag) ;
}


