#!/usr/local/bin/perl -w
# Copyright (c) BMK 2014/3/4
# Writer:         lium <lium@biomarker.com.cn>
# Program Date:   2014/3/4
# Last Modified:  2014/3/4
my $ver="1.0";

use strict;
use Cwd;
use Getopt::Long;
my $BEGIN=time();


# get opts
my %opts;
GetOptions(\%opts,"i=s","key=s","h" );
if(! defined($opts{i}) ||! defined($opts{key}) || defined($opts{h})){
	&USAGE;
}


###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";


# open file
my $type_file = "$opts{key}.type.stat";
my $region_file = "$opts{key}.region.stat";
open (IN,"$opts{i}") || die "Can't open file $opts{i}\n";
open (OUT,">$type_file") || die "open or create file $type_file is failed\n";
open (OUT2,">$region_file") || die "open or create file $region_file is failed\n";



# ignore no-used lines
my $ignore_ok = 0;
while(<IN>){
	chomp;
	next if (/^$/);		# ignore the blank lines
	# check individual info
	if(/Number of effects by type and region/){
		$ignore_ok = 1;
		last;
	}
}


# check ignore result
if($ignore_ok==0){
	print "error: Number of effects by type and region not found\n";
	exit;
}


# abstract type stat
my $thead_ok = 0;
my $th_count = 0;
while(<IN>){
	chomp;
	next if (/^$/);		# ignore the blank lines
	# check <thead>
	if(/<thead>/){
		$thead_ok++;
		next;
	}
	# check over condition
	last if(/<\/table>/);
	# check thread
	next if($thead_ok==0);

	# get header
	if( ($_=~/<th/) || ($_=~/<td/)){
		# get content
		my $con = $_;
		$con =~ s/<b>//g;
		$con =~ s/<\/b>//g;
		($con) = $con=~/>\s+([\S ]+)\s+</;
		# print
		#print "$th_count: $_\n";
		if( ($th_count%4==0) || ($th_count%4==2) ){
			print OUT "$con\t";
		}
		if( $th_count%4==3 ) {
			print OUT "$con\n";
		}
		# add count
		$th_count ++;
	}
}

close OUT;


# check 
if($thead_ok!=1){
	print "stat stable 2 not found\n";
	exit;
}


# abstract region stat
$th_count = 0;
$thead_ok = 0;
while(<IN>){
	chomp;
	next if (/^$/);		# ignore the blank lines
	# check <thead>
	if(/<thead>/){
		$thead_ok++;
		next;
	}
	# check over condition
	last if(/<\/table>/);
	# check thread
	next if($thead_ok==0);

	# get info 
	if( ($_=~/<th/) || ($_=~/<td/)){
		# get content
		my $con = $_;
		$con =~ s/<b>//g;
		$con =~ s/<\/b>//g;
		($con) = $con=~/>\s+([\S ]+)\s+</;
		# print
		if( ($th_count%4==0) || ($th_count%4==2) ){
			print OUT2 "$con\t";
		}
		if( $th_count%4==3 ) {
			print OUT2 "$con\n";
		}
		# add count
		$th_count ++;
	}
}


close OUT2;
close IN;


###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";
&Runtime($BEGIN);


sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub Runtime{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total elapsed time: ${t}s\n\n";
}




sub USAGE{
print <<"Usage End.";
Description: abstract stat inforamtion from snpEff summary file (html format)
	Version: $ver
Usage: perl abs_stat_snpEff_summary.pl -i snpEff_summary.html -key out
	-i	SNP rawdata file(vcf format)		must be given

	-key	the prefix for output			must be given

	-h	Help document
Usage End.
exit;
}


