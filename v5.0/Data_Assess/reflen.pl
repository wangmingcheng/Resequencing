#!/usr/bin/perl -w
# 
# Copyright (c) BMK 2009
# Writer:         Yangsh <ysh123.love@163.com>
# Program Date:   2010.
# Modifier:       Yangsh <ysh123.love@163.com>
# Last Modified:  2010.
my $ver="1.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

######################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作

my %opts;
GetOptions(\%opts,"fa=s","o=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{fa}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		
		Version: $ver

	Usage:

		-fa    infile     must be given

		-o    outfile    must be given

		-h    Help document

	Usage End.

	exit;
}

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
################
my $programe_dir=basename($0);
my $in=$opts{fa};
my $out=$opts{o};
$/=">";
open (OUT,">$out");
open (IN,"$in")||die"Can't open $in\n";
<IN>;
print OUT "#Chrom\tAGCTN\tAGCT\tN\n";
my ($ATCGN,$ATCG,$N,$num)=(0,0,0,0);
while (<IN>)
{
	chomp;
	my ($head,$seq)=(split/\n/,$_,2);
	$head=(split/\s+/,$head)[0];
	$num++;
	$seq=~s/\s+//g;
	$seq=uc($seq);
	my $len=length($seq);
	my $seq1=$seq;$seq1=~s/N//g;
	my $len1=length($seq1);
	my $len2=$len-$len1;
	print OUT "$head\t$len\t$len1\t$len2\n";
	#print  "$head\t$len\t$len1\t$len2\n";
	$ATCGN+=$len;
	$ATCG+=$len1;
	$N+=$len2;
}
close(IN);
print OUT "#Total_Len\t$ATCGN\t$ATCG\t$N\n";

print  "\n-------------Total info-------------\n";
print  "Number :$num\n";
print  "ATCGN  :$ATCGN\n";
print  "ATCG   :$ATCG\n";
print  "N      :$N\n";
$/="\n";
close(OUT);
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
