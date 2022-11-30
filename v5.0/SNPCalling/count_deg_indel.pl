#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my ($fIn,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"i:s"=>\$fIn,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);
my %indel_hash;
my @sample_array;
open (IN,$fIn) or die $!;
open (OUT,">$fOut") or die $!;
#$/=">";
while (<IN>) 
{
	chomp;
	if ($_=~/^\#/)
	{
		my @line=split/\t/,$_;
		for (my $i=0;$i<@line ;$i++)
		{
			if ($line[$i]=~/R\d+_base$/)
			{
				my $sample=(split/\_/,$line[$i])[0];
				push(@sample_array,$sample);
			}
		}
		next;
	}
	next if ($_=~/^\s+$/);
	my @data=split /\t/,$_;
	for (my $m=5;$m<@data-1 ;$m+=3)
	{
		my $sample1_num=($m-5)/3+1;#########µÈ±È
		for (my $n=5;$n<@data-1 ;$n+=3)
		{
			my $sample2_num=($n-5)/3+1;
			if ($data[$m] ne $data[$n] && $data[$m]!~/N/&&$data[$n]!~/N/)
			{
				$indel_hash{$sample_array[$sample1_num-1]}{$sample_array[$sample2_num-1]}++;
			}
			if ($sample_array[$sample1_num-1] eq $sample_array[$sample2_num-1])
			{
				$indel_hash{$sample_array[$sample1_num-1]}{$sample_array[$sample2_num-1]}=0;
			}
		}
	}
}
print OUT " \t";
foreach my $k1 (sort keys %indel_hash)
{
	print OUT "$k1\t"
}
print OUT "\n";
foreach my $sam1 (sort keys %indel_hash)
{
	print OUT "$sam1\t";
	foreach my $sam2(sort keys %{$indel_hash{$sam1}})
	{
		print OUT "$indel_hash{$sam1}{$sam2}\t";
	}print OUT "\n";
}

close (IN) ;
close (OUT) ;

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

#######################################################################################

sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Usage:
  Options:
  -i <file>  input file,forced 	eg:XXX.filter.anno.gatk.list
  
  -o <file>  output file,forced  
 
  -h         Help

USAGE
	print $usage;
	exit;
}
