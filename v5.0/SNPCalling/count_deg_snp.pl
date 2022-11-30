#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"i:s"=>\$fIn,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);
my %snp_hash;
open (IN,$fIn) or die $!;
open (OUT,">$fOut") or die $!;
#$/=">";
my @sample_num;
while (<IN>) 
{
	chomp;
	if ($_=~/\#/)
	{
		my @sample=split/\t/,$_;
		for (my $i=0;$i<@sample ;$i++)
		{
			if ($sample[$i]=~/^R\d+/)
			{
				push(@sample_num,$sample[$i]);
			}
		}
		next;
	}
	next if ($_=~/^\s*$/);
	my @data=split /\t/,$_;
	for (my $i=3;$i<@data;$i++)
	{
#		print "$data[$i]\n";
		my $sam1_num=$i-3;
		for (my $j=3;$j<@data;$j++)
		{
			my $sam2_num=$j-3;
#			print "$data[$j]\n";die;
			if ($data[$i] ne $data[$j]&&$data[$i]!~/N/&&$data[$j]!~/N/)
			{
				$snp_hash{$sample_num[$sam1_num]}{$sample_num[$sam2_num]}++;
			}
			if ($sample_num[$sam1_num] eq $sample_num[$sam2_num])
			{
				$snp_hash{$sample_num[$sam1_num]}{$sample_num[$sam2_num]}=0;
			}
		}
	}
}
print OUT " \t";
foreach my $k1 (sort keys %snp_hash)
{
	print OUT "$k1\t"
}
print OUT "\n";
foreach my $sam1 (sort keys %snp_hash)
{
	print OUT "$sam1\t";
	foreach my $sam2(sort keys %{$snp_hash{$sam1}})
	{
		print OUT "$snp_hash{$sam1}{$sam2}\t";
	}print OUT "\n";
}

close (IN) ;
close OUT;
#######################################################################################

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
  -i <file>  input file,forced :eg :xxx/SNP/result/cucumber.formal.snp.vqsr.filter.snp
  
  -o <file>  output file,forced  
 
  -h         Help

USAGE
	print $usage;
	exit;
}
