#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $Time_Start = &sub_format_datetime(localtime($BEGIN_TIME));
print "Program Starts Time:$Time_Start\n";
my $version="1.0.0";
#######################################################################################

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
my %vcf_hash;
open (IN,$fIn) or die $!;
open (OUT,">$fOut") or die $!;
#$/=">";
while (<IN>) 
{
	chomp;
	next if ($_=~/^\#\#/);
	my @data=split /\t/,$_;
	if ($_=~/^\#/)
	{
		print OUT "$data[0]\t$data[1]\t$data[3]\t$data[4]\t";
		for (my $i=9;$i<@data ;$i++)
		{
			print OUT "$data[$i]\t";
		}
		print OUT "\n";
	}
	else
	{
		$vcf_hash{$data[0]}{$data[1]}=$data[3]."\t".$data[4];
		for (my $j=9;$j<@data ;$j++)
		{
			my $type=(split/\:/,$data[$j])[0];
			$vcf_hash{$data[0]}{$data[1]}.="\t".$type;
		}
	}
}
foreach my $sca (sort keys %vcf_hash)
{
	foreach my $pos (sort keys %{$vcf_hash{$sca}})
	{
		print OUT "$sca\t$pos\t$vcf_hash{$sca}{$pos}\n";
	}
}








close (IN) ;
close (OUT) ;
#######################################################################################
my $Time_End   = sub_format_datetime(localtime(time()));
print STDOUT "Program Ends Time:$Time_End\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#######################################################################################

sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

#######################################################################################

sub max{#&max(lists or arry);
	#求列表中的最大值
	my $max=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$max=$max>$temp?$max:$temp;
	}
	return $max;
}

#######################################################################################

sub min{#&min(lists or arry);
	#求列表中的最小值
	my $min=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$min=$min<$temp?$min:$temp;
	}
	return $min;
}

#######################################################################################

sub revcom(){#&revcom($ref_seq);
	#获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return uc $seq;			  
}

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
ProgramName:
Version:	$version
Contact:	Wujianhu <wujh\@biomarker.com.cn> 
Program Date:   2014.7.10
Description:	this program is used to extract indel result
Usage:
  Options:
  -i <file>  input file,forced : eg :xxx/Annotation/indel_anno/cucumber.formal.indel.vqsr.filter.anno.gatk.vcf 
  
  -o <file>  output file,forced  
 
  -h         Help

USAGE
	print $usage;
	exit;
}
