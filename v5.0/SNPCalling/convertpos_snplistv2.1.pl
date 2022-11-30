#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($markerPos,$vcf,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$od,
				"v:s"=>\$vcf,
				"pos:s"=>\$markerPos,

				) or &USAGE;
&USAGE unless ($vcf and $markerPos and $od );

mkdir $od if (! -d $od);
$vcf=&ABSOLUTE_DIR($vcf);
$markerPos=&ABSOLUTE_DIR($markerPos);
$od=&ABSOLUTE_DIR($od);

#my $file=basename($vcf);
open OUT,">$od/ConvertPos.vcf" or die $!;
open OUT2,">$od/ConvertPos.table" or die $!;
open IN,$vcf or die $!;
my %SNP;
while (<IN>) {
	chomp;
	next if (/^$/);
	if (/^\#/) {
		print OUT "$_\n";
	}
	else {
		my ($chr,$pos)=(split /\s+/,$_)[0,1];
		#print $chr;die;
		$SNP{$chr}{$pos}=$_;
	}
}
close IN;

open IN,$markerPos or die $!;
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my ($link,$id,$start,$len)=split /\s+/,$_;
	#print $link;die;
	for (my $i=$start;$i<=$start+$len ;$i++) {
		if (exists $SNP{$link}{$i}) {
			my ($id_ref,$pos_ref, $pos_id, $str)=split /\s+/,$SNP{$link}{$i},4;
			my $pos_marker=$pos_ref-$start;
			$pos_id = join("__", $id, $pos_marker);
			print OUT "$id\t$pos_marker\t$pos_id\t$str\n";
			print OUT2 "$id_ref\t$pos_ref\t$id\t$pos_marker\n";
		}
	}
}
close IN;
close OUT;


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################

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

################################################################################################################

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

################################################################################################################

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

################################################################################################################

sub revcom(){#&revcom($ref_seq);
	#获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return uc $seq;			  
}

################################################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	Chen Xiang <chenx\@biomarker.com.cn> 
Description:	convert SNP and INDEL position through Marker position in reference sequence
Usage:
  Options:
  -v	<file>	input file,vcf format,forced
  -pos	<file>	input file,marker.pos,forced
  -o	<dir>	output dir,forced
 
  -h         Help

USAGE
	print $usage;
	exit;
}
