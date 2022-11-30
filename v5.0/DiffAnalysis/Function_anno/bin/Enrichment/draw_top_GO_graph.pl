#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Encode;
use Spreadsheet::WriteExcel;
use Spreadsheet::ParseExcel;  
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$deg,$key,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"deg:s"=>\$deg,
				"k:s"=>\$key,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($fIn and $deg and $key and $od);

mkdir $od unless -d $od;
mkdir "$od/Graph" unless -d "$od/Graph";
my $parser   = Spreadsheet::ParseExcel->new();
my $workbook = $parser->parse("$fIn");

if ( !defined $workbook ) {
	die $parser->error(), ".\n";
}

my $mid_file_list="$od/GO.list_File_0";
for (my $i=0; ;$i++) {
	if (-f $mid_file_list) {
		my $now=$i+1;
		$mid_file_list="$od/GO.list_File_".$now;
		next;
	}
	last;
}

for my $worksheet ( $workbook->worksheets() ) {
	my $True_sheet = $worksheet->get_name();
	if ($True_sheet!~/^GO\.list/) {
		next;
	}
	open (OUT,">$mid_file_list") or die $!;
	my $Val;
	my ( $row_min, $row_max ) = $worksheet->row_range();
	my ( $col_min, $col_max ) = $worksheet->col_range();
	for my $row ( $row_min .. $row_max ) {
		for my $col ( $col_min .. $col_max ) {
			my $cell = $worksheet->get_cell( $row, $col );
			next unless $cell;
			$Val.=$cell->value()."\t";
		}
		$Val=~s/\t$/\n/;
	}
	print OUT $Val;
	close OUT;
}



if (-f "$mid_file_list") {
	`perl $Bin/draw_GO_DAG_map.pl -All_GO $mid_file_list -DEG_list $deg -od $od/Graph -key $key`;
#	`rm $mid_file_list`;
}



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
Contact:	Zhang XueChuan <zhangxc\@biomarker.com.cn> 
Program Date:   2013.10.17
Usage:
  Options:
  -i     <file>  input file,All_Database_annotation.xls,forced 
  
  -deg   <file>  deg file,forced 
  
  -k     <str>   keywords of output file,forced 
  
  -od    <file>  output dir,forced 
  
  -h         Help

USAGE
	print $usage;
	exit;
}
