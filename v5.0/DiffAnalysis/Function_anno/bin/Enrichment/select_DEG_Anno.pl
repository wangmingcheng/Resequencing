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
my ($fIn,$deg,$o);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"deg:s"=>\$deg,
				"o:s"=>\$o,
				) or &USAGE;
&USAGE unless ($fIn and $deg and $o);


my $head_line;
my $limit=0;
my %DEG;
open (IN,"$deg") or die $!;
while (<IN>) {
	$_=~s/\s+$//;
	$head_line=$_ if /^\#/;
	next if /^\#/;
	my ($name,$val)=split/\t/,$_,2;
	$DEG{$name}=$val;
}
close IN;
open (OUT,">$o") or die $!;

my %Data;

my $parser   = Spreadsheet::ParseExcel->new();
my $workbook = $parser->parse("$fIn");

if ( !defined $workbook ) {
	die $parser->error(), ".\n";
}

for my $worksheet ( $workbook->worksheets() ) {
	my $True_sheet = $worksheet->get_name();
	if ($True_sheet!~/^Integrated_Function.anno/) {
		next;
	}
	my ( $row_min, $row_max ) = $worksheet->row_range();
	my ( $col_min, $col_max ) = $worksheet->col_range();
	for my $row ( $row_min .. $row_max ) {
		my $col_s=$col_min+1;
		if ($row==$row_min) {
			if ($limit==0) {
				foreach my $col ($col_s..$col_max) {
					my $cell = $worksheet->get_cell( $row, $col );
					next unless $cell;
					$head_line.="\t".$cell->value();
				}
				print OUT "$head_line\n";
				$limit++;
			}
		}
		else {
			my $gene=$worksheet->get_cell( $row, $col_min )->value();
			next unless exists $DEG{$gene};
			print OUT "$gene\t$DEG{$gene}";
			my $Val;
			for my $col ( $col_s .. $col_max ) {
				my $cell = $worksheet->get_cell( $row, $col );
				next unless $cell;
				$Val.="\t".$cell->value();
			}
			print OUT "$Val\n";
		}
	}
}

close(OUT);


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
  -i    <file>  input file,All_Database_annotation.xls,forced 
  
  -deg  <file>  deg file,forced

  -o    <file>  output file,forced 
  
  -h         Help

USAGE
	print $usage;
	exit;
}
