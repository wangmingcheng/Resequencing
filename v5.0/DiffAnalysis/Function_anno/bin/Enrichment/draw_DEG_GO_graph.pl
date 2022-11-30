#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Encode;
use Cwd qw(abs_path);
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
$od=abs_path($od);
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
	if ($True_sheet!~/^GO.list/) {
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


my $mid_file_tree="$od/GO_tree.stat_File_0";
for (my $i=0; ;$i++) {
	if (-f $mid_file_tree) {
		my $now=$i+1;
		$mid_file_tree="$od/GO_tree.stat_File_".$now;
		next;
	}
	last;
}

for my $worksheet ( $workbook->worksheets() ) {
	my $True_sheet = $worksheet->get_name();
	if ($True_sheet!~/^GO_tree\.stat/) {
		next;
	}
	open (OUT,">$mid_file_tree") or die $!;
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

if (-f "$mid_file_list" && -f "$mid_file_tree") {
    print "perl $Bin/KeggGo_enrich_map_web.pl -d $deg -k $key -i $od -o $od -func go\n";
	`perl $Bin/KeggGo_enrich_map_web.pl -d $deg -k $key -i $od -o $od -func go`;
	`rm $mid_file_list`;
	`rm $mid_file_tree`;
}



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

################################################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
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
