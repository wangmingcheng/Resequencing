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
mkdir "$od/Cog_Anno" unless -d "$od/Cog_Anno";

my %Class=(
	"J" => [1,"Translation, ribosomal structure and biogenesis"],
	"A" => [2,"RNA processing and modification"],
	"K" => [3,"Transcription"],
	"L" => [4,"Replication, recombination and repair"],
	"B" => [5,"Chromatin structure and dynamics"],
	"D" => [6,"Cell cycle control, cell division, chromosome partitioning"],
	"Y" => [7,"Nuclear structure"],
	"V" => [8,"Defense mechanisms"],
	"T" => [9,"Signal transduction mechanisms"],
	"M" => [10,"Cell wall/membrane/envelope biogenesis"],
	"N" => [11,"Cell motility"],
	"Z" => [12,"Cytoskeleton"],
	"W" => [13,"Extracellular structures"],
	"U" => [14,"Intracellular trafficking, secretion, and vesicular transport"],
	"O" => [15,"Posttranslational modification, protein turnover, chaperones"],
	"C" => [16,"Energy production and conversion"],
	"G" => [17,"Carbohydrate transport and metabolism"],
	"E" => [18,"Amino acid transport and metabolism"],
	"F" => [19,"Nucleotide transport and metabolism"],
	"H" => [20,"Coenzyme transport and metabolism"],
	"I" => [21,"Lipid transport and metabolism"],
	"P" => [22,"Inorganic ion transport and metabolism"],
	"Q" => [23,"Secondary metabolites biosynthesis, transport and catabolism"],
	"R" => [24,"General function prediction only"],
	"S" => [25,"Function unknown"],
);

my %DEG;
open (IN,"$deg") or die $!;
while (<IN>) {
	next if /^\#/;
	my $name=(split/\t/,$_)[0];
	$DEG{$name}=1;
}
close IN;


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
	for my $row ( $row_min .. $row_max ) {
		my $now_name = $worksheet->get_cell( $row, 0 );
		$now_name=$now_name->value();
		next unless $now_name;
		next unless exists $DEG{$now_name};
		my $cell = $worksheet->get_cell( $row, 1 );
		next unless $cell;
		my $val=$cell->value();
		my ($class)=$1 if($val=~/(\[[A-Z]+\])/);
		next if(!defined $class);
		$class=~s/[\[\]]//g;
		foreach my $elem (split "",$class) {
			$Data{$elem}++;
		}
	}
}

open (OUT,">$od/Cog_Anno/$key.Cog.classfy.stat") || die "can't open file  $od/Cog_Anno/$key.Cog.classfy.stat\n ";
&OutHash(\%Data);
close(OUT);

#`ssh compute-0-14 /share/nas2/genome/biosoft/R/3.1.1/bin/Rscript $Bin/cog_anno_plot.r $od/Cog_Anno/$key.Cog.classfy.stat $od/Cog_Anno/$key.Cog.classfy.png `;
`/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript $Bin/cog_anno_plot.r $od/Cog_Anno/$key.Cog.classfy.stat $od/Cog_Anno/$key.Cog.classfy.png `;


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
	#????????????????
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
	#????????????????
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
	#????????????????????????????????????????????????ATTCCC->GGGAAT
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

sub OutHash
{
	my ($hash)=@_;
	print OUT "#ID\tClass_Name\tNumbers\n";
	foreach my $sty (sort {$Class{$a}[0] <=> $Class{$b}[0]} keys %Class)
	{
		print OUT "$sty\t"."$Class{$sty}[1]\t".(${$hash}{$sty}||0)."\n";
	}
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

  -k    <str>   keywords of output file,forced 
  
  -od   <dir>   output dir,forced 
  
  -h         Help

USAGE
	print $usage;
	exit;
}
