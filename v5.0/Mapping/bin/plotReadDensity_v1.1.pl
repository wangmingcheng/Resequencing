#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
#===============================================================

my $fBam;
my $input_file_format = "depth"; ## can be sam, bam or depth

my $fKey;
my $dOut;

my $chrlistfile ;
my $assemble_level;
my $x_title;
my $y_title;
my $title;

my $window = 1000;

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fBam,
				"f:s"=>\$input_file_format,
				"o:s"=>\$dOut,
				"k:s"=>\$fKey,
				"l:s"=>\$assemble_level,
				"c:s"=>\$chrlistfile,

				"x:s"=>\$x_title,
				"y:s"=>\$y_title,
				"t:s"=>\$title,
				
		
				) or &USAGE;
&USAGE unless ($fBam and $fKey and $chrlistfile and $assemble_level);

$dOut||="./";
mkdir $dOut unless (-d $dOut) ;

$x_title||= "Chromsome position";
$y_title||= "Median of read deinsity(log(2))";
$title||= "Genomewide distribution of read coverage";

my $prefix = "$dOut/$fKey";
my $out_coverage_file = "$prefix.cov.txt";

open (OUT,">$out_coverage_file") or die $!;

if ($input_file_format eq 'depth') {
        open (IN,"<",$fBam) or die $!;
}

$/="\n";

my @tmp = ();
my $cur_pos_index = 0;
my $last_chr = "";
my $total_chro = 0;

while (<IN>) {
	chomp;
	s/\r//g;
	next if (/^$/ || /^;/ || /^\#/) ;
	
	my ($chr, $pos, $depth) = split;
	
	if ($chr ne $last_chr) {
		## output last chromosome's info any more
		if ($last_chr ne "") {
			my $midian_cov = mid(\@tmp);

			print OUT join("\t",
				$last_chr,
				$cur_pos_index,
				log($midian_cov+1)/log(2),
			),"\n";
		}

		$cur_pos_index = 1;
		@tmp = ();

		push @tmp, $depth;
		$last_chr = $chr;
		$total_chro++;
		
		next;
	}

	if ($pos > $cur_pos_index * $window) {

		my $midian_cov = mid(\@tmp);
		
		print OUT join("\t",
			$chr,
			$cur_pos_index,
			log($midian_cov+1)/log(2),
		),"\n";

		while ((++$cur_pos_index) * $window < $pos) {
			
			print OUT join("\t",
				$chr,
				$cur_pos_index,
				0,
			),"\n";
		}

		@tmp = ();
	}
	
	push @tmp, $depth;

	$last_chr = $chr;
	
}
close (IN) ;
close (OUT) ;

#------------------------------------------------------------------
# plot calling r scripts 
#------------------------------------------------------------------


my $cmd = "Rscript $Bin/genomeCoveragehorizontalArea.r --infile $out_coverage_file --idfile $chrlistfile --outfile $out_coverage_file.png --group.col 1 --x.col 2 --y.col 3 --x.lab \"$x_title\" --y.lab \"$y_title\" --title.lab \"$title\" --unit kb --level $assemble_level";

&process_cmd($cmd);

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub mid {
	my $arr = shift;

	return 0 if (@$arr == 0) ;
	my @list = sort @$arr;
	return $list[(@list-@list%2)/2]; 
} 

sub process_cmd {#
	my ($cmd) = @_;

	print STDERR "CMD: $cmd\n";
	my $start_time = time;

	my $ret = system($cmd);
	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}

	my $end_time = time;
	print STDERR "CMD Done, elapsed time ( ", $end_time - $start_time, "s )\n\n";

	return;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Usage:
  Options:
  -i <file>  input sam/bam/depth file, result of samtools, forced.
  -f <str>   input file format, depth|sam|bam, default [depth].
  -l <str>   reference assemble level(chr or scaffold)
  
  -c <int>   chr name list file, contain chr id need to draw, forced.
  -k <str>   key of output file , forced.
  -o <dir>   direction where restult produced.
  
  -x <str>   title of x axis, optional.
  -y <str>   titile of y axis, optional.
  -t <str>   titile of graph, optional.
  
USAGE
	print $usage;
	exit;
}

