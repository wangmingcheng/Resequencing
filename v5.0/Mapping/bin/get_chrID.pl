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
my ($fIn,$fOut,$fIn1);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"o:s"=>\$fOut,
				"n:s"=>\$fIn1,
				) or &USAGE;
&USAGE unless ($fIn and $fOut and $fIn1);
open (IN,$fIn) or die $!;
open (OUT,">$fOut") or die $!;

my $num;
my %chr_len;
my $key;
my $length;

while (<IN>) {
	chomp;
	if ($_ =~ /^>/) {
		$key = (split)[0] ;
		$key =~ s/>// ;
		$chr_len{$key} = 1;
		$length = 0;
		$num++;
	}else{
		$length += length($_);
		$chr_len{$key} = $length;
	}
}
 
my @keys = sort { $chr_len{$b} <=> $chr_len{$a} } keys %chr_len;


print "The chromosomes' number is :$num";

my %sort;

if ($fIn1 > $num) {
	print "Exceed the limit!!!";
}else{
	for (my $i = 0 ;$i <= ($fIn1 -1) ;$i++) {
		$sort{$keys[$i]} = $chr_len{$keys[$i]};
	}
}

my ($chr1, $chr2) ;
foreach my $chars (sort {(($chr1)=$a=~/(\d+)$/, $chr1) <=> (($chr2)=$b=~/(\d+)$/, $chr2)} keys %sort) {
	print OUT "$chars\t$sort{$chars}\n";
}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
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

sub GetTime {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:Cheng XinChao<chengxc\@biomarker.com.cn> 
Description:
Usage:
  Options:
  -i <input file>   
  -o <output file>
  -n <num>

  -h         Help

USAGE
	print $usage;
	exit;
}
