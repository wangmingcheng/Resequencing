#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($data,$detail,$od);
GetOptions(
	"help|?" =>\&USAGE,
	"data:s"=>\$data,
	"detail:s"=>\$detail,
	"od:s" =>\$od,
) or &USAGE;

&USAGE unless ($data or $detail or $od);

unless (-d "$od/DataAssess"){
    `mkdir -p $od/DataAssess`;
}
open (IN,$data) or die $!;
my %seq_hash;

my $sample;
while (<IN>){
    chomp;
    next if /^$/ || /^\#/;
    my @data=split /\s+/,$_;
    if ($data[0]=~/Sample/i){
        $sample=$data[1];                  
    }
    if ($data[0]=~/^fq1/){
        $seq_hash{$sample}{1}=$data[1];                            
    }
    if ($data[0]=~/^fq2/){
        $seq_hash{$sample}{2}=$data[1];
    }
}
open (FASTP,">$od/DataAssess/fastp.sh") or die $!;
foreach my $sample (sort keys %seq_hash){
    chomp $sample;
    `mkdir -p $od/DataAssess/$sample` unless (-d "$od/DataAssess/$sample");
    print FASTP "cd $od/DataAssess/$sample && fastp -i $seq_hash{$sample}{1} -I $seq_hash{$sample}{2} -o $sample\_clean.1.fq.gz -O $sample\_clean.2.fq.gz --thread 7 --json $sample.json --html $sample.html\n";
}
system ("parallel -j 10 <$od/DataAssess/fastp.sh");
##################################################
sub USAGE {
        my $usage=<<"USAGE";

Usage:
  -data		data config 
  -detail	detail parameter
  -od		output directory	
  -h         Help
    
USAGE
        print $usage;
        exit;
}
