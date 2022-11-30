#!/usr/bin/env perl
use warnings;
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);

my ($indelvcf,$od);
GetOptions(
	"i:s" => \$indelvcf,
	"o:s" => \$od,
    "h|?" => \&help,
) || &help;

&help unless ($indelvcf && $od );

sub help
{
        print <<"	Usage End.";
        
    Usage:
        -indelvcf 
        -od 
        
	Usage End.
	exit;
}

$indelvcf = abs_path($indelvcf);
mkdir($od,0755) unless -d $od;
$od = abs_path($od);

my $cmd = "cd $od && /share/nas2/genome/biosoft/plink/1.9/plink --vcf $indelvcf --genome full --allow-extra-chr --out indel\n";
system($cmd);
my %indel_hash;
open IN,"$od/indel.genome" or die $!;
open OUT,">$od/DEG_indel.stat" or die $!;
my %snp_hash;
while(<IN>){
	chomp;
	s/^\s+//g;
	next if /FID1/;
	my ($sam1,$sam2,$ibs0,$ibs1) = (split /\s+/)[0,2,14,15];
	my $diff = $ibs0 + $ibs1;
	$indel_hash{$sam1}{$sam2} = $diff;
	$indel_hash{$sam2}{$sam1} = $diff;
	$indel_hash{$sam1}{$sam1} = 0;
	$indel_hash{$sam2}{$sam2} = 0;
}
close IN;

print OUT " \t";
foreach my $k1 (sort keys %indel_hash){
	print OUT "$k1\t"
}
print OUT "\n";
foreach my $sam1 (sort keys %indel_hash){
	print OUT "$sam1\t";
	foreach my $sam2(sort keys %{$indel_hash{$sam1}}){
		if($sam1 eq $sam2){
			print OUT "0\t";
		}else{
			print OUT "$indel_hash{$sam1}{$sam2}\t";
		}
	}
	print OUT "\n";
}
close OUT;



