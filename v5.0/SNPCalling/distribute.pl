#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin $Script);
my (%gene,%sum);
my (@gene_id,@group);
my $i = 0;
my ($In,$Outfile,$dis_num);
GetOptions(
			"i:s"=>\$In,
			"o:s"=>\$Outfile,
			"n:i"=>\$dis_num);
if (!defined $In || !defined $Outfile || !defined $dis_num){
	print "\n\tUsage:\tperl $0 -i fa_length_file -o out_prefix -n number\n\n" ;
	exit(1);
}
open IN,"$In" or die $!;
while (<IN>) {
	next if (/^\#/); 
	my @string = split/\t/;
	$gene{$string[0]}=$string[1];
	push @gene_id,$string[0];
}
close IN;
my @f=glob("$Outfile.*.list");
if(@f >0 ){
	system("rm $Outfile.*.list");
}

my$s=0;
$i=1;

my$first=1;
foreach  (@gene_id) {

        if($gene{$_}>10000000 ){
                open OUT,">$Outfile.$i.list" or die $!;
                $i++;
                $s=0;
        }elsif($s>10000000 ){
                open OUT,">$Outfile.$i.list" or die $!;
                $i++;
                $s=0;
        }elsif($first and $gene{$_}<=10000000 and $i==1){
                open OUT,">$Outfile.$i.list" or die $!;
                $i++;
                $s=0;
                $first=0;

        }
        $s+=$gene{$_};
        print OUT "$_\n";



}
close OUT;

