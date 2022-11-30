#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($cfg,$fOut,$unit,$key);
GetOptions(
	"help|?" =>\&USAGE,
	"o:s"=>\$fOut,
	"cfg:s"=>\$cfg,
	"unit:s"=>\$unit,
	"key:s"=>\$key
) or &USAGE;
&USAGE unless ($cfg and $fOut);
my %cfg_hash;
my %snp_hash;
my %indel_hash;
$unit||=100000;
open (IN,$cfg) or die $!;

while (<IN>) 
{
	chomp;
	next if ($_=~/\#/);
	next if ($_=~/^\s*$/);
	my @data=split/\s+/,$_;
	$cfg_hash{$data[0]}=$data[1] if ($data[0]=~/^Ref/);
	$cfg_hash{$data[0]}=$data[1] if ($data[0]=~/^Chr/);
	$cfg_hash{$data[0]}=$data[1] if ($data[0]=~/^Circle_type/);
	$cfg_hash{$data[0]}=$data[1] if ($data[0]=~/^SNPvcf/);
	$cfg_hash{$data[0]}=$data[1] if ($data[0]=~/^indelvcf/);
	$cfg_hash{$data[0]}=$data[1] if ($data[0]=~/^SVinfo/);
	$cfg_hash{$data[0]}=$data[1] if ($data[0]=~/^CNVinfo/);

}
close(IN);
open (IN,$cfg_hash{SNPvcf}) or die $!;
while (<IN>)
{
	chomp;
	next if ($_=~/^\#/);
	next if ($_=~/^\s*$/);
	my @data=split/\t/,$_;
	my $areanum=int ($data[1]/$unit)+1;
	$snp_hash{$data[0]}{$areanum}++;
}
close(IN);
open (OUT,">$fOut/snpdensity.list") or die $!;
my @ref=split/\,/,$cfg_hash{Chr};
for (my $i=0;$i<@ref ;$i++)
{
	foreach my $area (keys %{$snp_hash{$ref[$i]}})
	{
		my $start=($area-1)*$unit;
		my $end=($area)*$unit-1;
		my $density=$snp_hash{$ref[$i]}{$area}/$unit;
		print OUT "$ref[$i]\t$start\t$end\t$density\n";
	}
}
close(OUT);
open (IN,$cfg_hash{indelvcf}) or die $!;
while (<IN>)
{
	chomp;
	next if ($_=~/^\#/);
	next if ($_=~/^\s*$/);
	my @data=split/\t/,$_;
	my $areanum=int ($data[1]/$unit)+1;
	$indel_hash{$data[0]}{$areanum}++;
}
close(IN);


open (OUT,">$fOut/indeldensity.list") or die $!;
for (my $i=0;$i<@ref ;$i++)
{
	foreach my $area (keys %{$indel_hash{$ref[$i]}})
	{
		my $start=($area-1)*$unit;
		my $end=($area)*$unit-1;
		my $density=$indel_hash{$ref[$i]}{$area}/$unit;
		print OUT "$ref[$i]\t$start\t$end\t$density\n";
	}
}
close(OUT);
open (INS,">$fOut/svheatmap.INS.list") or die $!;
open (DEL,">$fOut/svheatmap.DEL.list") or die $!;
open (INV,">$fOut/svheatmap.INV.list") or die $!;
open (ITXCTX,">$fOut/svheatmap.ITXCTX.list") or die $!;

if(exists $cfg_hash{SVinfo}){
	open (IN,$cfg_hash{SVinfo}) or die "Can't open $cfg_hash{SVinfo}, $!\n" ;
	while (<IN>)
	{
		chomp;
		next if ($_=~/^\#/);
		next if ($_=~/^\s*$/);
		my @data=split/\s+/,$_;
		for (my $j=0;$j<@ref ;$j++)
		{
			my$s=&max($data[1],$data[4]);
			my$e=&min($data[1],$data[4]);
			
			if ($ref[$j]=~/$data[0]/&&$data[6]=~/INS/)
			{
				if($data[0] eq $data[3]){print INS "$data[0]\t$s\t$e\n" ;}
			}
			if ($ref[$j]=~/$data[0]/&&$data[6]=~/DEL/)
			{
				if($data[0] eq $data[3]){print DEL "$data[0]\t$s\t$e\n";}
			}
			if ($ref[$j]=~/$data[0]/&&$data[6]=~/INV/)
			{
				if($data[0] eq $data[3]){print INV "$data[0]\t$s\t$e\n";}
			}
			if ($ref[$j]=~/$data[0]/&&$data[6]=~/ITX|CTX/)
			{
				print ITXCTX "$data[0]\t$data[1]\t$data[1]\t$data[3]\t$data[4]\t$data[4]\n";;
			}
		}
	}
	close (ITXCTX);
	close (INS);
	close (DEL);
	close (INV);
}
if (exists $cfg_hash{CNVinfo}){
	open (IN,"$cfg_hash{CNVinfo}") or die "Can't open $cfg_hash{CNVinfo}, $!\n" ;
	open (CNV,">$fOut/CNV.list") or die $!;
	while (<IN>){
		chomp;
		next if ($_=~/^\#/);
		next if ($_=~/^\s*$/);
		my @data=split/\s+/,$_;
		print CNV "$data[0]\t$data[1]\t$data[2]\n";
	}
	close(CNV);
	if(exists $cfg_hash{SVinfo}){
		system(qq(perl $Bin/circos_v1.pl --chr $fOut/karyotype.txt --circle $fOut/snpdensity.list --type line --circle $fOut/indeldensity.list --type line --circle $fOut/CNV.list --type highlight --circle $fOut/svheatmap.INS.list --type highlight --circle $fOut/svheatmap.DEL.list --type highlight --circle $fOut/svheatmap.INV.list --type highlight   --circle $fOut/svheatmap.ITXCTX.list --type link --circle_width 0.05 --od $fOut --outfile $key.circos));
	}else{
		system(qq(perl $Bin/circos_v1.pl --chr $fOut/karyotype.txt --circle $fOut/snpdensity.list --type line --circle $fOut/indeldensity.list --type line --circle $fOut/CNV.list --type highlight  --circle_width 0.1 --od $fOut --outfile $key.circos));
		
	}
}else{
	if(exists $cfg_hash{SVinfo}){
		system(qq(perl $Bin/circos_v1.pl --chr $fOut/karyotype.txt --circle $fOut/snpdensity.list --type line --circle $fOut/indeldensity.list --type line  --circle $fOut/svheatmap.INS.list --type highlight --circle $fOut/svheatmap.DEL.list --type highlight --circle $fOut/svheatmap.INV.list --type highlight   --circle $fOut/svheatmap.ITXCTX.list --type link --circle_width 0.05 --od $fOut --outfile $key.circos));
	}else{
		system(qq(perl $Bin/circos_v1.pl --chr $fOut/karyotype.txt --circle $fOut/snpdensity.list --type line --circle $fOut/indeldensity.list --type line   --circle_width 0.1 --od $fOut --outfile $key.circos));
		
	}	
}

#######################################################################################
sub max{
        my $max=shift;
        my $temp;
        while (@_) {
                $temp=shift;
                $max=$max>$temp?$max:$temp;
        }
        return $max;
}



sub min{
        my $min=shift;
        my $temp;
        while (@_) {
                $temp=shift;
                $min=$min<$temp?$min:$temp;
        }
        return $min;
}



sub USAGE {#
	my $usage=<<"USAGE";

Usage:

  -cfg <file>  input file,fasta format,forced
  
  -unit <file>  input file,fasta format,forced
  
  -o <file>  output file,forced  
 
  -h         Help

USAGE
	print $usage;
	exit;
}
