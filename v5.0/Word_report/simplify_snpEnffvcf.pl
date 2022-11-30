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
my ($snpEnff,$fOut);
GetOptions(
	"help|?" =>\&USAGE,
	"o:s"=>\$fOut,
	"snpEnff:s"=>\$snpEnff,
) or &USAGE;
&USAGE unless ($snpEnff and $fOut);
my $dir=dirname($snpEnff);
my %snpEnff_vcf;
open (IN,$snpEnff) or die $!;
open (OUT,">$dir/mid_snpEff.vcf") or die $!;
while (<IN>) 
{
	chomp;
	next if ($_=~/^\#\#/);
	if ($_=~/^\#/) 
	{
		print OUT "$_\n";
	}
	else
	{
		my @info=split/\t/,$_;
		$snpEnff_vcf{$info[0]}{$info[1]}=$info[3]."\t".$info[4]."\t".$info[5];
		for (my $i=9;$i<@info;$i++)
		{
			$snpEnff_vcf{$info[0]}{$info[1]}.="*".$info[$i];
		}
		my @anno=split/;/,$info[7];
		for (my $i=0;$i<@anno;$i++)
		{
			if ($anno[$i]=~/SNPEFF_EFFECT|SNPEFF_TRANSCRIPT_ID/)
			{
				$snpEnff_vcf{$info[0]}{$info[1]}.=";".$anno[$i];
			}
		}
	}
}
foreach my $chr (keys %snpEnff_vcf)
{
	foreach my $pos (keys %{$snpEnff_vcf{$chr}})
	{
		if($snpEnff_vcf{$chr}{$pos}=~/=START_LOST|=NON_SYNONYMOUS_START|=STOP_GAINED|=STOP_LOST|=RARE_AMINO_ACID|=NON_SYNONYMOUS_CODING/)
		{
			print OUT "$chr\t$pos\t$snpEnff_vcf{$chr}{$pos};NON_SYNONYMOUS\n";
		}
		elsif ($snpEnff_vcf{$chr}{$pos}=~/=SYNONYMOUS_CODING|=SYNONYMOUS_START|=SYNONYMOUS_STOP/)
		{
			print OUT "$chr\t$pos\t$snpEnff_vcf{$chr}{$pos};SYNONYMOUS\n";
		}
		else
		{
			print OUT "$chr\t$pos\t$snpEnff_vcf{$chr}{$pos}\n";
		}
	}
}
close(IN);
close(OUT);

open (IN,"$dir/mid_snpEff.vcf") or die $!;
open (OUT,">$fOut") or die $!;
my $tital=<IN>;
chomp($tital);
my @header=split/\t/,$tital;
print OUT "$header[0]\t$header[1]\t$header[3]\t$header[4]\t$header[5]\t";
for (my $i=9;$i<@header;$i++)
{
	print OUT "$header[$i]_base\t$header[$i]_alt_depth\t$header[$i]_reads_depth\t";
	
}
print OUT "SNPEFF_EFFECT\n";
while (<IN>)
{
        chomp;
        next if ($_=~/\#/);
        my @unit=split/\t/,$_;
        my @use_info=split/\;/,$unit[4],2;
        my @format=split/\*/,$use_info[0];
        print OUT "$unit[0]\t$unit[1]\t$unit[2]\t$unit[3]\t";
        print OUT "$format[0]\t";
        my $alt=join"\,",@unit[2,3];
        my @alt_type=split/\,/,$alt;
        for (my $i=1;$i<@format ;$i++)
        {
                if ($format[$i]=~/\.\/\.|[.]/) 
                {
                        print OUT "N\tN\tN\t"
                }
                elsif ($format[$i]=~/\d\/\d|\d/)
                {
                        my ($sampletype,$altdepth,$readsdepth)=split/:/,$format[$i];
                        if($sampletype =~ /\//){
                                my ($ref_type_num,$sample_type_num)=split/\//,$sampletype;
                                print OUT "$alt_type[$ref_type_num],$alt_type[$sample_type_num]\t$altdepth\t$readsdepth\t";
                        }else{
                                print OUT "$sampletype\t$altdepth\t$readsdepth\t";
                        }
                }
        }
        if ($use_info[1]) 
        {
                print OUT "$use_info[1]\n";
        }
        else
        {
                print OUT "\n";
        }
}
`rm $dir/mid_snpEff.vcf`;
close (IN) ;
close (OUT) ;

#######################################################################################


sub USAGE {#
	my $usage=<<"USAGE";
Usage:
  -snpEnff <file>  input file,snpEnff.vcf,forced
  
  -o <file>  output file,forced  
 
  -h         Help

USAGE
	print $usage;
	exit;
}
