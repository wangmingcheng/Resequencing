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
my ($cfg,$key,$fOut,$detail);
GetOptions(
	"help|?" =>\&USAGE,
	"o:s"=>\$fOut,
	"cfg:s"=>\$cfg,
	"key:s"=>\$key,
	"detail:s"=>\$detail,
) or &USAGE;
&USAGE unless ($cfg and $fOut and $key and $detail);

my %hash;
open (IN,$detail) or die $!;
while (<IN>) {
	chomp;
	next if /^$/ || /^\#/;
	my @data=split /\s+/,$_;
	$hash{$data[0]}=$data[1];
}
close IN;
my $circos_dir = $hash{circos_dir} ;

my %cfg_hash;
my %hash_ref;
my $total_len;

open (IN,$cfg) or die $!;
open (OUT,">$fOut/karyotype.txt") or die $!;

while (<IN>) 
{
	chomp;
	next if ($_=~/^\#/);
	next if ($_=~/^\s*$/);
	my @data=split/\s+/,$_;
	$cfg_hash{$data[0]}=$data[1] if ($data[0]=~/^Ref/);
	$cfg_hash{$data[0]}=$data[1] if ($data[0]=~/^Chr/);
	$cfg_hash{$data[0]}=$data[1] if ($data[0]=~/^Circle_type/);
}
close(IN);
%hash_ref = %{&Getrefinfo($cfg_hash{Ref})};
my @chrtype=split/\,/,$cfg_hash{Chr};
for (my $i=0;$i<@chrtype;$i++)
{
	if (exists($hash_ref{$chrtype[$i]}))
	{
		my $num=$i+1;
		print OUT "chr\t-\t$chrtype[$i]\t$chrtype[$i]\t0\t$hash_ref{$chrtype[$i]}\tchr$num\n";
	}
}
foreach my $len (keys %hash_ref)
{
	$total_len+=$hash_ref{$len};
}
close(OUT);

my $unit=int($total_len/3000);

print "perl $Bin/extract_result2circos.pl -cfg $cfg -o $fOut -unit $unit -key $key\n";
system("perl $Bin/extract_result2circos.pl -cfg $cfg -o $fOut -unit $unit -key $key");
close (IN) ;
close (OUT) ;

##############################################################################################################
sub Getrefinfo{
	my %ref_hash;
	my @files = @_ ;
	for (@files)
	{
		$/=">";
		open (IN, $_) or die "Can't open $_, $!\n";
		{
			while (<IN>)
			{
				chomp;
				next if ($_=~/^\s*$/);
				my ($id,$seq)=split/\n/,$_,2;
				my $real_id=(split/\s+/,$id)[0];
				$seq=~s/\r//g;
				my $chr_len=length($seq);
				$ref_hash{$real_id}=$chr_len;
			}
		}
	}
	return (\%ref_hash);
}


sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}



sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";

Usage:

  -cfg <file>  input file,fasta format,forced 
  
  -key <file>  input file,fasta format,forced

  -o <file>  output file,forced  
  
  -detail  <file> 
 
  -h         Help

USAGE
	print $usage;
	exit;
}
