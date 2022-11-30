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
my ($detail,$data,$od);
GetOptions(
	"help|?" =>\&USAGE,
	"detail:s"=>\$detail,
	"data:s"=>\$data,
	"od:s"=>\$od,
) or &USAGE;
&USAGE unless ($od or $detail or $data);
my %cfg_hash;
open (IN,$detail) or die $!;
while (<IN>){
	chomp;
	next if /^$/ || /^\#/;
	next if /^\%/;
	my @data=split /\s+/,$_;
	$cfg_hash{$data[0]}=$data[1];
}
close(IN);
my $chrnum = 0;
my $snpEff = $cfg_hash{snpEff};
my $project = $cfg_hash{Project};
my $java = "$cfg_hash{java}";
my $ref = $cfg_hash{ref};
my $gff = $cfg_hash{gff};
my $samtools = $cfg_hash{samtools};
###########################################################################################################
my %data_hash;
my $sam_num=0;
open (Rawdata,"$od/DataAssess/AllSample_GC_Q.stat") or die $!;
while (<Rawdata>){
	chomp;
	next if ($_=~/\#/);
	next if ($_=~/^\s*$/);
	next if ($_=~/^SampleID/);
	$sam_num++;
	my ($sample,$BaseSum,$Q30)=(split/\s+/,$_)[0,2,7];
	open (READ,$data) or die $!;
	while (my $reads=<READ>){
		chomp($reads);
		next if ($reads=~/\#/);
		next if ($reads=~/^\s*$/);
		my @rd=split/\s+/,$reads;
		if($rd[0] =~/^Sample/){
			my $sam = $rd[1];
			my $read1 = <READ>;
			chomp $read1;
			my @rd1=split/\s+/,$read1;
			$data_hash{$sam}{$rd1[0]} = $rd1[1];
			my $read2 = <READ>;
			chomp $read2;
			my @rd2=split/\s+/,$read2;
			$data_hash{$sam}{$rd2[0]} = $rd2[1];				
		}
	}
	close(READ);
}
close(Rawdata);
	
######################################################################################################################
my $ana_sam=0;
my $leninfo;
foreach my $samp (sort keys %data_hash){
	$ana_sam++;
}
if ($ana_sam!=$sam_num){
	print "the sample is not enough\n";
}

mkdir("$od/Analysis",0755) unless -d "$od/Analysis";
mkdir("$od/Analysis/Ref",0755) unless -d "$od/Analysis/Ref";
mkdir("$od/Analysis/Ref/$project",0755) unless -d "$od/Analysis/Ref/$project";
my $base = "$od/Analysis/Ref/$project";
`cp $ref $base/sequences.fa`;
`cp $gff $base/genes.gff`;
`cp $Bin/snpEff.config $od/Analysis/Ref/`;
`cp $ref $od/Analysis/Ref/sequences.fa`;
`cp $gff $od/Analysis/Ref/genes.gff`;
`$samtools faidx $base/sequences.fa`;
`$samtools dict $base/sequences.fa > $base/sequences.dict`;
`echo "data_dir = $od/Analysis/Ref" >> $od/Analysis/Ref/snpEff.config`;
`echo $project.genome:$base/sequences.fa >> $od/Analysis/Ref/snpEff.config`;
my $cmd = "$java -Xmx10240m -XX:ParallelGCThreads=5 -jar $snpEff build -gff3 -c $od/Analysis/Ref/snpEff.config -v $project >$base/build.log 2>$base/warning.log";
&run_or_die($cmd) if !-e "$base/snpEffectPredictor.bin";

$cmd = "perl $Bin/reflen.pl -fa $ref -o $od/Analysis/Ref/sequences.fa.len";
&run_or_die($cmd);
my $len = "$od/Analysis/Ref/sequences.fa.len";
open LEN,$len or die $!;
while(<LEN>){
	chomp;
	next if /^$/;
	my ($aa,$bb)= split /\s+/;
	if (/\#Total_Len/){
		$leninfo = $bb;
	}
}
close LEN;

open(ID, "$cfg_hash{chr_id}");
while(<ID>){
	chomp;
	next if (/^$/);
	next if (/^\#/);
	$chrnum++;
}
#######################################################################################################
open (CON,">$od/main.conf") or die $!;
my $type =<<"USAGE";
Ref1	$ref
Ref2	$od/Analysis/Ref/ref.genome.fa
Gff1	$gff
Gff2	$od/Analysis/Ref/ref.genome.gff
chr_id	$cfg_hash{chr_id}

assemble_level	chr
RefLen	$leninfo
ChrNum	$chrnum
Sample	$ana_sam
Quality		33
platform	$cfg_hash{platform}
maxcpu		$cfg_hash{maxcpu}
maxnum		$cfg_hash{maxnum}
queue		$cfg_hash{queue}
link		$cfg_hash{link}
ud		$cfg_hash{ud}
type		$cfg_hash{type}
Ploidy	$cfg_hash{ploidy}
delete  $cfg_hash{delete}

picard_dir	$cfg_hash{picard_dir}
samtools	$cfg_hash{samtools}
bwa		$cfg_hash{bwa}
gatk_dir	$cfg_hash{gatk_dir}
vcfutils	$cfg_hash{vcfutils}
freec		$cfg_hash{freec}
svg2xxx		$cfg_hash{svg2xxx}
snpEff		$cfg_hash{snpEff}
circos_dir	$cfg_hash{circos_dir}
java		$cfg_hash{java}

USAGE
print CON "$type";

foreach my $sam (sort keys %data_hash){
	my $othersam=$sam;
	$othersam=~/(.)(\d*)/;
	unless ($1=~/R/){
		if ($2<10){
			$othersam="R0".$2;
		}else{
			$othersam="R".$2;
		}
	}
	foreach my $read (sort keys %{$data_hash{$sam}}){
		if ($read=~/fq1/){
			my $pre=$othersam."-1";
			print CON "$pre\t$data_hash{$sam}{$read}\n";
		}if ($read=~/fq2/){
			my $pre=$othersam."-2";
			print CON "$pre\t$data_hash{$sam}{$read}\n";
		}
	}
}
my $type1=<<"USAGE1";
Species		$project
USAGE1
print CON "$type1";


sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub show_log()
{
        my ($txt) = @_ ;
        my $time = time();
        my $Time = &sub_format_datetime(localtime($time));
        print "$Time:\t$txt\n" ;
        return ($time) ;
}

sub run_or_die()
{
        my ($cmd) = @_ ;
        &show_log($cmd);
        my $flag = system($cmd) ;
        if ($flag != 0){
                &show_log("Error: command fail: $cmd");
                exit(1);
        }
        &show_log("done.");
        return ;
}


sub USAGE {#
	my $usage=<<"USAGE";
Usage:
  -data <file>  input file,fasta format,forced 
  -detail <file>
  -od 	       output
  -h         Help

USAGE
	print $usage;
	exit;
}



