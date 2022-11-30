#!/usr/bin/perl -w
# 
# Copyright (c) BIO_MARK 2014
# Writer:         Dengdj <dengdj@biomarker.com.cn>
# Program Date:   2014
# Modifier:       Dengdj <dengdj@biomarker.com.cn>
# Last Modified:  2014.
my $ver="1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

####################

my %opts;
GetOptions(\%opts,"i=s","o=s","s=s","r=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{i}) || !defined($opts{s}) || !defined($opts{r}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		
		Version: $ver

	Usage:

		-i           vcf file                    <infile>         must be given

		-o           out file                    <outfile>        must be given

		-r           ref file                    <infile>         must be given

		-s           species                     <string>         must be given

		-h           Help document

	Usage End.

	exit;
}

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
################

# get parameter
my $vcffile = $opts{i} ;
my $outfile = $opts{o} ;
my $outdir = dirname($outfile);
mkdir $outdir ;
$outdir = &ABSOLUTE_DIR($outdir) ;
$outfile = "$outdir/".basename($outfile);
my $species = $opts{s} ;
my $reffile = $opts{r} ;
my $java_dir = "/share/nas2/genome/biosoft/java/jdk/current/bin/" ;
my $gatk = "/share/nas2/genome/biosoft/GenomeAnalysisTK/current/GenomeAnalysisTK.jar" ;
my $snpEff = "/share/nas2/genome/biosoft/snpEff/current/snpEff.jar" ;

# annotation for variants
my $cmd = "cd $outdir && java -XX:ParallelGCThreads=5 -jar $snpEff $species -v $vcffile -o gatk > $outfile" ;
print "$cmd\n" ;
print `$cmd` ;
(my $gatk_anno = $outfile) =~ s/.vcf$// ;
$gatk_anno .= ".gatk.vcf" ;
$cmd = "$java_dir/java -XX:ParallelGCThreads=5 -jar $gatk -T VariantAnnotator -R $reffile -A SnpEff --variant $vcffile --snpEffFile $outfile -L $vcffile -o $gatk_anno" ;
print "$cmd\n" ;
print `$cmd` ;


###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";

###############Subs
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;
	$cur_dir =~ s/\n$//;
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;
		$dir =~ s/\n$// ;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;
		$return =~ s/\n$// ;
	}
	else
	{
		warn "Warning just for file and dir\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}

