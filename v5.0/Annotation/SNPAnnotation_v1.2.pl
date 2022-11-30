#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);

my %opts;
GetOptions(\%opts,"i=s","o=s","r=s","od=s","mode=s","c=s","h" );

if(!defined($opts{i}) || !defined($opts{r}) || !defined($opts{od}) || !defined($opts{mode}) || !defined($opts{c}) || defined($opts{h}))
{
	print <<"	Usage End.";

	Usage:

		-i           vcf file                               <infile>                   must be given
		-od          outdir                                 <outdir>                   must be given
		-r           ref file                               <infile>                   must be given
		-mode        SNP/Indel                              <string>                   must be given
		-c	     config
		-h           Help document

	Usage End.

	exit;
}

#####################################################################################################
my $config = $opts{c};
my %hash;
open CFG,$config or die $!;
while(<CFG>){
        chomp;
        next if /^$/ || /^\#/;
        my ($key,$value) = split /\s+/;
        $hash{$key} = $value;
}
close CFG;

my $outprefix = $hash{Species};
my $ud = $hash{ud};
my $maxproc = $hash{maxcpu};
my $queue = $hash{queue};
my $java = $hash{java};
my $gatk = $hash{gatk_dir};
my $snpEff = $hash{snpEff};

my $vcffile = $opts{i} ;
$vcffile = abs_path($vcffile) ;
my $outdir = $opts{od} ;
mkdir $outdir unless(-d $opts{od});
$outdir = abs_path($outdir) ;
my $reffile = $opts{r} ;
$reffile = abs_path($reffile);
my $mode = $opts{mode} ;

my $min_length ||= -10 ;
my $max_length ||= 10 ;


my $shdir = "$outdir/sh_dir" ;
mkdir $shdir ;
my $anno_file = &annotation_for_variants($vcffile);

my @sample_files = ();
my @samples = ();
my $sample_list_file = "$outdir/../../Ref/sample.list" ;
$sample_list_file = abs_path($sample_list_file);
my $asamples = &get_sample_files(\@sample_files, $anno_file, $sample_list_file);
@samples = @{$asamples} ;

# result statistic
&statistic_annotation_result(\@sample_files, $anno_file, $mode);

############################################################################################
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
	my $start_time = &show_log($cmd);
	my $flag = system($cmd) ;
	if ($flag != 0){
		my $end_time = &show_log("Error: command fail: $cmd");
		&show_log("Elaseped time: ".($end_time - $start_time)."s\n");
		exit(1);
	}
	my $end_time = &show_log("done.");
	&show_log("Elaseped time: ".($end_time - $start_time)."s\n");
	return ;
}

sub qsub()
{
	my ($shfile) = @_ ;
	my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc $maxproc --queue $queue --reqsub $shfile --independent" ;
	&run_or_die($cmd);

	return ;
}

sub get_sample_files()
{
	my ($asample_files, $vcffile, $sample_list_file) = @_ ;
	my $shfile = "$shdir/get_sample_vcf.sh" ;
	my @samples = ();
	open (IN, $sample_list_file) || die "Can't open $sample_list_file, $!\n" ;
	open (SH, ">$shfile") || die "Can't creat $shfile, $!\n" ;
	while(<IN>){
		chomp ;
		next if (m/^\#/ || m/^\s*$/);
		my $sample = (split)[0] ;
		my $dir = "$outdir/$sample" ;
		mkdir "$dir" ;
		open OUT,">$dir/$sample.txt" or die $!;
		print OUT "$sample\n";
		close OUT;
		my $basename = basename($vcffile);
		(my $outfile = "$dir/$basename") =~ s/.vcf$/.$sample/ ;
		#print SH "$java -jar $gatk/GenomeAnalysisTK.jar -T SelectVariants -R $reffile --variant $vcffile -o $outfile -sn $sample -env \n" ;
		print SH "$Bin/vcftools --vcf $vcffile --recode --out $outfile --keep $dir/$sample.txt --recode-INFO-all --max-missing 1 --non-ref-ac 1\n";
		push @{$asample_files}, "$outfile.recode.vcf" ;
		push @samples, $sample ;
	}
	close(IN);
	close(SH);
	&qsub($shfile);
	return(\@samples) ;
}

sub annotation_for_variants()
{
	my ($vcffile) = @_ ;
	open (IN,"$vcffile") or die "Can't open $vcffile,$!\n";
	my $basename = basename($vcffile) ;
	(my $outfile = "$outdir/$basename") =~s/.vcf$/.anno.vcf/ ;
	my $cmd = "" ;
	if (!-f $outfile){
		$cmd = "ln -s $vcffile $outdir" ;
		&run_or_die($cmd);
	}
	(my $vcffile_filter = $basename)=~s/.vcf$/_filter.vcf/ ;
	open (VCF,">$outdir/$vcffile_filter") or die "Can't creat $vcffile_filter, $!\n";
	while (<IN>)
	{
		chomp;
		if (m/^\#/ || m/^\s*$/)
		{
			print VCF "$_\n";
		}
		else
		{
			my @data=split/\t/,$_;
			next if ($data[4]=~/\*/);
            print VCF join "\t",@data;
            print VCF "\n";
		}
	}
	close IN;
	close VCF;
        $cmd = "cd $outdir && $java -Xmx10240m -jar $snpEff $outprefix -v $vcffile_filter -c $outdir/../../Ref/snpEff.config -o gatk -ud $ud > $outfile" ;
	&run_or_die($cmd);

	(my $gatk_anno = $outfile) =~ s/.vcf$/.gatk.vcf/ ;
	$cmd = "perl $Bin/extract_oneEff_anno.pl -i $outfile -o $gatk_anno" ;
	&run_or_die($cmd);
	return ($gatk_anno);
}

sub statistic_annotation_result()
{
	my ($aanno_files, $anno_file, $mode) = @_ ;
	my $dir = "$outdir/anno_stat" ;
	mkdir $dir ;
	&link_to_dir($aanno_files, $dir);
	mkdir "$outdir/result" ;
	my $outfile = "$outdir/result/$mode.anno_stat.xls" ;
	my $cmd = "perl $Bin/anno_stat_v1.2.pl -id $dir -o $outfile -m $mode" ; #huangls 2015-8-10
	&run_or_die($cmd);
	$cmd = "python $Bin/plot_ann.py -i $outfile  -o $outdir/result   -p $mode.anno" ;
	&run_or_die($cmd);

	if ($mode =~ /Indel/i){
		my $out_stat_file = "$outdir/result/$outprefix.combine.$mode" ;
		$cmd = "perl $Bin/indel_stat_v1.1.pl -i $anno_file -o $out_stat_file -m $min_length -x $max_length" ;
		&run_or_die($cmd);
		my @disfile=();
		for my $i (@{$aanno_files}){
			(my $out = $i) =~ s/.vcf$/.indel.length/ ;
			my$bname=basename($i);
			my($sample)=$bname=~/\.(R\d+)\./;
			push @disfile,"$out.dis";
			$cmd = "perl $Bin/indel_stat_v1.3.pl -i $i -o $out -m $min_length -x $max_length -s $sample" ;
			&run_or_die($cmd);
		}
		
		&run_or_die("cat @disfile >$outdir/result/all.sample.indel.length.stat");
		&run_or_die("python $Bin/plot_indel_range.py -i $outdir/result/all.sample.indel.length.stat -o $outdir/result/ -n all.sample.indel.length.distribution");
	}

	return ;
}

sub link_to_dir()
{
	my ($aanno_files, $dir) = @_ ;
	for (my $i=0; $i<@{$aanno_files}; $i++){
		my $cmd = "ln -sf ".$aanno_files->[$i]." $dir/" ;
		&run_or_die($cmd);
	}

	return ;
}

