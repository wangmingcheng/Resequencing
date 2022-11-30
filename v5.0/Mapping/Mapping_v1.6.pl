#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my %opts;
GetOptions(\%opts,"c=s","od=s","h" );

if(!defined($opts{c}) || !defined($opts{od}) || defined($opts{h}))
{
	print <<"	Usage End.";

	Usage:

		-c           config file	<infile>
		-od          outdir		<outdir>
		-h           Help document

	Usage End.

	exit;
}


## init parameter
my $configfile = $opts{c} ;
my $outdir = $opts{od} ;
mkdir($outdir,0755) unless -d $outdir;
$outdir = abs_path($outdir);
mkdir("$outdir/tmp",0755) unless -d "$outdir/tmp";
mkdir("$outdir/bwa_mid",0755) unless -d "$outdir/bwa_mid";
mkdir("$outdir/sh_dir",0755) unless -d "$outdir/sh_dir";
my $sh_dir = "$outdir/sh_dir";
my $tmp_dir = "$outdir/tmp";

my %hash;
my $data = {};
open CFG,$configfile or die $!;
while(<CFG>){
	chomp;
	next if /^$/ || /^\#/;
	my ($key,$value) = split /\s+/;
	$hash{$key} = $value;
	if($key =~ /^(R\d+)-1/){
		$data->{$1}->[0]=$value;
	}
	if($key =~ /^(R\d+)-2/){
		$data->{$1}->[1]=$value;
	}
}
close CFG;
my $ref = $hash{Ref2};
my $reflen = $hash{RefLen};
my $outprefix = $hash{Species};
my $bwa = $hash{bwa};
my $platform = $hash{platform};
my $java = $hash{java};
my $picard_dir = $hash{picard_dir};
my $chr_num = $hash{ChrNum};
my $chr_list = $hash{chr_id};
my $maxcpu = $hash{maxcpu};
&index_reference($ref,$reflen);
&mapping_reads();
&mapping_result_stat();


#################################################################################################
# index reference
sub index_reference(){
	my ($ref, $reflen) = @_ ;
	my $bwa_index_sh="$sh_dir/$outprefix.bwa.index.sh" ;
	my $ref_dir = dirname($ref);
	my $ref_ori = "$ref_dir/sequences.fa";
	open(SH, ">$bwa_index_sh") or die $! ;
	if ($reflen >= 1000000000){
		print SH "$bwa index -a bwtsw $ref\n" ;
		print SH "$bwa index -a bwtsw $ref_ori\n";
	}else{
		print SH "$bwa index -a is $ref\n" ;
		print SH "$bwa index -a is $ref_ori\n";
	}
	close(SH);
	system("sh $bwa_index_sh");
}

# mapping reads
sub mapping_reads(){
	my $map_mem_sh = "$sh_dir/$outprefix.map.mem.sh" ;
	open (SH, ">$map_mem_sh") or die $! ;
	foreach my $sample (sort keys %{$data}){
		my $fq1name = basename($data->{$sample}->[0]); 
		my $fq2name = basename($data->{$sample}->[1]); 
		print SH "$bwa mem $ref ", $data->{$sample}->[0], " ", $data->{$sample}->[1], " -t 6 -M -R \"\@RG\\tID:$sample\\tLB:$sample\\tPL:ILLUMINA\\tSM:$sample\" > $outdir/bwa_mid/$outprefix.$sample.sam\n" if $platform eq "Illumina";
		print SH "$bwa mem $ref ", $data->{$sample}->[0], " ", $data->{$sample}->[1], " -t 6 -M -R \"\@RG\\tID:$sample\\tLB:$sample\\tPL:BGISEQ-500\\tSM:$sample\\tPU:$sample\\tCN:BGI\" > $outdir/bwa_mid/$outprefix.$sample.sam\n" if $platform eq "BGI";
	}
	close(SH);
	system("parallel -j 20 <$map_mem_sh") ;

	my $map_bam_sh = "$sh_dir/$outprefix.map.bam.sh" ;
	my $out_bam_dir = "$outdir/result" ;
	mkdir($out_bam_dir,0755) unless -d $out_bam_dir;
	open (SH, ">$map_bam_sh") or die $! ;
	for my $sample (sort keys %{$data}){
		print SH "$java -Xmx30G -Djava.io.tmpdir=$tmp_dir -jar $picard_dir/SortSam.jar VALIDATION_STRINGENCY=LENIENT INPUT=$outdir/bwa_mid/$outprefix.$sample.sam OUTPUT=$out_bam_dir/$outprefix.$sample.sort.bam SORT_ORDER=coordinate TMP_DIR=$tmp_dir && rm $outdir/bwa_mid/$outprefix.$sample.sam\n" ;
	}
	close(SH);
	&qsub($map_bam_sh,$queue);
}

sub mapping_result_stat()
{
	my $indir = "$outdir/result" ;
	my $dir = "$outdir/Map_stat" ;
	my $ref_dir = dirname($ref) ;
	my $ref_len_file = "$ref_dir/ref.genome.fa.len";
	my $cmd = "perl $Bin/Map_stat_v1.1.pl -id $indir -od $dir -rl $ref_len_file -c  $configfile\n";
	&run_or_die($cmd);
	return ;
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
