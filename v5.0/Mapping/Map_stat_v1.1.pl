#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my %opts;
GetOptions(\%opts,"id=s","rl=s","od=s","c=s");

if(!defined($opts{id}) || !defined($opts{c}) || !defined($opts{rl}) || !defined($opts{od}) || defined($opts{h}) )
{
	print <<"	Usage End.";
	
	Usage:

		-id          indir of bam files                           <dir>         must be given
		-rl          reference length file                        <file>        must be given
		-od          outdir of stat result                        <dir>         must be given
		-c	     config
		-h           Help document

	Usage End.

	exit;
}

##############################################
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

my $indir = $opts{id} ;
$indir = abs_path($indir);

my $reflen_file = $opts{rl};
$reflen_file = abs_path($reflen_file);

my $outdir = $opts{od} ;
mkdir($outdir,0755) unless -d $outdir;
$outdir = abs_path($outdir);

my $outprefix = $hash{Species} ;
my $maxproc = $hash{maxcpu};
my $queue = $hash{queue};
my $samtools = $hash{samtools};
my $picard_dir = $hash{picard_dir};
my $min_insert = 0;
my $max_insert = 800;
my $chrlistfile = $hash{chr_id};
my $sh_dir = "$outdir/sh_dir" ;
mkdir($sh_dir,0755) unless -d $sh_dir;

####################################################################################################
my @bamfiles = ();
&get_bamfiles(\@bamfiles, $indir);

# get depth file use samtools
my $adepth_files = &generate_depth_file(\@bamfiles);

# static map ratio
my ($map_stat_file) = &generate_map_ratio(\@bamfiles);

# static depth and coverage
my ($adepth_dis_fils, $depth_stat_file) = &stat_depth_and_coverage($adepth_files);

# link stat result to result dir
&link_stat_file_to_result_dir($map_stat_file, $depth_stat_file);

# draw reads depth distribution
&draw_reads_depth_distribution($adepth_dis_fils);

# draw coverage depth distribution
&draw_coverage_depth_distribution($adepth_files);

# draw insert size
&draw_insert_size_distribution(\@bamfiles);

###############Subs
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

## qsub
sub qsub()
{
	my ($shfile) = @_ ;
	my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc $maxproc --queue $queue --reqsub $shfile --independent" ;
	&run_or_die($cmd);
	return ;
}

## get bam files
sub get_bamfiles()
{
	my ($abamfiles, $dir) = @_ ;
	@{$abamfiles} = glob("$dir/*.bam");
	if (scalar (@{$abamfiles}) == 0){
		&show_log("No bam file found, please check") ;
		exit (1);
	}
	return ;
}

sub get_chrlist_file()
{
	my ($reflen_file, $chr_num) = @_ ;
	my $chrlistfile = "$outdir/chr.list" ;
	my $cmd = "perl $Bin/bin/generate_chrlist.pl -i $reflen_file -c $chr_num -o $chrlistfile" ;
	&run_or_die($cmd);

	return($chrlistfile);
}

sub generate_depth_file()
{
	my ($abamfiles) = @_ ;
	my $dir = "$outdir/depth" ;
	mkdir $dir ;
	my @depth_files = ();
	my $shfile = "$sh_dir/samtools.depth.sh" ;
	open (SH, ">$shfile") || die "Can't creat $shfile, $!\n" ;
	for my $bamfile (@{$abamfiles}){
		(my $basename = basename($bamfile)) =~ s/bam$/depth/ ;
		my $outfile = "$dir/$basename" ;
		print SH "$samtools depth $bamfile > $outfile\n" ;
		push @depth_files, $outfile ;
	}
	close(SH);
	&qsub($shfile);
	return(\@depth_files);
}

sub generate_map_ratio()
{
	my ($abamfiles) = @_ ;
	my $dir = "$outdir/map_stat" ;
	mkdir $dir ;
	my $shfile = "$sh_dir/samtools.stat.sh" ;
	open (SH, ">$shfile") || die "Can't creat $shfile, $!\n" ;
	for my $bamfile (@{$abamfiles}){
		(my $basename = basename($bamfile)) =~ s/bam$/stat/ ;
		my $outfile = "$dir/$basename" ;
		print SH "$samtools flagstat $bamfile > $outfile\n" ;
	}
	close(SH) ;
	&qsub($shfile) ;

	# merge stat result
	my $outfile = "$dir/$outprefix.map_stat.xls" ;
	my $cmd = "perl $Bin/bin/map_ratio_stat_v1.1.pl -id $dir -o $outfile" ;
	&run_or_die($cmd);
	return ($outfile) ;
}

sub stat_depth_and_coverage()
{
	my ($adepth_files) = @_ ;
	my @depth_dis_files = ();
	my $dir = "$outdir/depth_stat" ;
	mkdir $dir ;
	my $shfile = "$sh_dir/depth.stat.sh" ;
	open (SH, ">$shfile") || die "Can't open $shfile, $!\n" ;
	for my $depth_file (@{$adepth_files}){
		my $basename = basename($depth_file) ;
		my $outfile = "$dir/$basename" ;
		print SH "perl $Bin/bin/depth_coverage_stat_v1.1.pl -i $depth_file -l $reflen_file -o $outfile\n" ;
		push @depth_dis_files, "$outfile.dis.cut" ;
	}
	close(SH);
	&qsub($shfile);

	# merge
	my $outfile = "$dir/$outprefix.depth_stat.xls" ;
	my $cmd = "perl $Bin/bin/depth_coverage_result_merge.pl -id $dir -o $outfile" ;
	&run_or_die($cmd);

	return (\@depth_dis_files, $outfile);
}


sub link_stat_file_to_result_dir()
{
	my ($map_stat_file, $depth_stat_file) = @_ ;
	my $dir = "$outdir/result" ;
	mkdir $dir ;
	my $map_basename = basename($map_stat_file);
	my $map_result = "$dir/$map_basename" ;
	my $depth_basename = basename($depth_stat_file);
	my $depth_result = "$dir/$depth_basename" ;
	if (!-f $map_result){
		my $cmd = "ln -s $map_stat_file $map_result" ;
		&run_or_die($cmd) ;
	}
	if (!-f $depth_result){
		my $cmd = "ln -s $depth_stat_file $depth_result" ;
		&run_or_die($cmd) ;
	}

	return ;
}

sub draw_reads_depth_distribution()
{
	my ($adepth_dis_fils) = @_ ;
	my $dir = "$outdir/depth_distr_png" ;
	mkdir $dir ;
	my $shfile = "$sh_dir/depth_distr.sh" ;
	open (SH, ">$shfile") || die "Can't creat $shfile, $!\n" ;
	for my $depth_dis_file (@{$adepth_dis_fils}){
		my $basename = basename($depth_dis_file) ;
		my $outfile = "$dir/$basename.png" ;
		print SH "Rscript $Bin/bin/dualAxis_v1.r -i $depth_dis_file -o $outfile --x.col 1 --y1.col 2 --y2.col 3 --x.lab \"Sequencing depth\" --y1.lab \"Percent of base\" --y2.lab \"Percent of cumulative base\" --legend.xpos 0.7 --legend.ypos 0.9\n" ;
	}
	close(SH);
	&qsub($shfile);
	return ;
}

sub draw_coverage_depth_distribution()
{
	my ($adepth_files) = @_ ;
	my $dir = "$outdir/coverage_distr_png" ;
	mkdir $dir ;
	my $shfile = "$sh_dir/coverage_distr.sh" ;
	open (SH, ">$shfile") || die "Can't creat $shfile, $!\n" ;
	for my $depth_file (@{$adepth_files}){
		my $basename = basename($depth_file) ;
		print SH "perl $Bin/bin/plotReadDensity_v1.1.pl -i $depth_file -k $basename -o $dir -c $chrlistfile -l chr\n" ;
	}
	close(SH);
	&qsub($shfile);
	return ;
}

sub draw_insert_size_distribution()
{
	my ($abamfiles) = @_ ;
	my $dir = "$outdir/insert_distr_pdf" ;
	mkdir $dir ;
	my $shfile = "$sh_dir/insert_distr.sh" ;
	open (SH, ">$shfile") || die "Can't creat $shfile, $!\n" ;
	for my $bamfile (@{$abamfiles}){
		(my $basename = basename($bamfile)) =~ s/.bam$// ;
		my $outfile = "$dir/$basename.insert" ;
		my $cmd = "perl $Bin/bin/insert_size_distribution.pl -i $bamfile -o $outfile -m $min_insert -x $max_insert -c $config" ;
		print SH $cmd,"\n" ;
	}
	close(SH);
	&qsub($shfile);
	return ;
}


