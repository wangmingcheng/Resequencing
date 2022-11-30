#!/usr/bin/perl -w
# 
# Copyright (c) BIO_MARK 2014
# Writer:         Dengdj <dengdj@biomarker.com.cn>
# Program Date:   2014
# Modifier:       Dengdj <dengdj@biomarker.com.cn>
# Last Modified:  2014.
my $ver="1.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

#=============================================================
use newPerlBase;
my $Title="generate_chrlist";
my $version=$ver;
#my %config=%{readconf("$Bin/../../reseq_script.cfg")};
my %config=%{selectconf("$Bin/../../confdir/")};  
#==============================================================
my %opts;
GetOptions(\%opts,"i=s","o=s","c=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{i}) || !defined($opts{c}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		
		Version: $ver

	Usage:

		-i           chr length file                      <infile>     must be given
		-o           result file                          <outfile>    must be given
		-c           chr number for draw, 0 for all       <int>        must be given
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
my $chrlen_file = $opts{i} ;
my $outfile = $opts{o} ;
my $chr_num = $opts{c} ;

# reading chr length file
my %hlen = ();
&reading_chr_length_file($chrlen_file, \%hlen);

# sort chr and out put result
&output_result($outfile, $chr_num, \%hlen);


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

# &show_log("txt")
sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &sub_format_datetime(localtime($time));
	print "$Time:\t$txt\n" ;
	return ($time) ;
}

#&run_or_die($cmd);
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

## qsub
=pop
sub qsub()
{
	my ($shfile, $maxproc) = @_ ;
	if (`hostname` =~ /cluster/){
#		my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc $maxproc --reqsub $shfile --independent" ;
		my $cmd = "sh $config{qsub} --maxproc $maxproc --reqsub $shfile --independent" ;
		&run_or_die($cmd);
	}
	else{
#		my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc $maxproc --reqsub $shfile --independent" ;
		my $cmd = "sh $config{qsub} --maxproc $maxproc --reqsub $shfile --independent" ;
		&run_or_die($cmd);
	}

	return ;
}
=cut
#&reading_chr_length_file($chrlen_file, \%hlen);
sub reading_chr_length_file()
{
	my ($chrlen_file, $ahlen) = @_ ;
	open (IN, $chrlen_file) || die "Can't open $chrlen_file, $!\n" ;
	while (<IN>){
		chomp ;
		next if (m/^\#/ || m/^\s*$/);
		my ($chr, $len) = split ;
		$ahlen->{$chr} = $len ;
	}
	close(IN);

	return ;
}

#&output_result($outfile, $chr_num, \%hlen);
sub output_result()
{
	my ($outfile, $chr_num, $ahlen) = @_ ;
	my @select_chr = sort {$ahlen->{$b} <=> $ahlen->{$a}} keys %{$ahlen};
	if ($chr_num == 0){
		$chr_num = @select_chr ;
	}
	my @final_chr = @select_chr[0..($chr_num - 1)] ;
	open (OUT, ">$outfile") || die "Can't creat $outfile, $!\n" ;
	my ($chr1, $chr2) ;
	my @chr_list = sort {(($chr1)=$a=~/(\d+)/,$chr1) <=> (($chr2)=$b=~/(\d+)/,$chr2)} @final_chr ;
	for my $chr (@chr_list) {
		print OUT $chr, "\n" ;
	}
	close(OUT);

	return ;
}


