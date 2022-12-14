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
#===============================================================
use newPerlBase;
my $Title="merge_stat";
my $version=$ver;
#my %config=%{readconf("$Bin/../../reseq_script.cfg")};
my %config=%{selectconf("$Bin/../../confdir/")}; 

#===========请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作

my %opts;
GetOptions(\%opts,"a=s","b=s","o=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{a}) || !defined($opts{b}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		
		Version: $ver

	Usage:

		-a           first stat file               <infile>           must be given

		-b           second stat file              <infile>           must be given

		-o           merged stat file              <outfile>          must be given

		-h           Help document

	Usage End.

	exit;
}

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
################

# get parameters
my $statfile1 = $opts{a} ;
my $statfile2 = $opts{b} ;
my $outfile = $opts{o} ;

# read stat files
my %hstat = ();
my @type1 = &reading_stat_file($statfile1, \%hstat);
my @type2 = &reading_stat_file($statfile2, \%hstat);
my @type = (@type1, @type2);

# output
&output_stat_result($outfile, \%hstat, \@type);

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

#&reading_stat_file($statfile1, \%hstat);
sub reading_stat_file()
{
	my ($statfile, $ahstat) = @_ ;
	open (IN, $statfile) || die "Can't open $statfile, $!\n" ;
	my $head = <IN> ;
	my ($sample, @type) = split /\s+/, $head ;
	while (<IN>){
		my ($sample, @values) = split ;
		for (my $i=0; $i<@values; $i++){
			$ahstat->{$sample}->{$type[$i]} = $values[$i] ;
		}
	}
	close(IN);

	return (@type);
}

#&output_stat_result($outfile, \%hstat, \@type);
sub output_stat_result()
{
	my ($outfile, $ahstat, $atype) = @_ ;
	open (OUT, ">$outfile") || die "Can't creat $outfile, $!\n" ;
	print OUT "#Sample\t", join("\t", @{$atype}), "\n" ;
	for my $sample (sort keys %{$ahstat}){
		print OUT $sample ;
		for my $type (@{$atype}){
			print OUT "\t", $ahstat->{$sample}->{$type} ;
		}
		print OUT "\n" ;
	}
	close(OUT) ;

	return ;
}


