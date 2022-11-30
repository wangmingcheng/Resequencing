#!/usr/bin/perl -w
# 
# Copyright (c) BMK 2014
# Writer:         lium <lium@biomarker.com.cn>
# Program Date:   2014/4/10.
# Modifier:       lium <lium@biomarker.com.cn>
# Last Modified:  2014/4/10.
#
# ------------------------------------------------------------------
# add package dir
# ------------------------------------------------------------------
unshift(@INC, "/share/nas2/genome/bmksoft/tool/bmkPerlBase/current/");


# ------------------------------------------------------------------
# setting version number
# ------------------------------------------------------------------
my $version="1.0";


# ------------------------------------------------------------------
# load package
# ------------------------------------------------------------------
use strict;
use Cwd;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require bmkPerlBase;



# ------------------------------------------------------------------
# new bmkPerlBase
# ------------------------------------------------------------------
our $bmk = bmkPerlBase->new();


# ------------------------------------------------------------------
# main pipeline
# ------------------------------------------------------------------
# Time start
my $BEGIN=time();
$bmk->timeLog(-info=>"Start Time");

# pipeline
main_pipeline();

# Time end
$bmk->timeLog(-info=>"End Time");
$bmk->totalRunTime(-start_time=>"$BEGIN");



# ------------------------------------------------------------------
# the define of main_pipeline
# it contains main step in this program
# ------------------------------------------------------------------
sub main_pipeline{
	# init variable
	my %opts;
	my %RunArgv = ();

	# get option and check
	get_options(\%opts);

	# check input file
	para_check(\%opts, \%RunArgv);
	$bmk->timeLog(-info=>"para_check is done");

	# check exp and create awk exp
	create_qsub_shell(\%opts, \%RunArgv);
	$bmk->timeLog(-info=>"create_qsub_shell is done");

}


# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
sub get_options {
	my ($opts) = @_;
	GetOptions($opts, "anno=s", "deg=s", "outdir=s", "key=s", "queue=s","h" ) or USAGE();

	# check
	USAGE() if(!defined($opts->{"anno"}) || !defined($opts->{"deg"}) || !defined($opts->{"outdir"}) 
		|| !defined($opts->{"key"}) || defined($opts->{"h"}));
}

# ------------------------------------------------------------------
# the regular expression create and check 
# ------------------------------------------------------------------
sub create_qsub_shell {
	# get parameter
	my ($opts, $run_argv) = @_;

	# open file
	my $outdir = $run_argv->{outdir};
	open(OUT, ">$outdir/enrich.sh") 
		or $bmk->logAndDie(-info=>"open or create file ($outdir/enrich.sh) is failed");

	# print enrich cmd
	my $enrich = $run_argv->{anno};
	my $deg = $run_argv->{deg};
	my $key = $opts->{key};
	print OUT "perl $Bin/bin/Enrichment/select_DEG_Anno.pl -i $enrich -deg $deg -o $outdir/$key.annotation.xls\n";
	print OUT "perl $Bin/bin/Enrichment/draw_DEG_GO_graph.pl -i $enrich -deg $deg -k $key -od $outdir && perl $Bin/bin/Enrichment/draw_top_GO_graph.pl -i $enrich -deg $deg -k $key -od $outdir\n";
	#print OUT "perl $Bin/bin/Enrichment/draw_top_GO_graph.pl -i $enrich -deg $deg -k $key -od $outdir\n";
	print OUT "perl $Bin/bin/Enrichment/draw_DEG_KEGG_graph.pl -i $enrich -deg $deg -k $key -od $outdir\n";
	print OUT "perl $Bin/bin/Enrichment/draw_DEG_COG_graph.pl -i $enrich -deg $deg -k $key -od $outdir\n";
	print OUT "perl $Bin/bin/Enrichment/draw_DEG_KOG_graph.pl -i $enrich -deg $deg -k $key -od $outdir\n";
	close(OUT); 

	# calling qsub
	system"sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc 100 --queue $opts->{queue} --resource vf=15G --reqsub $outdir/enrich.sh --independent";
}




# ------------------------------------------------------------------
# check parameter 
# ------------------------------------------------------------------
sub para_check {
	# get parameter
	my ($opts, $run_argv) = @_;

	# check input file
	if(!-f $opts->{"anno"}){
		$bmk->logAndDie(-info=>"Error: file \"$opts->{anno}\" not exist");
	}
	if(!-f $opts->{"deg"}){
		$bmk->logAndDie(-info=>"Error: file \"$opts->{deg}\" not exist");
	}

	# get absolute path
	# anno
	$run_argv->{"anno"} = $bmk->absolutePath(-in=>$opts->{"anno"});
	# deg
	$run_argv->{"deg"} = $bmk->absolutePath(-in=>$opts->{"deg"});

	# outdir
	$bmk->mkdirOrDie(-dir=>$opts->{"outdir"});
	$run_argv->{"outdir"} = $bmk->absolutePath(-in=>$opts->{"outdir"});
}



	



# ------------------------------------------------------------------
# print usage information and exit
# ------------------------------------------------------------------
sub USAGE {
	my $usage=<<"USAGE";
Version: 	Version$version
Writer: 	lium <lium\@biomarker.com.cn>
Program Date: 	2014/4/10.
Modifier: 	lium <lium\@biomarker.com.cn>
Last Modified:  2014/4/10.

Description: 	do enrichment analysis for part of gene
Usage:
  Example: 	1) perl enrichment_analysis.pl -anno All_Database_annotation.xls -deg T1_vs_T2.DEG_final.xls \\
			-key T1_vs_T2 -outdir ./outdir
  Options:
  forced parameters:
  -anno 	<str> 	the file name of all gene anno
  -deg 		<str> 	the file name of part gene
  -key 		<int> 	the prefix of output
  -outdir 	<int> 	the dir for output
  optional parameters:
  -m		max prosess for qsub,default 100
  -h 		<none> 	Help

USAGE
	print $usage;
	exit(1);
}






