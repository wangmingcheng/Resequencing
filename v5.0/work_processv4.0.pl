#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
#================================================================================
my @argv = @ARGV;
my ($data,$detail,$assign);
GetOptions(
	"help|?" =>\&USAGE,
	"data:s" =>\$data,
	"detail:s"=>\$detail,
	"assign:s"=>\$assign,
) or &USAGE;

&USAGE unless ($data || $detail  || $assign);

$data = abs_path($data);
$detail = abs_path($detail);

####################################################################################################
open DETAIL,$detail or die $!;
my %config;
while(<DETAIL>){
	chomp;
	next if /^$/ || /^\#/;
	next if /^\%/;
	my ($key,$value) = split /\s+/ ;
	$config{$key} = $value;
}
close DETAIL;
my $ploidy = $config{ploidy};
my $link = $config{link};
my $ud = $config{ud};
my $maxcpu = $config{maxcpu};
my $max = $config{max};
my $type = $config{type};
my $platform = $config{platform};
my $Project = $config{Project};
my $od = $config{analysis_dir};
mkdir($od,0755) unless -d $od;
mkdir("$od/Work_sh",0755) unless -d "$od/Work_sh";

####################################################################################################
$assign ||= "1,2,3,4,5,6,7,8,9,10";
my @step=split/,/,$assign;
print "assign step for run : @step\n";
my ($step2,$step3,$step4,$step5,$step6,$step7,$step8,$ana_step,@reseq_step);
my $pre_step ||= 100;
my $draw_step ||= 100;
my $diff_step ||= 100;
my $extr_step ||= 100;
my $word_step ||= 100;
my @all=sort {$a<=>$b} @step;
for (my $st=0;$st<@step ;$st++)
{
	if ($step[$st]==1)
	{
		$pre_step=shift (@all);
	}
	if ($step[$st]==2)
	{
		$step2=shift (@all)-1;
		push @reseq_step,$step2;
	}
	if ($step[$st]==3)
	{
		$step3=shift (@all)-1;
		push @reseq_step,$step3;
	}
	if ($step[$st]==4)
	{
		$step4=shift (@all)-1;
		push @reseq_step,$step4;
	}
	if ($step[$st]==5)
	{
		$step5=shift (@all)-1;
		push @reseq_step,$step5;
	}
	if ($step[$st]==6)
	{
		$step6=shift (@all)-1;
		push @reseq_step,$step6;
	}
	if ($step[$st]==7)
	{
		$draw_step=shift (@all);
	}
	if ($step[$st]==8)
	{
		$diff_step=shift (@all);
	}
	if ($step[$st]==9)
	{
		$extr_step=shift (@all);
	}
	if ($step[$st]==10)
	{
		$word_step=shift (@all);
	}
}

if (@reseq_step){
	$ana_step=join "\,",@reseq_step;
}

`echo "perl $0 @argv" >$od/Work_sh/Reseq_pipeline.sh`;

#############################################################################step1
if ($pre_step==1){
	my $cmd="perl $Bin/Data_Assess/preparation_analysis.pl -data $data -detail $detail -od $od\n";
	open OUT ,">$od/Work_sh/step1_Data_preparation.sh";
	print OUT $cmd;
	close OUT;
	&run_or_die($cmd);
}

if ($ana_step){
	mkdir("$od/Analysis",0755) unless -d "$od/Analysis";
	my $cmd="perl $Bin/Resequncing_v1.4.pl -c $od/main.conf -od $od/Analysis -assign $ana_step\n";
	open OUT,">$od/Work_sh/Reseq_subpipeline.sh";
	print OUT $cmd;
	close OUT;
	&run_or_die($cmd);
}

if ($draw_step==7)
{
	my $cmd="perl $Bin/Draw_Circos/process_Circos.pl -cfg $detail \n";
	open OUT,">$od/Work_sh/step7_Circos.sh";
	print OUT $cmd;
	close OUT;
	&run_or_die($cmd);
}

if ($diff_step==8)
{
	{
		my $cmd="perl $Bin/DiffAnalysis/snp_indel_stats.pl -cfg $detail\n";
		open OUT,">$od/Work_sh/step8_1_snp_indel_stats.sh";
		print OUT $cmd;
		close OUT;
		&run_or_die($cmd);
	}
	{
		my $cmd="perl $Bin/DiffAnalysis/gene.stat.pl -cfg $detail\n";
		open OUT,">$od/Work_sh/step8_2_gene.stat.sh";
		print OUT $cmd;
		close OUT;
		&run_or_die($cmd);
	}
	{
		my $cmd="perl $Bin/DiffAnalysis/Function_anno/extract_DEG_gene_v1.1.pl -cfg $detail\n";
		open OUT,">$od/Work_sh/step8_3_extract_DEG_gene.sh";
		print OUT $cmd;
		close OUT;
		&run_or_die($cmd);
	}
}

if ($extr_step==9)
{
	my $cmd="perl $Bin/Word_report/result.extract_v1.2.pl -detail $detail -data $data\n";
	open OUT,">$od/Work_sh/step9_extract_result.sh";
	print OUT $cmd;
	close OUT;
	&run_or_die($cmd);
}

if ($word_step==10)
{
   	my $cmd_html="python $Bin/Word_report/reseq_xml_report.py -detail $detail -data $data\n";
	open OUT,">$od/Work_sh/step10_web_report.sh";
	print OUT $cmd_html;
	close OUT;
	&run_or_die($cmd_html);
}

#######################################################################################
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
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

sub qsub()
{
	my ($shfile, $queue, $ass_maxproc) = @_ ;
	$queue ||= $queue ;
	$ass_maxproc ||= $maxcpu ;
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc $ass_maxproc --queue $queue --resource vf=15G --reqsub $shfile --independent" ;
	&run_or_die($cmd);
}
#########################################################################
sub USAGE {
	my $usage=<<"USAGE";

--------------------------------------------------------------------------------------------------------------------

Usage

	-data		data config				must be given
	-detail		detail config				must be given
	
	-assign		assign some steps to run		optional [default: 1,2,3,4,5,6,7,8,9,10]
		
	    step number:     
			Data_Preparation                           1
			mapping                                    2
			calling SNP and short Indel                3
			calling SV                                 4
			detect CNV                                 5
			snp annotation                             6
			Draw Circos                                7
			Different Analysis                         8
			Result extract                             9
			web  report                                10

--------------------------------------------------------------------------------------------------------------------

USAGE
	print $usage;
	exit;
}
