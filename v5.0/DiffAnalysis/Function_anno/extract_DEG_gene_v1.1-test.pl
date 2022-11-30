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
my $cfg;
GetOptions(
	"help|?" =>\&USAGE,
	"cfg:s"=>\$cfg,
) or &USAGE;

&USAGE unless ($cfg);

my %hash;
open(CFG,$cfg) or die $!;
while (<CFG>)
{
	chomp;
	next if ($_=~/\#/);
	next if ($_=~/^\s*$/);
	my @data=split/\s+/,$_;
	$hash{$data[0]} = $data[1];
}
close(CFG);

my %gene_hash;
my %anno_hash;
my $indir="$hash{analysis_dir}";
unless (-d "$indir/Diff_analysis"){
		`mkdir $indir/Diff_analysis`;
}
if (-d "$indir/Analysis/Annotation/") {
	my @snpsample=glob "$indir/Analysis/Annotation/snp_anno/R*";
	for (my $i=0;$i<@snpsample ;$i++){
        next if !-d $snpsample[$i];
		my $bname=basename($snpsample[$i]);
		next if ($bname !~ /^R\d+/);
		unless (-d "$indir/Diff_analysis/$bname"){
			`mkdir $indir/Diff_analysis/$bname`;
		}
		`cp $snpsample[$i]/*.list $indir/Diff_analysis/$bname/`;
	}
	my @indelsample=glob "$indir/Analysis/Annotation/indel_anno/R*";
	for (my $j=0;$j<@indelsample ;$j++){
        next if !-d $indelsample[$j];
		my $bname=basename($indelsample[$j]);
		next if ($bname !~ /^R\d+/);
		unless (-d "$indir/Diff_analysis/$bname")
		{
			`mkdir $indir/Diff_analysis/$bname`;
		}
		`cp $indelsample[$j]/*.list $indir/Diff_analysis/$bname/`;
	}
	if(-d "$indir/Analysis/SV/"){
		my @svsample=glob "$indir/Analysis/SV/*.list";
		for (my $n=0;$n<@svsample ;$n++){
			my $bname=basename($svsample[$n]);
			my $sam=(split/\./,$bname)[0];
			unless (-d "$indir/Diff_analysis/$sam"){
				`mkdir $indir/Diff_analysis/$sam`;
			}
			`cp $svsample[$n] $indir/Diff_analysis/$sam`;
		}
	}
}
my $annodb="$hash{anno}";
open(ANNO,"$annodb/Integrated_Function.annotation.xls") || die "Can't open $annodb/Integrated_Function.annotation.xls, $!\n" ;
my $head = "" ;
while (<ANNO>){
	chomp;
	if ($_=~/^\#/){
		$head .= $_."\n" ;
		next ;
	}
	next if ($_=~/^\s*$/);
	my @data=split/\t/,$_;
    print "$data[0]\n";
	$anno_hash{$data[0]}=$_;
}
close(ANNO);

print "$head\n";
my @allsample=glob "$indir/Diff_analysis/R*";
for (my $i=0;$i<@allsample ;$i++){
	my $bname=basename($allsample[$i]);
	my @gene=glob "$allsample[$i]/*";
	for (my $j=0;$j<@gene ;$j++){
		open(GENE,$gene[$j]) or die $!;
		while (<GENE>){
			chomp;
			next if ($_=~/\#/);
			next if ($_=~/^\s*$/);
			my $gene=(split/\s+/,$_)[0];
			$gene_hash{$bname}{$gene}=$gene;
           # print "$bname\t$gene\n";
		}
	}
}
open(ALLDE,">$indir/Diff_analysis/Allsample_Diffgene.list");
print ALLDE "Sample\tGene_Number\n";
foreach my $sam (sort keys %gene_hash){
	my $all_deg_num=0;
	my $all_deg_anno_num=0;
	open(AN,">$indir/Diff_analysis/$sam/Integrated_Function.annotation.xls");
	open(DE,">$indir/Diff_analysis/$sam/$sam.vs.ref.gene.list");
	print DE "#GeneID\t\t\tMarker\n";
	#print AN "#GeneID\tCOG_class\tCOG_class_annotation\tGO_annotation\tKEGG_annotation\tSwissprot_annotation\tnr_annotation\n";
	print AN "$head" ;
	foreach my $id (sort keys %{$gene_hash{$sam}}){
		print DE "$id\t.\t.\tup\n";
		$all_deg_num++;
        #print "$id\n";
		if (exists $anno_hash{$id}){
			$all_deg_anno_num++;
			print AN "$anno_hash{$id}\n";
		}
	}
	close(AN);
	close(DE);
	print ALLDE "$sam\t";
	print ALLDE "$all_deg_anno_num\n";
}
close(ALLDE);
unless (-d "$indir/Diff_analysis/work_sh"){
	`mkdir $indir/Diff_analysis/work_sh`;
}
open(SH,">$indir/Diff_analysis/work_sh/deg_anno.sh");
for (my $i=0;$i<@allsample ;$i++){
	my $bname=basename($allsample[$i]);
	my $deggenelist="$allsample[$i]/$bname.vs.ref.gene.list";
	my $datafile="$hash{anno}/All_Database_annotation.xls";
	print SH "perl $Bin/enrichment_analysis.pl -queue $hash{queue} -anno $datafile -deg $deggenelist -key $hash{Project} -outdir $allsample[$i]/Diff_anno\n";
}
close(SH);
&qsub("$indir/Diff_analysis/work_sh/deg_anno.sh");


#######################################################################################
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
	my ($shfile) = @_ ;
	my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc $hash{maxcpu} --queue $hash{queue} --resource vf=15G --reqsub $shfile --independent" ;
	&run_or_die($cmd);
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
  -h         Help
  
USAGE
	print $usage;
	exit;
}


