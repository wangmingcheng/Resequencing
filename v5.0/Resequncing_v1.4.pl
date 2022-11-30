#! /usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);

my %opts;
GetOptions(\%opts,"c=s","od=s","start=s","end=s","assign=s","h" );

if( !defined($opts{c}) || !defined($opts{od}) || defined($opts{h}))
{
	print <<"	Usage End.";

	Usage:

		-c           config file                                must be given
		-od          outdir                                     must be given

		-start       start step                                 optional [1]
		-end         end step                                   optional [6]
		-assign      assign some steps to run                   optional [n1,n2,n3]
	
	step number:     
		             mapping                                    1
		             calling SNP and short Indel                2
		             calling SV                                 3
		             detect CNV                                 4
		             snp annotation                             5
	
	Other parameter:
		-h           Help document

	Usage End.

	exit;
}


###############################################################################
my $outdir = $opts{od} ;
mkdir($outdir,0755) unless -d $outdir;
$outdir = abs_path($outdir);
my $config = $opts{c};
$config = abs_path($config);
my $maxcpu;

################################prepare ref
&prepare_ref();

## run process
my @steps = &get_process($opts{assign}, $opts{start}, $opts{end});
for my $step (@steps){
	&run_process($step) ;
}

## get process
sub get_process(){
	my ($assign, $start, $end) = @_ ;
	my @process = ();
	if (defined $assign){
		@process = split /,/, $opts{assign} ;
		if (defined $start && defined $end){
			for (my $step=$start; $step<=$end; $step+=1){
				push @process, $step ;
			}
		}
		@process = sort {$a<=>$b} @process ;
	}
	else{
		$start ||= 1 ;
		$end ||= 5 ;
		@process = ($start..$end);
	}
	return(@process);
}

## run each process
sub run_process(){
	my ($num) = @_ ;
	SWITCH: {
		$num == 1 && do { &Mapping(); last SWITCH; };
		$num == 2 && do { &CallingSNP(); last SWITCH; };
		$num == 3 && do { &CallingSV(); last SWITCH; };
		$num == 4 && do { &DetectCNV(); last SWITCH; };
		$num == 5 && do { &Annotation(); last SWITCH; };
		print "step number: $num not between 1 to 5\n" ;
		exit (1);
	}

	return ;
}

######################################################################################################
sub prepare_ref()
{
	open CFG,$config or die $!;
	mkdir("$outdir/Ref",0755) unless -d "$outdir/Ref";
	open OUT,">$outdir/Ref/sample.list" or die $!;
	my ($link,$samtools);
	while(<CFG>){
		chomp;
		next if /^$/ || /^\#/;
		my ($key,$value) = split /\s+/;
		next if !defined $key;
		if($key eq "Ref1"){
			`cp $value $outdir/Ref` if(!-e "$outdir/Ref/sequences.fa");
		}
		if($key eq "Gff1"){
			`cp $value $outdir/Ref` if(!-e "$outdir/Ref/genes.gff");
			`cp $value $outdir/Ref/ref.genome.gff`  if(!-e "$outdir/Ref/ref.genome.gff");
		}
		if($key eq "chr_id"){
			`cp $value $outdir/Ref/chr_id.txt`;
		}
		if($key =~ /^(R\d+)-1/){
			print OUT "$1\n";
		}
		$link = $value if $key eq "link";
		$samtools = $value if $key eq "samtools";
		$maxcpu = $value if $key eq "maxcpu";
	}
	close CFG;
	close OUT;
	
	my $ref = "$outdir/Ref/sequences.fa";
	my $link_ref = "$outdir/Ref/ref.genome.fa";
	my $cmd = "perl $Bin/Ref/refine_ref_v2.1.pl -in $ref -od $outdir/Ref -of ref.genome.fa -chr $outdir/Ref/chr_id.txt -link $link\n" ;
	&run_or_die($cmd) if !-e $link_ref;
	
	my $refine_ref_len = "$outdir/Ref/ref.genome.fa.len";
	$cmd = "perl $Bin/Ref/ref_GC_len.pl -ref $link_ref -od $outdir/Ref\n" ;
	&run_or_die($cmd) if !-e $refine_ref_len;
	
	$cmd = "$samtools faidx $ref" ;
	&run_or_die($cmd) if !-e "$ref.fai";
	$cmd = "$samtools faidx $link_ref" ;
	&run_or_die($cmd) if !-e "$link_ref.fai";
}


###################################################################################################
sub Mapping()
{
	my $map_outdir = "$outdir/Mapping" ;
	mkdir($map_outdir,0755) unless -d $map_outdir;
	my $cmd = "perl $Bin/Mapping/Mapping_v1.6.pl -c $config -od $map_outdir\n" ;
	open OUT,">$outdir/../Work_sh/step2_Mapping.sh";
	print OUT $cmd;
	close OUT;
	&qsub("$outdir/../Work_sh/step2_Mapping.sh");
	my $num = `ls $outdir/Mapping/Map_stat/depth/* | wc -l`;chomp $num;
	if($num != 0){
		my $cmd_rd = "rm $outdir/Mapping/Map_stat/depth/*";
		&run_or_die($cmd_rd);
	}
}

###################################################################################################
sub CallingSNP()
{
	my $snp_outdir = "$outdir/SNP" ;
	mkdir($snp_outdir,0755) unless -d $snp_outdir;
	my $indir = "$outdir/Mapping/result" ;
	my $cmd = "perl $Bin/SNPCalling/SNPCalling_v2.0.pl -id $indir -od $snp_outdir -c $config\n" ;
	open OUT,">$outdir/../Work_sh/step3_CallingSNP.sh";
	print OUT $cmd;
	close OUT;
	&qsub("$outdir/../Work_sh/step3_CallingSNP.sh");
	my $bam = `ls $outdir/SNP/duplicate_marking/*.dedup.bam |wc -l`;chomp $bam;
	if($bam != 0){
		my $cmd_rd = "rm $outdir/SNP/duplicate_marking/*.dedup.bam && rm $outdir/SNP/duplicate_marking/*.dedup.bam.bai";
	        &run_or_die($cmd_rd);
	}
	$bam = `ls $outdir/Mapping/result/*.sort.bam |wc -l`;chomp $bam;
	if($bam != 0){
		my $cmd_rd = "rm $outdir/Mapping/result/*.sort.bam && rm $outdir/Mapping/result/*.sort.bam.bai";
	        &run_or_die($cmd_rd);
	}	
}

###################################################################################################
sub CallingSV()
{
	my $sv_outdir = "$outdir/SV" ;
	mkdir($sv_outdir,0755) unless -d $sv_outdir;
	my $indir = "$outdir/SNP/local_realignment" ;
	my $cmd = "perl $Bin/SVCalling/SVCalling_v1.2.pl -id $indir -od $sv_outdir -c $config\n" ;
	open OUT,">$outdir/../Work_sh/step4_CallingSV.sh";
	print OUT $cmd;
	close OUT;
	&qsub("$outdir/../Work_sh/step4_CallingSV.sh");
}

###################################################################################################
sub DetectCNV ()
{
	my$refGC="$outdir/Ref/ref.genome.fa.GC";
	if(!-f $refGC){die "can't find $refGC file ,please check\n";}
	open IN,"$refGC" or die "$!";
	my$gc=0;
	while(<IN>){
		chomp;
		my@tmp=split(/\t/);
		$gc=$tmp[1] if($tmp[0]=~/Total_GC/);
	}
	close(IN);
	my $GCrange=($gc/100-$gc*0.1/100)."-".($gc/100+$gc*0.1/100);
	my $indir = "$outdir/SNP/local_realignment";
	mkdir "$outdir/CNV" unless (-d "$outdir/CNV");
	my $cmd = "perl $Bin/CNV/cnvDetectByfreec.pl -c $config -id $indir -od $outdir/CNV -GCrange $GCrange" ;
	open OUT,">$outdir/../Work_sh/step5_DetectCNV.sh";
	print OUT $cmd;
	close OUT;
	&qsub("$outdir/../Work_sh/step5_DetectCNV.sh");

}


###################################################################################################
sub Annotation ()
{
	my $snp_vcf_file = glob("$outdir/SNP/result/*filter.snp.vcf") ;
	my $indel_vcf_file = glob("$outdir/SNP/result/*filter.indel.vcf");
	chomp($snp_vcf_file, $indel_vcf_file) ;
	my $sample_list_file = "$outdir/Ref/sample.list" ;
	mkdir "$outdir/Annotation/" unless -d "$outdir/Annotation/";
	mkdir "$outdir/Annotation/snp_anno" unless -d "$outdir/Annotation/snp_anno";
	my $dir = "$outdir/Annotation/snp_anno" ;
	&annotation_and_stat_snp($snp_vcf_file, $sample_list_file, $dir);
	
	mkdir "$outdir/Annotation/indel_anno" unless -d "$outdir/Annotation/indel_anno";	
	$dir = "$outdir/Annotation/indel_anno" ;
	&annotation_and_stat_indel($indel_vcf_file, $sample_list_file, $dir);
}

sub annotation_and_stat_snp()
{
	my ($vcffile, $sample_list_file, $dir) = @_ ;
	my $ref = "$outdir/Ref/sequences.fa" ;
	my $cmd = "perl $Bin/Annotation/SNPAnnotation_v1.2.pl -i $vcffile -od $dir -r $ref -mode SNP -c $config\n" ;
    	open OUT,">$outdir/../Work_sh/step6_1_Annotation_snp.sh";
	print OUT $cmd;
	close OUT;
	
	&qsub("$outdir/../Work_sh/step6_1_Annotation_snp.sh");
	return ;
}

sub annotation_and_stat_indel()
{
	my ($vcffile, $sample_list_file, $dir) = @_ ;
	my $ref = "$outdir/Ref/sequences.fa" ;
	my $cmd = "perl $Bin/Annotation/SNPAnnotation_v1.2.pl -i $vcffile -od $dir -r $ref -mode Indel -c $config\n" ;
    	open OUT,">$outdir/../Work_sh/step6_2_Annotation_indel.sh";
	print OUT $cmd;
	close OUT;
	
	&qsub("$outdir/../Work_sh/step6_2_Annotation_indel.sh");
	return ;
}

sub qsub()
{
	my ($shfile) = @_ ;
	my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc $maxcpu --reqsub $shfile --independent --queue $queue" ;
	&run_or_die($cmd);

	return ;
}

sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &sub_format_datetime(localtime($time));
	print "$Time:\t$txt\n" ;
	return ($time) ;
}

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
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
