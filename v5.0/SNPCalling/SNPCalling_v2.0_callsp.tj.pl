#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

my %opts;
GetOptions(\%opts,"od=s","id=s","c=s","h" );

if(!defined($opts{id}) || !defined($opts{od})|| !defined($opts{c}) || defined($opts{h})){
	print <<"	Usage End.";

	Forced parameters:
		-od          outdir                       <dir>        must be given
		-id          indir of bam file            <dir>        must be given
		-c	     config
		-h           Help document

	Usage End.

	exit(1);
}

## init parameters
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
my $ref = $hash{Ref2};
my $outprefix = $hash{Species};
my $sample_num = $hash{Sample};
my $ploidy = $hash{Ploidy};
my $queue = $hash{queue};
my $maxcpu = $hash{maxcpu};
my $maxnum = $hash{maxnum};
my $link = $hash{link};
my $chr_num = $hash{ChrNum};
my $picard_dir = $hash{picard_dir};
my $java = $hash{java};
my $gatk_dir = $hash{gatk_dir};
my $bcftools_dir = $hash{bcftools_dir}; 
my $samtools = $hash{samtools};
my $vcfutils = $hash{vcfutils};

my $thread_num = 4;

my $outdir = $opts{od} ;
mkdir($outdir,0755) unless -d $outdir;
my $indir = $opts{id} ;
my $tmp_dir = "$outdir/tmp" ;
mkdir($tmp_dir,0755) unless -d $tmp_dir;
##########################################################################
#准备bam文件，建立索引等
=cut
my @bamfiles = () ;
&prepare(\@bamfiles);

#bam文件去冗余
my @dup_bam_files=();
@dup_bam_files = &duplicate_marking(\@bamfiles);
#bam文件局部重比对
my @realn_bam_files=();
@realn_bam_files = &local_realignment(\@dup_bam_files);
=cut
#变异信号检测
my @realn_bam_files = glob ("$indir/*bam");
my ($vcf) = &variant_calling(\@realn_bam_files);
#过滤
my ($filter_snp_vcf, $filter_indel_vcf) = &filter_variants($vcf);
#位置还原
if($link eq 'link'){
	($filter_snp_vcf, $filter_indel_vcf) = &filter_pos_convert($filter_snp_vcf, $filter_indel_vcf);
}
#软连接
my ($final_snp_vcf, $final_indel_vcf) = &link_result($filter_snp_vcf, $filter_indel_vcf);
#转换成snplist文件
my $snp_list_file = &convert_vcf_to_snplist($final_snp_vcf) ;

if ($sample_num > 1){
	&convert_snpvcf_to_DEG_snpvcf($final_snp_vcf) ;
	&convert_indelvcf_to_DEG_indelvcf($final_indel_vcf) ;
}

##结果统计
&static_snp_result($snp_list_file);

###############Subs
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub prepare()
{
	my ($abamfiles) = @_ ;
	my $dir = "$outdir/prepare" ; 
	mkdir $dir ;
	&get_bamfiles($abamfiles, $indir);
	&index_fa();
	&index_bam_file($abamfiles, $dir);
	return ;
}

sub get_bamfiles()
{
	my ($abamfiles, $dir) = @_ ;
	@{$abamfiles} = glob("$dir/*.bam");
	if (scalar (@{$abamfiles}) == 0){
		print "No bam file found, please check\n" ;
		exit (1);
	}
	return ;
}


## index ref fa
sub index_fa()
{
	my $ref_dir = dirname($ref);
	my $ref_ori = "$ref_dir/sequences.fa";
	if (!-f "$ref.fai"){
		my $cmd = "$samtools faidx $ref" ;
		&run_or_die($cmd);
		$cmd = "$samtools faidx $ref_ori";
		&run_or_die($cmd);
	}
	(my $dict_file = $ref) =~ s/\.fa$|\.fasta$/.dict/ ;
	if (!-f "$dict_file"){
		my $cmd = "$java -Xmx50G -Djava.io.tmpdir=$tmp_dir -jar $picard_dir/CreateSequenceDictionary.jar REFERENCE=$ref OUTPUT=$dict_file" ;
		&run_or_die($cmd);
	}
	(my $dict_file_ori = $ref_ori) =~ s/\.fa$|\.fasta$/.dict/;
	if (!-f "$dict_file_ori"){
		my $cmd = "$java -Xmx50G -Djava.io.tmpdir=$tmp_dir -jar $picard_dir/CreateSequenceDictionary.jar REFERENCE=$ref_ori OUTPUT=$dict_file_ori" ;
		&run_or_die($cmd);
	}
	my $fa_length_file = "$ref.len" ;
	my $dir = dirname($ref);
	if (!-f $fa_length_file){
		my $cmd = "perl $Bin/ref_GC_len.pl -ref $ref -od $dir" ;
		&run_or_die($cmd);
	}
	return ;
}

## index for bam
sub index_bam_file()
{
	my ($abamfiles, $dir) = @_ ;
	my $shfile = "$dir/$outprefix.bam.idx.sh" ;
	open (SH, ">$shfile") || die "Can't creat $shfile, $!\n" ;
	for my $bamfile (@{$abamfiles}){
		(my $tmpfile = $bamfile) =~ s/.bam$/.bai/ ;
		if (!-f $tmpfile && !-f "$bamfile.bai"){
			print SH "$samtools index $bamfile\n" ;
		}
	}
	close(SH);
	&qsub($shfile,$queue, $maxcpu);
	return ;
}

sub get_sample_bam_files()
{
	my ($abamfiles, $ahsamples_bam) = @_ ;
	for my $bamfile (@{$abamfiles}){
		my $basename = basename($bamfile);
		my $sample_id = $basename ;
		if ($basename =~ /\.(R\d+)\.sort.bam$/ || $basename =~ /\.(R\d+)L\d+\.sort.bam$/){
			$sample_id = $1 ;
		}
		push @{$ahsamples_bam->{$sample_id}}, $bamfile ;
	}
	return ;
}

## Duplicate marking
sub duplicate_marking()
{
	my ($abamfiles) = @_ ;
	my $dir = "$outdir/duplicate_marking" ;
	mkdir $dir ;
	my @dupbamfiles = ();
	my $shfile = "$dir/$outprefix.dedup.sh" ;
	open (SH, ">$shfile") || die "Can't creat $shfile, $!\n" ;
	my %hsamples_bam = ();
	&get_sample_bam_files($abamfiles, \%hsamples_bam);
	for my $sample (sort keys %hsamples_bam){
		my $basename = "$outprefix.$sample.dedup.bam" ;
		my $outfile = "$dir/$basename" ;
		(my $metrics_file = "$dir/$basename") =~ s/bam$/metrics/ ;
		print SH "$java -Xmx20G -Djava.io.tmpdir=$tmp_dir -jar $picard_dir/MarkDuplicates.jar MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=512 VALIDATION_STRINGENCY=LENIENT" ;
		for my $bamfile (@{$hsamples_bam{$sample}}){
			print SH " INPUT=$bamfile" ;
		}
		print SH " OUTPUT=$outfile METRICS_FILE=$metrics_file\n" ;
		push @dupbamfiles, $outfile ;
	}
	close(SH) ;
	&qsub($shfile,$queue, $maxcpu);
	&index_bam_file(\@dupbamfiles, $dir);
	return (@dupbamfiles) ;
}

## Local realignment
sub local_realignment()
{
	my ($abamfiles) = @_ ;
	my @realnbamfiles = ();
	my $dir = "$outdir/local_realignment" ;
	mkdir "$dir" unless (-d $dir);
	my $shfile1 = "$dir/$outprefix.intervals.sh" ;
	open (SH1, ">$shfile1") || die "Can't creat $shfile1, $!\n" ;
	my $shfile2 = "$dir/$outprefix.realn.sh" ;
	open (SH2, ">$shfile2") || die "Can't creat $shfile2, $!\n" ;
	
	for my $bamfile (@{$abamfiles}){
		my $basename = basename($bamfile) ;
		(my $outfile = "$dir/$basename") =~ s/bam$/realn.bam/ ;
		my $intervals_file = "$bamfile.realn.intervals" ;
		print SH1 "$java -Xmx30G -Djava.io.tmpdir=$tmp_dir -jar $gatk_dir/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref -o $intervals_file -I $bamfile\n" ;
		print SH2 "$java -Xmx30G -Djava.io.tmpdir=$tmp_dir -jar $gatk_dir/GenomeAnalysisTK.jar -T IndelRealigner -R $ref -targetIntervals $intervals_file -o $outfile -I $bamfile \n" ;
		push @realnbamfiles, $outfile ;
	}
	
	close(SH1);
	close(SH2);
	&qsub($shfile1,$queue, $maxcpu);
	&qsub($shfile2,$queue, $maxcpu);
	return (@realnbamfiles) ;
}

sub variants_calling_by_gatk() 
{
	my ($abamfiles, $gatk_vcf, $dir) = @_ ;
	my $shfile = "$dir/$outprefix.0.gatk_var_call.sh";
	my $sample_num = @{$abamfiles};	
	my @list_files = &get_chr_list_file($dir,50);
	my $cmd = "";
	open (SH, ">$shfile") || die "Can't open $shfile, $!\n" ;
	my $combinesh = "$dir/$outprefix.1.CombineGVCFs.sh" ;
	open (SH1, ">$combinesh") || die "Can't open $combinesh, $!\n" ;
	my $genotypesh = "$dir/$outprefix.2.GenotypeGVCFs.sh" ;
	open (SH2, ">$genotypesh") || die "Can't open $genotypesh, $!\n" ;
	my $rmsh = "$dir/bye_gvcf.sh" ;
	open (SH3, ">$rmsh") || die "Can't open $rmsh, $!\n" ;
	my @vcf_files = () ;
	for (my $i=1; $i<=@list_files; $i++){
		my $list_file = $list_files[$i-1] ;
		my @comGVCF=();
		my @GVCF=();
		for (my $j=0;$j< @{$abamfiles};$j++){
		my$bamfile=$$abamfiles[$j];
		my $basename = "$outprefix.$i.$j.g.vcf";
		my $outfile = "$dir/$basename" ;
		print SH "$java  -XX:ParallelGCThreads=5  -Djava.io.tmpdir=$tmp_dir -Xmx20G -jar $gatk_dir/GenomeAnalysisTK.jar -T HaplotypeCaller -R $ref --indelSizeToEliminateInRefModel 50 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 --sample_ploidy $ploidy  " ;
		print SH "  -nct $thread_num -o $outfile  -L $list_file" ;
		push @GVCF, " --variant $outfile " ;
		print SH3 "rm $outfile\n";
		print SH3 "rm $outfile.idx\n";
		print SH " -I $bamfile \n";
		}
		print SH1 "$java -XX:ParallelGCThreads=5  -Djava.io.tmpdir=$tmp_dir -Xmx30G -jar $gatk_dir/GenomeAnalysisTK.jar -T CombineGVCFs -R $ref --disable_auto_index_creation_and_locking_when_reading_rods -o $dir/cohort.$i.g.vcf @GVCF \n";
		push @comGVCF," --variant $dir/cohort.$i.g.vcf ";
		print SH3 "rm $dir/cohort.$i.g.vcf\n";
		print SH3 "rm $dir/cohort.$i.g.vcf.idx\n";
		
		print SH2 "$java  -XX:ParallelGCThreads=5  -Djava.io.tmpdir=$tmp_dir -Xmx50G -jar $gatk_dir/GenomeAnalysisTK.jar -T GenotypeGVCFs -nt 1 -R $ref --disable_auto_index_creation_and_locking_when_reading_rods -o $dir/combine.$i.snp.indel.vcf @comGVCF \n";
		push @vcf_files, "$dir/combine.$i.snp.indel.vcf";
		print SH3 "rm $dir/combine.$i.snp.indel.vcf\n";
		print SH3 "rm $dir/combine.$i.snp.indel.vcf.idx\n";
	}
	close(SH);
	close(SH1);
	close(SH2);
	close(SH3);
	&qsub($shfile, $queue, $maxnum); 
	&qsub($combinesh, $queue, $maxcpu);
	&qsub($genotypesh, $queue, $maxcpu);
		
	my $combine_vcf_file = &combine_vcf_files_cat(\@vcf_files, $dir,$gatk_vcf);
	return;
}

sub combine_vcf_files_cat()
{
	my ($ahgvcf_files, $dir,$vcf) = @_ ;
	my $shfile = "$dir/$outprefix.combine.sample.vcf.sh" ;
	open (SH, ">$shfile") || die "Can't creat $shfile, $!\n" ;	
	print SH "$java -Xmx50G -Djava.io.tmpdir=$tmp_dir -cp $gatk_dir/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R $ref -out $vcf  -assumeSorted ";	
	for my $sample (@{$ahgvcf_files}){
		print SH " -V $sample " ;
	}
	print SH "\n" ;
	close(SH);
	&qsub($shfile, $queue, $maxcpu);
	return($vcf);
}


sub get_chr_list_file()
{
	my ($dir,$num) = @_ ;
	$num ||= 50;
	my $fa_len_file = "/share/nas1/hum/project/SLAF/BMK190630-S797-total-mangcao/callsnp_analysis/refined_ref/Miscanthus.fa.len"; #"$ref.len" ;
	my $outfile = "$dir/$outprefix.chr_allocation" ;
	my $cmd = "perl $Bin/distribute.pl -i $fa_len_file -o $outfile -n $num" ;
	if (!-f $outfile){
		&run_or_die($cmd) ;
	}
	my @list_files = ();
	@list_files = glob ("$outfile.*.list") ;
	my @orderList=();
	for (my $i=1; $i<=@list_files; $i++){
		push @orderList,"$outfile.".$i.".list";
	}
	return (@orderList);
}

sub variant_calling()
{
	my ($abamfiles) = @_ ;
	my $ref_dir = dirname($ref);
	my $dir = "$outdir/variant_calling" ;
	mkdir($dir,0755) unless -d $dir;
	my $gatk_vcf = "$dir/$outprefix.raw.vcf" ;
	&variants_calling_by_gatk($abamfiles, $gatk_vcf, $dir);
	return ($gatk_vcf);
}

sub extract_SNP_and_Indel()
{
	my ($variant_file, $dir) = @_ ;
	my $basename = basename($variant_file) ;
	(my $snp_vcf_file = "$dir/$basename") =~ s/vcf$/snp.vcf/ ;
	(my $indel_vcf_file = "$dir/$basename") =~ s/vcf$/indel.vcf/ ;
	my $cmd = "$java -Xmx50G -Djava.io.tmpdir=$tmp_dir -jar $gatk_dir/GenomeAnalysisTK.jar -T SelectVariants -R $ref -V $variant_file -selectType SNP -o $snp_vcf_file" ;
	if (!-f $snp_vcf_file){
		&run_or_die($cmd) ;
	}
	
	$cmd = "$java -Xmx50G -Djava.io.tmpdir=$tmp_dir -jar $gatk_dir/GenomeAnalysisTK.jar -T SelectVariants -R $ref -V $variant_file -selectType INDEL -o $indel_vcf_file" ;
	if (!-f $indel_vcf_file){
		&run_or_die($cmd) ;
	}
	return ($snp_vcf_file, $indel_vcf_file);
}

sub filter_variants()
{
	
	my ($vcf) = @_ ;
	my $dir = "$outdir/filter_variants" ;
	mkdir "$dir" ;
	my $snp_basename = basename($vcf);
	my $ref_dir = dirname($ref);
	(my $filter_snp_vcf = "$dir/$snp_basename") =~ s/vcf$/filter.vcf/ ;
	my $cmd="perl $vcfutils varFilter -w 5 -W 10 $vcf >$vcf.tmp && $java -Xmx50G -Djava.io.tmpdir=$tmp_dir -jar $gatk_dir/GenomeAnalysisTK.jar -T VariantFiltration -R $ref -V $vcf.tmp --filterExpression \"QUAL<30||QD < 2.0 || FS > 60.0 || MQ < 40.0  \"  --clusterWindowSize 5 --clusterSize 2 --filterName my_snp_filter -o $filter_snp_vcf.tmp && awk \'\$7 ==\"PASS\"|| \$0~/#/{print \$0}\' $filter_snp_vcf.tmp > $filter_snp_vcf &&rm $filter_snp_vcf.tmp && rm $vcf.tmp";
	if (!-f $filter_snp_vcf){
		&run_or_die($cmd);	
	}
	my($snp_vcf, $indel_vcf) = &extract_SNP_and_Indel($filter_snp_vcf,$dir);
	return ($snp_vcf, $indel_vcf);
}

sub result_pos_convert()
{
	my ($original_vcf, $convert_dir) = @_;
	my $ref_dir = dirname($ref);
	my $id_table = "$ref_dir/id_table.txt";
	my $cmd = "perl $Bin/convertpos_snplistv2.1.pl -v $original_vcf -pos $id_table -o $convert_dir";
	&run_or_die($cmd);
	my $converted_vcf = "$convert_dir/ConvertPos.vcf";
	return $converted_vcf;
}

sub filter_pos_convert()
{
	my ($snp_vcf, $indel_vcf) = @_;
	my $filter_snp_dir = dirname($snp_vcf);
	my $snp_convert_dir = "$filter_snp_dir/snp_posconvert";
	my $indel_convert_dir = "$filter_snp_dir/indel_posconvert";
	my $converted_snp_vcf = "$snp_convert_dir/ConvertPos.vcf";
	my $converted_indel_vcf = "$indel_convert_dir/ConvertPos.vcf";
	if (!-f $converted_snp_vcf){
		$converted_snp_vcf = &result_pos_convert($snp_vcf, $snp_convert_dir);
	}else{
		print "File exists: $converted_snp_vcf, position converting not do again\n";
	}
	if (!-f $converted_indel_vcf){
		$converted_indel_vcf = &result_pos_convert($indel_vcf, $indel_convert_dir);
	}else {
		print "File exists: $converted_indel_vcf, position converting not do again\n";
	}
	my $idx_snp = "$snp_vcf\.idx";
	if (-f $idx_snp) {
		my $cmd = "rm $idx_snp";
		&run_or_die($cmd);
	}
	my $idx_indel = "$indel_vcf\.idx";
	if (-f $idx_indel){
		my $cmd = "rm $idx_indel";
		&run_or_die($cmd);
	}
	if (-f "$converted_snp_vcf"){
		my $cmd = "mv $converted_snp_vcf $snp_vcf";
		&run_or_die($cmd);
	}
	if (-f "$converted_indel_vcf"){
		my $cmd = "mv $converted_indel_vcf $indel_vcf";
		&run_or_die($cmd);
	}
	return($snp_vcf, $indel_vcf);
}

sub link_result()
{
	my ($filter_snp_vcf, $filter_indel_vcf) = @_ ;
	my $dir = "$outdir/result" ;
	mkdir $dir ;

	my $snp_basename = basename($filter_snp_vcf) ;
	my $final_snp_vcf = "$dir/$snp_basename" ;
	my $indel_basename = basename($filter_indel_vcf) ;
	my $final_indel_vcf = "$dir/$indel_basename" ;

	my $cmd = "ln -s $filter_snp_vcf $filter_indel_vcf $dir/" ;
	if (!-f $final_snp_vcf || !-f $final_indel_vcf){
		&run_or_die($cmd);
	}else{
		my $txt = "Outfile is existed:\n[$final_snp_vcf] [$final_indel_vcf]\nSo this command not perform agian:\n[$cmd]" ;
		&show_log($txt);
	}
	
	return ($final_snp_vcf, $final_indel_vcf);
}

sub convert_vcf_to_snplist()
{
	my ($vcf_file) = @_ ;
	(my $outfile = $vcf_file) =~ s/vcf$/snp/ ;
	my $cmd = "perl $Bin/vcf_to_snplist_v1.2.pl -ref 1 -i $vcf_file -o $outfile" ;
	if (!-f $outfile){
		&run_or_die($cmd) ;
	}
	else{
		my $txt = "Outfile is existed:\n[$outfile]\nSo this command not perform agian:\n[$cmd]" ;
		&show_log($txt);
	}
	
	return($outfile) ;	
}

sub convert_snpvcf_to_DEG_snpvcf()
{
	my ($vcf_file) = @_;
	(my $DEG_vcf_file = $vcf_file)=~s/vcf$/DEG.vcf/;
	my $cmd = "perl $Bin/diffsnp_select_from_vcf.pl -i $vcf_file -o $DEG_vcf_file" ;
	if (!-f $DEG_vcf_file){
		&run_or_die($cmd);
	}
	else{
		my $txt = "Outfile is existed:\n[$DEG_vcf_file]\nSo this command not perform agian:\n[$cmd]" ;
		&show_log($txt);
	}
	my $DEG_snp=&convert_vcf_to_snplist($DEG_vcf_file);
	my $DEGdir=dirname($DEG_vcf_file);
	my $cmd1="perl $Bin/count_deg_snp_v1.1.pl -i $DEG_vcf_file -o $DEGdir";
	if (!-f "$DEGdir/DEG_snp.stat"){
		&run_or_die($cmd1);
	}
	else{
		my $txt = "Outfile is existed:\n[$DEGdir/DEG_snp.stat]\nSo this command not perform agian:\n[$cmd1]" ;
		&show_log($txt);
	}
}

sub convert_indelvcf_to_DEG_indelvcf()
{
	my ($vcf_file) = @_;
	(my $DEG_vcf_file = $vcf_file)=~s/vcf$/DEG.vcf/;
	my $cmd = "perl $Bin/diffsnp_select_from_vcf.pl -i $vcf_file -o $DEG_vcf_file" ;
	if (!-f $DEG_vcf_file){
		&run_or_die($cmd);
	}
	else{
		my $txt = "Outfile is existed:\n[$DEG_vcf_file]\nSo this command not perform agian:\n[$cmd]" ;
		&show_log($txt);
	}
	my $DEGdir=dirname($DEG_vcf_file);
	(my $DEG_vcf_list = $DEG_vcf_file)=~s/vcf$/list/;
	my $cmd1= "perl $Bin/simplify_snpEnffvcf.pl -snpEnff $DEG_vcf_file -o $DEG_vcf_list";
	if (!-f $DEG_vcf_list){
		&run_or_die($cmd1);
	}
	else{
		my $txt = "Outfile is existed:\n[$DEG_vcf_list]\nSo this command not perform agian:\n[$cmd1]" ;
		&show_log($txt);
	}
	my $cmd2= "perl $Bin/count_deg_indel_v1.1.pl -i $DEG_vcf_file -o $DEGdir";
	if (!-f "$DEGdir/DEG_indel.stat"){
		&run_or_die($cmd2);
	}
	else{
		my $txt = "Outfile is existed:\n[$DEGdir/DEG_indel.stat]\nSo this command not perform agian:\n[$cmd2]" ;
		&show_log($txt);
	}
}

sub static_snp_result()
{
	my ($snp_list_file) = @_ ;
	(my $outfile = $snp_list_file) =~ s/$/.stat/ ;
	my $cmd = "perl $Bin/snp_stat_v1.1.pl -i $snp_list_file -o $outfile" ;
	if (!-f $outfile){
		&run_or_die($cmd);
	}
	else{
		my $txt = "Outfile is existed:\n[$outfile]\nSo this command not perform agian:\n[$cmd]" ;
		&show_log($txt);
	}

	return ;
}

## show log
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
	print "$cmd\n";
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
	my ($shfile, $queue, $ass_maxproc) = @_ ;
	$ass_maxproc ||= $maxcpu ;
	$queue ||= $queue;
	my $dir=dirname($shfile);
	my $cmd = "cd $dir && /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc $ass_maxproc --queue $queue --resource vf=15G --reqsub $shfile --independent" ;
	&run_or_die($cmd);
	return ;
}

