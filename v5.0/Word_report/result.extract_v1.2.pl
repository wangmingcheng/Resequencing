#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $Time_Start = &sub_format_datetime(localtime($BEGIN_TIME));
print "Program Starts Time:$Time_Start\n";

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($data,$detail);
GetOptions(
	"help|?" =>\&USAGE,
	"detail:s"=>\$detail,
	"data:s"=>\$data,
	) or &USAGE;
&USAGE unless ($data || $detail);
my %cfg_hash;
my %sample_in_config;

open (IN,$detail) or die $!;
while (<IN>){
	chomp;
	next if /^$/ || /^\#/;
	my @data=split /\s+/,$_;
	$cfg_hash{$data[0]} = $data[1];
}
close(IN);

my $j=1;
open IN,$data or die $!;
while(<IN>){
	chomp;
	next if /^$/ || /^\#/;
	if(/Sample/){
		my ($bmk) = (split /\s+/)[1];
		$sample_in_config{$j-1} = $bmk;
		$j++;
	}
}
close IN;
	
###################################################################################################
my $indir=$cfg_hash{analysis_dir};
mkdir("$indir/Up_load",0755) unless -d "$indir/Up_load";
mkdir("$indir/BMK_Backup",0755) unless -d "$indir/BMK_Backup";
`cp $detail $data $indir/BMK_Backup`;
`cp $indir/Analysis/Ref/sequences.fa $indir/BMK_Backup`;
`cp $indir/Analysis/Ref/ref.genome.gff	$indir/BMK_Backup/genes.gff`;
`cp $indir/Analysis/Ref/sample.list	$indir/BMK_Backup/`;
`cp $indir/Work_sh/Reseq_pipeline.sh $indir/BMK_Backup/`;

open README,">$indir/Up_load/reseq_readme.txt" or die $!;
my $tmp = "
+------------------------------------------------------------------------------+
|                    重测序callsnp结果说明文档                                  |
+------------------------------------------------------------------------------+

目录结构及文件说明：
********************************************************************************";
print README $tmp;

mkdir("$indir/Up_load/Dataassess",0755) unless -d "$indir/Up_load/Dataassess";
open IN,"$indir/DataAssess/AllSample_GC_Q.stat" or die $!;
open (OUT,">$indir/Up_load/Dataassess/sample_data_assess.list") or die $!;
print OUT "BMK_ID\tClean_Reads\tClean_Base\tQ20(%)\tQ30(%)\tGC(%)\n";	
while(<IN>){
	next if /^$/ || /^\#/;
	next if /^SampleID/;
	my ($ID,$ReadSum,$BaseSum,$GC,$Q20,$Q30) = (split /\s+/)[0,1,2,3,5,7];
	print OUT "$ID\t$ReadSum\t$BaseSum\t$Q20\t$Q30\t$GC\n";
}
close IN;
close OUT;

if (-d "$indir/DataAssess/"){	
	mkdir("$indir/Up_load/Dataassess",0755) unless -d "$indir/Up_load/Dataassess";
	`cp $indir/DataAssess/AllSample_GC_Q.stat $indir/Up_load/Dataassess`;
	my @sample=glob"$indir/DataAssess/R*";
	for (my $i=0;$i<@sample;$i++){
		`cp $sample[$i]/*.png $indir/Up_load/Dataassess`;
		`cp $sample[$i]/*.pdf $indir/Up_load/Dataassess`;
	}
my $tmp = "
Dataassess/
├── AllSample_GC_Q.stat				#样品数据评估统计
├── R*.acgtn.pdf				#样品的各碱基比例分布图（pdf格式）
├── R*.acgtn.png				#样品的各碱基比例分布图（png格式）
├── R*.quality.pdf				#样品的质量分布图（pdf格式）
├── R*.quality.png				#样品的质量分布图（png格式）
└── sample_data_assess.list			#各个样品数据统计信息";
print README $tmp;
	}else{
		print "No such file or directory\n$indir/Rawdata/Dataassess\n";
}
##################################################################################################
if(-d "$indir/Analysis/Mapping"){
	mkdir("$indir/Up_load/Mapping",0755) unless -d "$indir/Up_load/Mapping";
	`cp $indir/Analysis/Mapping/Map_stat/coverage_distr_png/*.png* $indir/Up_load/Mapping/` if -d "$indir/Analysis/Mapping/Map_stat/coverage_distr_png/";
	`cp $indir/Analysis/Mapping/Map_stat/depth_distr_png/*.png $indir/Up_load/Mapping/` if -d "$indir/Analysis/Mapping/Map_stat/depth_distr_png/"; 
	`cp $indir/Analysis/Mapping/Map_stat/insert_distr_pdf/*.png* $indir/Up_load/Mapping/` if -d "$indir/Analysis/Mapping/Map_stat/insert_distr_pdf/";
	`cp $indir/Analysis/Mapping/Map_stat/result/*.xls $indir/Up_load/Mapping/` if -d "$indir/Analysis/Mapping/Map_stat/result/";
my $tmp = "
Mapping/
├── *.depth_stat.xls				#各样品覆盖深度和覆盖度比例统计
├── *.map_stat.xls				#各样品比对结果统计表
├── *.R*.sort.depth.cov.txt.png			#各样品测序数据在各染色体上的分布（pdf格式）
├── *.R*.sort.depth.cov.txt.png.pdf		#各样品测序数据在各染色体上的分布（png格式）
├── *.R*.sort.depth.dis.cut.png			#各样品测序深度和累计深度分布图（png格式）
├── *.R*.sort.insert.cut.png			#各样品插入片段分布图（pdf格式）
└── *.R*.sort.insert.cut.png.pdf		#各样品插入片段分布图（png格式）";
print README $tmp;		
	}else{
		print "No such file or directory\n$indir/Analysis/Mapping\n";
}
##############################################################################################
if (-d "$indir/Analysis/SNP/result/"){
	mkdir("$indir/Up_load/SNP",0755) unless -d "$indir/Up_load/SNP";
	mkdir("$indir/Up_load/INDEL",0755) unless -d "$indir/Up_load/INDEL";
	`cp $indir/Analysis/SNP/result/*.snp $indir/Up_load/SNP`;
	`cp $indir/Analysis/SNP/result/*snp.stat $indir/Up_load/SNP`;
	`cp $indir/Analysis/SNP/result/DEG_snp.stat $indir/Up_load/SNP`;
	`cp $indir/Analysis/SNP/result/DEG_indel.stat $indir/Up_load/INDEL`;
	`cp $indir/Analysis/SNP/result/*.snp.vcf $indir/Up_load/SNP`;
	`cp $indir/Analysis/SNP/result/*.indel.vcf $indir/Up_load/INDEL`;
}else{
	print "No such file or directory\n$indir/Up_load/SNP\n";
}
#####################################################################################
if (-d "$indir/Analysis/SV"){
	mkdir("$indir/Up_load/SV",0755) unless -d "$indir/Up_load/SV";
	my @maxSV=glob("$indir/Analysis/SV/*.max");
	for my $i(@maxSV){
		my $svname=basename($i);
		system("sed '1d' $i >$indir/Up_load/SV/$svname");
	}
	`cp $indir/Analysis/SV/*.info $indir/Up_load/SV`;
	`cp $indir/Analysis/SV/*.stat* $indir/Up_load/SV`;
	`cp $indir/Analysis/SV/*.list $indir/Up_load/SV`;
my $tmp = "
SV/
├── *.R*.max					#各样品的SV详细信息
├── *.sv.anno.info				#SV注释信息文件
├── *.sv.anno.stat				#SV注释统计文件
├── *.Diff.gene.list				#发生INS、DEL、INV的区域所在的基因
└── *.sv.stat.xls				#SV统计文件";
print README $tmp;				
	}else{
		print "No such file or directory\n$indir/Up_load/SV\n";
}
############################################################################################
if (-d "$indir/Analysis/CNV"){
	mkdir("$indir/Up_load/CNV",0755) unless -d "$indir/Up_load/CNV";
	my @CNVs = glob("$indir/Analysis/CNV/*/*CNVs");
	my @CNVann = glob("$indir/Analysis/CNV/*/*.xls");	
	for my $i (@CNVs){
		my $cnvname=basename($i);
		system("sed  '1i chr\tstart\tend\tpredicted copy number\ttype of alteration' $i >$indir/Up_load/CNV/$cnvname");
	}
	system("cp @CNVann $indir/Up_load/CNV");
my $tmp = "
CNV/
├── R*.CNV.gene.anno.xls			#CNV注释信息文件
└── *.R*.dedup.realn.bam_CNVs			#CNV结果文件";
print README $tmp;
	}else{
		print "No such file or directory\n$indir/Up_load/CNV\n";
}
##############################################################################################
if (-d "$indir/Analysis/Annotation/snp_anno/"){
	my @gatkvcf=glob "$indir/Analysis/Annotation/snp_anno/*.gatk.vcf";
	for (my $i=0;$i<@gatkvcf ;$i++){
		(my $outfile=basename($gatkvcf[$i]))=~s/\.vcf$/\.list/;
		my $cmd="perl $Bin/simplify_snpEnffvcf.pl -snpEnff $gatkvcf[$i] -o $indir/Up_load/SNP/$outfile ";
		&run_or_die($cmd) if (!-f "$indir/Up_load/SNP/$outfile");
		$cmd ="perl $Bin/vcf_to_snplist_v1.3.pl  -i $gatkvcf[$i] -o $indir/Up_load/SNP/snp.genotype.info.cloud -ref 1";
		&run_or_die($cmd) if (!-f "$indir/Up_load/SNP/snp.genotype.info.cloud");
	}
	`cp $indir/Analysis/Annotation/snp_anno/result/SNP.anno_stat.xls $indir/Up_load/SNP/`;
	`cp $indir/Analysis/Annotation/snp_anno/result/*.png $indir/Up_load/SNP/`;
	`cp $indir/Analysis/Annotation/snp_anno/result/*.pdf $indir/Up_load/SNP/`;
	`cp $indir/Analysis/Annotation/snp_anno/*.snp.anno.vcf $indir/Up_load/SNP/`;
	`cp $indir/Analysis/Annotation/snp_anno/*.anno.gatk.vcf $indir/Up_load/SNP/`;
	`cp $indir/Analysis/Annotation/snp_anno/*/*.list $indir/Up_load/SNP/`;
my $tmp = "
SNP/
├── DEG_snp.stat				#各样品间差异SNP数量统计
├── *.raw.filter.snp.anno.gatk.list		#过滤后所有样品的SNP注释文件		
├── *.raw.filter.snp.anno.vcf			#过滤后所有样品的SNP注释vcf文件
├── *.raw.filter.snp.anno.gatk.vcf		#过滤后所有样品的SNP注释vcf文件(只保留注释上的SNP变异位点)
├── *.raw.filter.snp.vcf			#基于硬过滤后的vcf文件（后续所有分析均基于该文件）
├── *.raw.filter.snp.DEG.snp			#各样品间有差异的SNP
├── *.raw.filter.snp.snp			#SNP文件
├── *.raw.filter.snp.snp.stat			#各样品各类型SNP统计
├── SNP.anno.R*.pie.pdf				#SNP注释结果统计饼图（pdf格式）
├── SNP.anno.R*.pie.png				#SNP注释结果统计饼图（png格式）
├── SNP.anno_stat.xls				#SNP注释结果统计
├── snp.genotype.info.cloud			#各样品SNP注释以及深度信息统计
├── SNP.mutation.distribution.png		#SNP各样品转换颠换图
├── SNP.quality.distribution.png		#SNP质量分布图（深度和距离累积图）
├── SNP.venn.pdf				#样品间SNP统计venn图（pdf格式），5个样品以内才会展示
└── SNP.venn.png				#样品间SNP统计venn图（png格式），5个样品以内才会展示";
print README $tmp;
	}else{
		print "No such file or directory\n$indir/Analysis/Annotation/snp_anno/\n";
}
if (-d "$indir/Analysis/Annotation/indel_anno/"){
	my @gatkvcf=glob "$indir/Analysis/Annotation/indel_anno/*.gatk.vcf";
	for (my $i=0;$i<@gatkvcf ;$i++){
		my $aa=basename($gatkvcf[$i]);
		(my $outfile = $aa) =~ s/\.vcf$/\.list/;
		my$cmd ="perl $Bin/vcf_to_indellist_v1.4.pl  -i $gatkvcf[$i] -o $indir/Up_load/INDEL/indel.genotype.info.cloud -ref 1";
		&run_or_die($cmd) if (!-f "$indir/Up_load/INDEL/indel.genotype.info.cloud");
		$cmd = "perl $Bin/simplify_snpEnffvcf.pl -snpEnff $gatkvcf[$i] -o $indir/Up_load/INDEL/$outfile \n" ;
		&run_or_die($cmd) if (!-f "$indir/Up_load/INDEL/$outfile");
	}
	`cp $indir/Analysis/Annotation/indel_anno/result/*.png $indir/Up_load/INDEL/`;
	`cp $indir/Analysis/Annotation/indel_anno/result/*.pdf $indir/Up_load/INDEL/`;
	`cp $indir/Analysis/Annotation/indel_anno/result/*.stat $indir/Up_load/INDEL/`;
	`cp $indir/Analysis/Annotation/indel_anno/result/Indel.anno_stat.xls $indir/Up_load/INDEL/`;
	`cp $indir/Analysis/Annotation/indel_anno/*.indel.anno.vcf $indir/Up_load/INDEL/`;
	`cp $indir/Analysis/Annotation/indel_anno/*.indel.anno.gatk.vcf $indir/Up_load/INDEL/`;
	`cp $indir/Analysis/Annotation/indel_anno/*/*.list $indir/Up_load/INDEL/`;
my $tmp = "
INDEL/
├── all.sample.indel.length.distribution.pdf				#所有样品INDEL长度分布图（pdf格式）
├── all.sample.indel.length.distribution.png				#所有样品INDEL长度分布图（png格式）
├── all.sample.indel.length.stat					#所有样品INDEL长度统计
├── DEG_indel.stat							#差异INDEL数量统计
├── Indel.anno.R*.pie.pdf						#注释结果统计饼图（pdf格式）
├── Indel.anno.R*.pie.png						#注释结果统计饼图（png格式）
├── Indel.anno_stat.xls							#INDEL注释结果统计
├── indel.genotype.info.cloud						#各样品INDEL注释以及深度信息统计
├── INDEL.venn.pdf							#样品间INDEL统计venn图（pdf格式），5个样品以内才会展示
├── INDEL.venn.png							#样品间INDEL统计venn图（pdf格式），5个样品以内才会展示
├── *.combine.Indel.dis.png						#INDEL长度统计分布图（总）
├── *.combine.Indel.stat						#INDEL在编码区以及非编码区的分布统计
├── *.raw.filter.indel.anno.gatk.list					#所有样品的INDEL注释文件
├── *.raw.filter.indel.anno.gatk.vcf					#各样品的INDEL注释vcf文件(只保留注释上的SNP变异位点)
├── *.raw.filter.indel.vcf						#基于硬过滤后的vcf文件（后续部分INDEL的注释均基于该文件）
└── *.raw.filter.indel.anno.vcf						#所有样品的INDEL注释vcf文件";
print README $tmp;		
	}
if (-d "$indir/Analysis/Circos"){
	mkdir("$indir/Up_load/Circos",0755) unless -d "$indir/Up_load/Circos";
	my @sample=glob "$indir/Analysis/Circos/R*";
	for (my $i=0;$i<@sample ;$i++){
		my $b_sam=basename($sample[$i]);
		`cp $sample[$i]/$b_sam.circos.png $indir/Up_load/Circos/`;
		`cp $sample[$i]/$b_sam.circos.svg $indir/Up_load/Circos/`;
	}
my $tmp = "
Circos/
├── R*.circos.png					#样品各类型变异在染色体的分布（png格式）
└── R*.circos.svg					#样品各类型变异在染色体的分布（svg格式）";
print README $tmp;	
}
if (-d "$indir/Diff_analysis"){	
	mkdir("$indir/Up_load/Diff_analysis",0755) unless -d "$indir/Up_load/Diff_analysis";
	`cp $indir/Diff_analysis/Diff_gene.stat $indir/Up_load/Diff_analysis/`;
	`cp $indir/Diff_analysis/Allsample_Diffgene.list $indir/Up_load/Diff_analysis/`;
	my @desample=glob "$indir/Diff_analysis/R*";
	for (my $i=0;$i<@desample ;$i++){
		my $bname=basename($desample[$i]);
		unless (-d "$indir/Up_load/Diff_analysis/$bname"){
			`mkdir $indir/Up_load/Diff_analysis/$bname`;
		}
		`cp $desample[$i]/Integrated_Function.annotation.xls $indir/Up_load/Diff_analysis/$bname`;
		`cp -r $desample[$i]/Diff_anno/Cog_Anno/ $indir/Up_load/Diff_analysis/$bname`;
		`cp -r $desample[$i]/Diff_anno/go_enrichment/ $indir/Up_load/Diff_analysis/$bname`;
		`cp -r $desample[$i]/Diff_anno/Graph/ $indir/Up_load/Diff_analysis/$bname`;
		`cp -r $desample[$i]/Diff_anno/pathway/ $indir/Up_load/Diff_analysis/$bname`;
	}
my $tmp = "
Diff_analysis/
├── Allsample_Diffgene.list					#所有样本变异基因数目统计
├── Diff_gene.stat						#不同类型变异基因注释数目
└── R*
    ├── Cog_Anno						#Cog注释结果
    │   ├── *.Cog.classfy.png
    │   └── *.Cog.classfy.stat
    ├── go_enrichment						#GO注释结果
    │   ├── *.GO.Biological.stat
    │   ├── *.GO.Biological.xls
    │   ├── *.GO.Cellular.stat
    │   ├── *.GO.Cellular.xls
    │   ├── *.GO_enrichment.stat.xls
    │   ├── *.GO.list
    │   ├── *.GO.list.txt
    │   ├── *.GO.Molecular.stat
    │   ├── *.GO.Molecular.xls
    │   └── *.GO.png
    ├── Graph							#topGO结果
    │   ├── *.KEGG.list
    │   ├── *.KEGG.Phase.png
    │   ├── *.topGO_BP.pdf
    │   ├── *.topGO_BP.png
    │   ├── *.topGO_BP.xls
    │   ├── *.topGO_CC.pdf
    │   ├── *.topGO_CC.png
    │   ├── *.topGO_CC.xls
    │   ├── *.topGO_MF.pdf
    │   ├── *.topGO_MF.png
    │   ├── *.topGO_MF.xls
    │   ├── topGO.list
    │   └── topGO.map
    ├── Integrated_Function.annotation.xls			#样本R*基因注释列表
    └── pathway							#KEGG结果
        ├── kegg_enrichment
        │   ├── *.Kegg.ko
        │   ├── *.KEGG.stat
        │   ├── *.KEGG.svg
        │   ├── *.KEGG.tree.stat
        │   └── *.KEGG.xls
        └── kegg_map
            ├── ko00010.html
            ├── ko00010.png
            ├── ……
            ├── ko04141.html
            └── ko04141.png";
print README $tmp;
}

`rm $indir/Up_load/SNP/all.sample.indel.length.stat` if(-f "$indir/Up_load/SNP/all.sample.indel.length.stat");
if (-d "$indir/Analysis/Annotation/snp_anno/Samples_diffsnp"){
	`cp -r $indir/Analysis/Annotation/snp_anno/Samples_diffsnp $indir/Up_load/SNP`;
}

$tmp = "
注：
1. 较小文本文件（小于50M）件建议老师用editplus文本文件编辑器打开或者在文件名后面加.xls后缀用excel打开方便查看。较大的文件（大于50M）建议用http://www.pilotedit.com/ 软件打开以避免电脑死机。
2. 结果文件，除图片外，其他均为文本文件，都可以用文本文件软件打开。
基本文件格式说明：
********************************************************************************
1. AllSample_GC_Q.stat 结果解释说明
SampleID				#样品名称
ReadSum					#样品reads数目
BaseSum					#样品碱基数目
GC(%)					#GC含量
N(%)					#N含量
Q20(%)					#Q20值
CycleQ20(%)				#CycleQ20值
Q30(%)					#Q30值

2. *.depth_stat.xls 结果解释说明
#BMK ID					#样品BMK编号
Ave_depth				#样品的平均覆盖深度
Cov_ratio_1X(%)				#深度为1X时对应的基因组覆盖比例
Cov_ratio_5X(%)				#深度为5X时对应的基因组覆盖比例
Cov_ratio_10X(%)			#深度为10X时对应的基因组覆盖比例

3. *.map_stat.xls 结果解释说明
#BMK ID					#样品BMK编号
Total_reads				#Clean Reads数，双端分别统计，即read1和read2记为2条reads
Mapped(%)				#定位到参考基因组的Clean Reads数占所有Clean Reads数的百分比
Properly mapped(%)			#双端测序序列均定位到参考基因组上且距离符合测序片段的长度分布

4. *.raw.filter.snp.anno.gatk.list 结果解释说明
#CHROM					#染色体编号
POS					#物理位置
REF					#参考基因组碱基类型
ALT					#突变碱基
QUAL					#质量值
R*_base					#样品R*的所有碱基类型
R*_alt_depth				#样品R*中的各碱基深度
R*_reads_depth				#样品R*中的reads在该位点的深度
SNPEFF_EFFECT				#样品R*结构注释结果

5. *.raw.filter.snp.anno.vcf 结果解释说明
#CHROM					#染色体编号
POS					#物理位置
ID					#无意义
REF					#参考基因组碱基类型
ALT					#突变碱基
QUAL					#质量值
FILTER					#GATK能使用其它的方法来进行过滤，过滤结果中通过则该值为”PASS”；若variant不可靠，则该项不为”PASS”或为”.”
INFO					#每一个参数的意义参考head部分以INFO起始的行；AC（Allele Count）：表示该Allele的数目；AF（Allele Frequence）：表示Allele的频率；AN（Allele Number）：表示Allele的总数目。对于1个二倍体样本来说：基因型 0/1 表示样本为杂合子，Allele数（AC）为1（双倍体的sample在该位点只有1个等位基因发生了突变），Allele的频率（AF）为0.5(二倍体的样本在该位点只有50%的等位基因发生了突变)，总的Allele（AN）为2； 基因型 1/1 则表示样本为纯合的，Allele数（AC）为2，Allele的频率（AF）为1，总的Allele（AN）为2；DP：reads覆盖度，其中有些reads被过滤掉；FS：使用Fisher’s精确检验来检测strand bias而得到的Phred格式的p值，该值越小越好。
FORMAT					#GT（Genotype）：表示样本基因型，数字之间用“/”隔开;AD（Allele Depth）：样本中每一种allele的reads覆盖度，中间以逗号（，）隔开，第一个数字表示ref基因型覆盖度，第二个数字表示第一个variant覆盖度，第三个数字表示第二个variant覆盖度（如果有第二种variant的话，以此类推）;DP（Depth）：表示样本在该位点的覆盖度;GQ（Genotype Quality）：基因型的质量值。Phred格式(Phred_scaled)的质量值，表示在该位点该基因型存在的可能性；该值越高，则Genotype的可能性越大；计算方法：Phred值 = -10 * log (1-p) ，其中p为基因型存在的概率。PL：指定的三种基因型的质量值(provieds the likelihoods of the given genotypes)。这三种指定的基因型为(0/0,0/1,1/1)，这三种基因型的概率总和为1。和之前不一样的是，该值越大，表明为该种基因型的可能性越小。 Phred值 = -10 * log (p)，其中p为基因型存在的概率。
R*					#样品R*的具体基因型与其它信息（0/0：无变异；0/1：杂合突变；1/1：纯合突变）

-------------------------------------------------------------------------------
copyright (c) BMK 2020
last update 2019.12.30";

print README $tmp;
close README;




#######################################################################################
my $Time_End   = sub_format_datetime(localtime(time()));
print STDOUT "Program Ends Time:$Time_End\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";

#######################################################################################

sub sub_format_datetime {
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
#########################################################################################
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
#########################################################################################
sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &sub_format_datetime(localtime($time));
	print "$Time:\t$txt\n" ;
	return ($time) ;
}
############################################################################################
sub USAGE {#
	my $usage=<<"USAGE";

Usage:

  -data		<file>  data file
  -detail	<file>  detail file
 
  -h         Help

USAGE
	print $usage;
	exit;
}

