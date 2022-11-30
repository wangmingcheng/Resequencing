use strict;
use warnings;
use File::Basename qw(basename dirname);
use Cwd qw(abs_path getcwd);
use FindBin qw($Bin $Script);
use Getopt::Long;
use lib "$Bin/config/";
use newPerlBase;

my %config=%{readconf("$Bin/config/db_file.cfg")};

my $od=$ENV{PWD};
`mkdir -p $od/00.data_clean` unless (-d "$od/00.data_clean");
`mkdir -p $od/01.genome_index` unless (-d "$od/01.genome_index");
`mkdir -p $od/02.bwa_samtools` unless (-d "$od/02.bwa_samtools");
`mkdir -p $od/03.gatk_MarkDuplicates` unless (-d "$od/03.gatk_MarkDuplicates");
`mkdir -p $od/04.gatk_HaplotypeCaller` unless (-d "$od/04.gatk_HaplotypeCaller");
`mkdir -p $od/05.gatk_CombineGVCFs` unless (-d "$od/05.gatk_CombineGVCFs");
`mkdir -p $od/06.gatk_GenotypeGVCFs` unless (-d "$od/06.gatk_GenotypeGVCFs");
`mkdir -p $od/07.gatk_SelectVariants` unless (-d "$od/07.gatk_SelectVariants");
`mkdir -p $od/08.gatk_VariantFiltration` unless (-d "$od/08.gatk_VariantFiltration");
`mkdir -p $od/09.Exclude_Filtered_Variants` unless (-d "$od/09.Exclude_Filtered_Variants");
`mkdir -p $od/Phylogenetic_tree` unless (-d "$od/Phylogenetic_tree");
`mkdir -p $od/Population_structure` unless (-d "$od/Population_structure");
`mkdir -p $od/Principal_components_analysis/` unless (-d "$od/Principal_components_analysis/");

my $data=shift;
open (IN,$data) or die $!;
my %seq_hash;

my $sample;
my $genome;
my @sample_name;
while (<IN>){
    chomp;
    next if /^$/ || /^\#/;
    my @data=split /\s+/,$_;
    if ($data[0]=~/Genome/i){
        open (INDEX, ">$od/01.genome_index/bwa_index.sh");
        $genome=basename($data[2]);
        my $name=(split/\.fa|\.fasta|\.fa\.gz|\.fasta\.gz/,$genome)[0];
        print INDEX "cd $od/01.genome_index && `cp $data[2] .` && $config{bwa} index $genome && $config{samtools} faidx $genome &&  $config{samtools} dict -o $name.dict $genome\n";
    }
    if ($data[0]=~/Sample/i){
        $sample=$data[1];
        push @sample_name, $sample;
    }
    if ($data[0]=~/^fq1/){
        $seq_hash{$sample}{1}=$data[1];                            
    }
    if ($data[0]=~/^fq2/){
        $seq_hash{$sample}{2}=$data[1];                    
    }            
}
open (FASTP,">$od/00.data_clean/fastp.sh") or die $!;
open (ALN,">$od/02.bwa_samtools/bwa_samtools.sh");
open (MD,">$od/03.gatk_MarkDuplicates/gatk_MarkDuplicates.sh");
open (HC,">$od/04.gatk_HaplotypeCaller/gatk_HaplotypeCaller.sh");
my @gvcf;
foreach my $samp (sort keys %seq_hash){
    chomp $samp;
    `mkdir -p $od/00.data_clean/$samp` unless (-d "$od/00.data_clean/$samp");
    print FASTP "cd $od/00.data_clean/$samp && fastp -i $seq_hash{$samp}{1} -I $seq_hash{$samp}{2} -o $samp\_clean.1.fq.gz -O $samp\_clean.2.fq.gz --thread 7 --json $samp.json --html $samp.html\n";                
    print ALN "$config{bwa} mem -t 10 -M -R \'\@RG\\tID\:$samp\\tSM\:$samp\\tPL:illumina\' $od/01.genome_index/$genome $od/00.data_clean/$samp/$samp\_clean.1.fq.gz $od/00.data_clean/$samp/$samp\_clean.2.fq.gz | $config{samtools} view -@ 10 -bS | $config{samtools} sort -@ 10 -o $od/02.bwa_samtools/$samp.sorted.bam\n";
    print MD "$config{gatk} MarkDuplicates -I $od/02.bwa_samtools/$samp.sorted.bam -M $od/03.gatk_MarkDuplicates/$samp\_duplication.metrics -O $od/03.gatk_MarkDuplicates/$samp.MarkDups\_sorted.bam && $config{samtools} index -@ 10 $od/03.gatk_MarkDuplicates/$samp.MarkDups\_sorted.bam\n";
    print HC "$config{gatk} HaplotypeCaller --native-pair-hmm-threads 10 -R $od/01.genome_index/$genome -I $od/03.gatk_MarkDuplicates/$samp.MarkDups\_sorted.bam -O $od/04.gatk_HaplotypeCaller/$samp.gvcf.gz -ERC GVCF\n";
    push @gvcf, "-V $od/04.gatk_HaplotypeCaller/$samp.gvcf.gz";
}
open (CG,">$od/05.gatk_CombineGVCFs/gatk_CombineGVCFs.sh");
print CG "$config{gatk} CombineGVCFs -R $od/01.genome_index/$genome @gvcf -O $od/05.gatk_CombineGVCFs/combined_cohort.gvcf.gz\n";

open (GG,">$od/06.gatk_GenotypeGVCFs/gatk_GenotypeGVCFs.sh");

print GG "$config{gatk} GenotypeGVCFs -R $od/01.genome_index/$genome -V $od/05.gatk_CombineGVCFs/combined_cohort.gvcf.gz -O $od/06.gatk_GenotypeGVCFs/combined_vcf.gz\n";

open (SV,">$od/07.gatk_SelectVariants/gatk_SelectVariants.sh");
print SV "$config{gatk} SelectVariants -R $od/01.genome_index/$genome -V $od/06.gatk_GenotypeGVCFs/combined_vcf.gz --select-type-to-include SNP -O $od/07.gatk_SelectVariants/SNP.vcf.gz\n";
print SV "$config{gatk} SelectVariants -R $od/01.genome_index/$genome -V $od/06.gatk_GenotypeGVCFs/combined_vcf.gz --select-type-to-include INDEL -O $od/07.gatk_SelectVariants/INDEL.vcf.gz\n";

open (VF,">$od/08.gatk_VariantFiltration/gatk_VariantFiltration.sh");
print VF "$config{gatk} VariantFiltration --variant $od/07.gatk_SelectVariants/SNP.vcf.gz --output $od/08.gatk_VariantFiltration/filtered_SNP.vcf.gz -R $od/01.genome_index/$genome --filter-name \"QD_filter\" -filter \"QD < 2.0\" --filter-name \"FS_filter\" -filter \"FS > 60.0\" --filter-name \"MQ_filter\" -filter \"MQ < 40.0\" --filter-name \"SOR_filter\" --filter \"SOR > 4.0\"\n";
print VF "$config{gatk} VariantFiltration --variant $od/07.gatk_SelectVariants/INDEL.vcf.gz --output $od/08.gatk_VariantFiltration/filtered_INDEL.vcf.gz -R $od/01.genome_index/$genome --filter-name \"QD_filter\" -filter \"QD < 2.0\" -filter-name \"FS_filter\" -filter \"FS > 200.0\" --filter-name \"SOR_filter\" -filter \"SOR > 10.0\"\n"; 

#print VF "$config{gatk} VariantFiltration -V $od/07.gatk_SelectVariants/SNP.vcf.gz -O $od/08.gatk_VariantFiltration/filtered_SNP.vcf.gz -R $od/01.genome_index/$genome -filter  \"(vc.hasAttribute(\'FS\'&&FS>60.0)||vc.hasAttribute(\'SOR\'&& SOR>4.0)||vc.hasAttribute(\'QD\')&&QD<2.0)||MQ<40.0||(vc.hasAttribute(\'ReadPosRankSum\')&&ReadPosRankSum<-8.0)||(vc.hasAttribute(\'MQRankSum\')&&MQRankSum<-12.5\)\" --filter-name filtered\n";
#print VF "$config{gatk} VariantFiltration -V $od/07.gatk_SelectVariants/INDEL.vcf.gz -O $od/08.gatk_VariantFiltration/filtered_INDEL.vcf.gz -R $od/01.genome_index/$genome -filter \"(vc.hasAttribute(\'FS\'&&FS>200.0)||vc.hasAttribute(\'SOR\'&& SOR>10.0)||vc.hasAttribute(\'QD\')&&QD<2.0)\" --filter-name filtered\n";

#Exclude Filtered Variants
open (EFV,">$od/09.Exclude_Filtered_Variants/exclude_filtered_variants.sh");
print EFV "$config{gatk} SelectVariants --exclude-filtered -R $od/01.genome_index/$genome -V $od/08.gatk_VariantFiltration/filtered_SNP.vcf.gz -O $od/09.Exclude_Filtered_Variants/final_snp.vcf.gz\n";
print EFV "$config{gatk} SelectVariants --exclude-filtered -R $od/01.genome_index/$genome -V $od/08.gatk_VariantFiltration/filtered_INDEL.vcf.gz -O $od/09.Exclude_Filtered_Variants/final_indel.vcf.gz\n";

`mkdir -p $od/work.sh` unless (-d "$od/work.sh");
`cd $od/work.sh/ && ln -s ../*/*sh .`;

#Phylogenetic tree
open (PT,">$od/Phylogenetic_tree/Phylogenetic_tree.sh");
print PT "$config{VCF2Dis} -InPut $od/09.Exclude_Filtered_Variants/final_snp.vcf.gz -OutPut $od/Phylogenetic_tree/final_SNP_dis.mat\n"; 
print PT "$config{fastme} -i $od/Phylogenetic_tree/final_SNP_dis.mat -o $od/Phylogenetic_tree/final_SNP_dis.mat.tree\n";

#population structure
open (PS,">$od/Population_structure/Admixture.sh");
print PS "$config{vcftools} --gzvcf $od/09.Exclude_Filtered_Variants/final_snp.vcf.gz --plink --out $od/Population_structure/plink\n$config{plink} --noweb --file $od/Population_structure/plink --geno 0.05 --maf 0.05 --hwe 0.0001 --make-bed --out $od/Population_structure/admixture\n";
for my $k (2..15){ 
    print PS "$config{admixture} --cv -j10 $od/Population_structure/admixture.bed $k > $od/Population_structure/log$k.out\n"; 
}

print PS "cd $od/Population_structure && awk \'/CV/ \{print \$3,\$4\}\' *out | cut -c 4,7-20 > admixture.cv.error\n";

print PS "cd $od/Population_structure && Rscript $Bin/bin/plotADMIXTURE.r -p admixture -i admixture.nosex -k 10 -l ",join(",",@sample_name),"\n";

#principal components analysis

open (PC,">$od/Principal_components_analysis/PCA.sh");
print PC "cd $od/Principal_components_analysis && $config{plink} --vcf $od/09.Exclude_Filtered_Variants/final_snp.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids \@:\# --indep-pairwise 100 10 0.1 --out plink\n";
print PC "cd $od/Principal_components_analysis && $config{plink} --vcf $od/09.Exclude_Filtered_Variants/final_snp.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids \@:\# --extract plink.prune.in --make-bed --pca --out pca\n";
print PC "cd $od/Principal_components_analysis && Rscript $Bin/bin/pca.r\n";

#PCA methods 2
print PC "\#methods 2\n\#cd $od/Principal_components_analysis && $config{plink} --vcf $od/09.Exclude_Filtered_Variants/final_snp.vcf.gz --allow-extra-chr --out plink\n";

print PC "\#cd $od/Principal_components_analysis && $config{plink} --bfile plink --pca 5 --out pca --allow-extra-chr\n";
print PC "\#cd $od/Principal_components_analysis && Rscript $Bin/bin/pca.r\n";

#system("$config{parallel} -j 10 <$od/00.data_clean/fastp.sh");
#system("sh $od/01.genome_index/bwa_index.sh");
#system("$config{parallel} -j 10 <$od/02.bwa_samtools/bwa_samtools.sh");
