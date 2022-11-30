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

my $data=shift;
open (IN,$data) or die $!;
my %seq_hash;

my $sample;
my $genome;
while (<IN>){
    chomp;
    next if /^$/ || /^\#/;
    my @data=split /\s+/,$_;
    if ($data[0]=~/Genome/i){
        open (INDEX, ">$od/01.genome_index/bwa_index.sh");
        $genome=basename($data[2]);
        my $name=(split/\.fa|\.fasta|\.fa\.gz|\.fasta\.gz/,$genome)[0];
        print INDEX "cd $od/01.genome_index && `cp $data[2] .` && $config{bwa} index $genome && $config{samtools} faidx $genome &&  $config{samtools} dict -a $genome -o $name.dict\n";
    }
    if ($data[0]=~/Sample/i){
        $sample=$data[1];                                      
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

open (SNP,">$od/07.gatk_SelectVariants/gatk_SelectVariants.sh");
print SNP "$config{gatk} SelectVariants -R $od/01.genome_index/$genome -V $od/06.gatk_GenotypeGVCFs/combined_vcf.gz --select-type-to-include SNP -O $od/07.gatk_SelectVariants/SNP.vcf.gz\n";

open (VF,">$od/08.gatk_VariantFiltration/gatk_VariantFiltration.sh");
print VF "$config{gatk} VariantFiltration --variant $od/07.gatk_SelectVariants/SNP.vcf.gz --output $od/08.gatk_VariantFiltration/filtered_SNP.vcf.gz -R $od/01.genome_index/$genome -filter \"(vc.hasAttribute(\'QD\')&&QD<2.0)||MQ<40.0||(vc.hasAttribute(\'ReadPosRankSum\')&&ReadPosRankSum<-8.0)||(vc.hasAttribute(\'MQRankSum\')&&MQRankSum<-12.5\)\" --filter-name filtered\n";

`mkdir -p $od/work.sh` unless (-d "$od/work.sh");
`cd $od/work.sh/ && ln -s ../*/*sh .`;

#Phylogenetic tree

#system("$config{parallel} -j 10 <$od/00.data_clean/fastp.sh");
#system("sh $od/01.genome_index/bwa_index.sh");
#system("$config{parallel} -j 10 <$od/02.bwa_samtools/bwa_samtools.sh");
