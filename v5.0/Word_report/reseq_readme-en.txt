
+------------------------------------------------------------------------------+
|             Genome Resequencing SNP calling Results Content                                  |
+------------------------------------------------------------------------------+

Result list and description：
********************************************************************************
Dataassess/
├── AllSample_GC_Q.stat				#Statistic table of data quality
├── R*.acgtn.pdf				#Distribution map of base ratio(pdf format)
├── R*.acgtn.png				#Distribution map of base ratio(png format)
├── R*.quality.pdf				#Distribution map of sample data quality(pdf format)
├── R*.quality.png				#Distribution map of sample data quality(png format)
└── sample_data_assess.list			#Statistic table of each sample data
Mapping/
├── *.depth_stat.xls				#Statistic table of sequencing depth and sequencing coverage of each sample
├── *.map_stat.xls				#Statistic table of alignment results of each sample
├── *.R*.sort.depth.cov.txt.png			#Distribution of sequencing data of each sample on each chromosome(pdf format)
├── *.R*.sort.depth.cov.txt.png.pdf		#Distribution of sequencing data of each sample on each chromosome(png format)
├── *.R*.sort.depth.dis.cut.png			#Distribution map of sequencing depth and cumulative depth(png format)
├── *.R*.sort.insert.cut.png			#Distribution map of insertion fragment(pdf format)
└── *.R*.sort.insert.cut.png.pdf		#Distribution map of insertion fragment(png format)
SV/
├── *.R*.max					#Detailed information of SV
├── *.sv.anno.info				#Annotation information of SV
├── *.sv.anno.stat				#Aoonation statistics of SV
├── *.Diff.gene.list				#Gene list of INS,DEL and INV
└── *.sv.stat.xls				#Statistic table of SV
CNV/
├── R*.CNV.gene.anno.xls			#Annotation information of CNV
└── *.R*.dedup.realn.bam_CNVs			#CNV results
SNP/
├── DEG_snp.stat				#Differential SNP statistics
├── *.raw.filter.snp.anno.gatk.list		#SNP annotation results of all filtered data		
├── *.raw.filter.snp.anno.vcf			#SNP annotation vcf results of all filtered data
├── *.raw.filter.snp.anno.gatk.vcf		#SNP annotation vcf results of all filtered data(only include annotated SNP mutation sites)
├── *.raw.filter.snp.vcf			#Filtered VCF file(all subsequent analysis is based on this file)
├── *.raw.filter.snp.DEG.snp			#Differential SNPs among samples
├── *.raw.filter.snp.snp			#SNP results
├── *.raw.filter.snp.snp.stat			#Statistics of SNP type of each sample
├── *.raw.filter.snp.anno.gatk.*.list		#Gene list of non-synonymous SNP mutations
├── SNP.anno.R*.pie.pdf				#Pie chart of SNP annotation statistics results(pdf format)
├── SNP.anno.R*.pie.png				#Pie chart of SNP annotation statistics results(png format)
├── SNP.anno_stat.xls				#Statistic table of SNP annotation reuslts
├── snp.genotype.info.cloud			#Statistics of SNP annotation and depth information of each sample
├── SNP.mutation.distribution.png		#Transition and transversion diagram of each sample
├── SNP.quality.distribution.png		#Distribution map of SNP quality(cumulative diagram of depth and distance)
├── SNP.venn.pdf				#Venn diagram of SNPs(pdf format), only if there are no more than five samples
└── SNP.venn.png				#Venn diagram of SNPs(png format), only if there are no more than five samples
INDEL/
├── all.sample.indel.length.distribution.pdf				#Distribution map of INDEL length(pdf format)
├── all.sample.indel.length.distribution.png				#Distribution map of INDEL length(png format)
├── all.sample.indel.length.stat					#INDEL length statistics
├── DEG_indel.stat							#Differential INDEL statistics
├── Indel.anno.R*.pie.pdf						#Pie chart of annotation statistics results(pdf format)
├── Indel.anno.R*.pie.png						#Pie chart of annotation statistics results(png format)
├── Indel.anno_stat.xls							#Statistic table of INDEL annotation results
├── indel.genotype.info.cloud						#Statistics of INDEL annotation and depth information of each sample
├── INDEL.venn.pdf							#Venn diagram of INDEL(pdf format), only if there are no more than five samples
├── INDEL.venn.png							#Venn diagram of INDEL(png format), only if there are no more than five samples
├── *.combine.Indel.dis.png						#Distribution map of INDEL length(total)
├── *.combine.Indel.stat						#INDEL distribution statistics in coding and non-coding areas
├── *.raw.filter.indel.anno.gatk.list					#INDEL annotation results of all sample
├── *.raw.filter.indel.anno.gatk.*.list					#Gene list of non-synonymous INDEL mutations
├── *.raw.filter.indel.anno.gatk.vcf					#INDEL annotation vcf results(only include annotated INDEL mutation sites)
├── *.raw.filter.indel.vcf						#Filtered VCF file(all subsequent analysis is based on this file)
└── *.raw.filter.indel.anno.vcf						#INDEL annotation vcf results of all samples
Circos/
├── R*.circos.png					#Distribution of each type of mutations in chromosome(png format)
└── R*.circos.svg					#Distribution of each type of mutations in chromosome(svg format)
Diff_analysis/
├── Allsample_Diffgene.list					#Statistics of mutated genes of all samples
├── Diff_gene.stat						#Statistics of mutant genes annotation of each type
└── R*
    ├── Cog_Anno						#Cog annotation results
    │   ├── *.Cog.classfy.png
    │   └── *.Cog.classfy.stat
    ├── go_enrichment						#GO annotation results
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
    ├── Graph							#topGO results
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
    ├── Integrated_Function.annotation.xls			#Gene annotation list of sample R*
    └── pathway							#KEGG results
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
            └── ko04141.png
Tips：
1. Files with size less than 50 Mb can be opened by text editors like Editplus/Ultraedit and also can be open with Microsoft Excel by change the name to "XXX.xls" format. Files with size larger than 50 Mb are recomonded to be handled with Pilotedit (http://www.pilotedit.com).
2. All files except figures(png/giff/jpeg/jpg...) can be opened with text editors.
Description of basic file：
********************************************************************************
1. AllSample_GC_Q.stat results description
SampleID				#Sample name
ReadSum					#Read number of samples
BaseSum					#Base number of samples
GC(%)					#GC content
N(%)					#N content
Q20(%)					#Q20 value
CycleQ20(%)				#CycleQ20 value
Q30(%)					#Q30 value

2. *.depth_stat.xls results description
#BMK ID					#BMK sample ID
Ave_depth				#Average coverage depth of samples
Cov_ratio_1X(%)				#Corresponding genomic coverage ratio at a depth of 1X
Cov_ratio_5X(%)				#Corresponding genomic coverage ratio at a depth of 5X
Cov_ratio_10X(%)			#Corresponding genomic coverage ratio at a depth of 10X

3. *.map_stat.xls results description
#BMK ID					#BMK sample ID
Total_reads				#Clean Reads number counted with two ends respectively, which read1 and read2 are counted as 2 reads
Mapped(%)				#Percentage of clean reads mapped to reference genome to total clean reads
Properly mapped(%)			#Both pair-end sequencing reads are mapped to reference genome and the distance matched the length of sequencing fragments

4. *.raw.filter.snp.anno.gatk.list results description
#CHROM					#Chromosome number
POS					#Position
REF					#Reference genome base type
ALT					#Mutant genes
QUAL					#Q value
R*_base					#All base types of sample R*
R*_alt_depth				#Sequencing depth of each base in sample R*
R*_reads_depth				#Reads depth in sample R* at this site
SNPEFF_EFFECT				#Structure annotation results of sample R*

5. *.raw.filter.snp.anno.vcf results description
#CHROM					#Chromosome number
POS					#Position
ID					#ID
REF					#Reference genome base type
ALT					#Mutant genes
QUAL					#Q value
FILTER					#Quality assessment result by GATK, if the result is reliable, the message is "PASS", the others will be ".".
INFO					#The variation information; AC: Allele count in genotypes, for each ALT allele, in the same order as listed; AF: Allele Frequency, for each ALT allele, in the same order as listed; AN: Total number of alleles in called genotypes; For diploid species, GT="0/1" means the allele is heterozygous, only one allelic loci mutated. Its Allele Frequency(AF) is 0.5and AN is 2; If the genotype is "1/1", means the allele in this sample is homozygous, its AC is 2, AF is 1 and AN is 2; DP: Approximate read depth; some reads may have been filtered; FS: Phred-scaled p-value using Fisher's exact test to detect strand bias, the smaller, the better; 
FORMAT					#GT: The sample genotype and numbers are separated by "/"; AD: Allelic depths for the ref and alt alleles in the order listed, the first number is the depth of the reference allele, the second number shows the depth for the 1st variant; DP: The read depth in this allele; GQ: Genotype Quality, phred format, indicating the possibility of the existence of the genotype at the site, the higher, the more possibility of this genotype. Phred value = -10 * log (1-p) ，p means the probability of this allele. PL="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification", the three specified genotypes are (0/0,0/1,1/1), and the total probability of the three genotypes is 1. Unlike before, the higher the value, the less likely it is to be the genotype. Phred = -10 * log (p), where p is the probability of genotype existence.
R*					#The genotype of the sample of R*, "0/0" means no mutation; "0/1" means heterozygous mutation；"1/1" means  homozygous mutation.

-------------------------------------------------------------------------------
copyright (c) BMK 2020
last update 2019.12.30
