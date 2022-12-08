# Resequencing
## GATK4教程
https://hpc.nih.gov/training/gatk_tutorial/preproc.html#preproc-tools<br>

## variant calling pipeline
VarDict: https://github.com/AstraZeneca-NGS/VarDict<br>

## 人类基因组结构变异数据库（GATK4质控）
### 1000 Genomes：https://www.internationalgenome.org/home<br>
The 1000 Genomes Project created a catalogue of common human genetic variation, using openly consented samples from people who declared themselves to be healthy. The reference data resources generated by the project remain heavily used by the biomedical science community.
### dbSNP：https://www.ncbi.nlm.nih.gov/snp/<br>
dbSNP contains human single nucleotide variations, microsatellites, and small-scale insertions and deletions along with publication, population frequency, molecular consequence, and genomic and RefSeq mapping information for both common variations and clinical mutations.
### hapmap：https://www.genome.gov/10001688/international-hapmap-project
HapMap (short for "haplotype map") is the nickname of the International HapMap Project, an international project that seeks to relate variations in human DNA sequences with genes associated with health. A haplotype is a set of DNA variations, or polymorphisms, that tend to be inherited together. A haplotype can refer to a combination of alleles or to a set of single nucleotide polymorphisms (SNPs) found on the same chromosome. The HapMap describes common patterns of genetic variation among people. 
### ClinVar: https://www.ncbi.nlm.nih.gov/clinvar<br>
ClinVar 数据库整合了基因组变异及其与人类健康关系的信息
### dbVar：https://www.ncbi.nlm.nih.gov/dbvar/
dbVar is NCBI's database of human genomic Structural Variation — large variants >50 bp including insertions, deletions, duplications, inversions, mobile elements, translocations, and complex variants.
### ClinGen：https://clinicalgenome.org/
ClinGen is a National Institutes of Health (NIH)-funded resource dedicated to building a central resource that defines the clinical relevance of genes and variants for use in precision medicine and research.
### GnomAD: http://www.gnomad-sg.org/
The Genome Aggregation Database (gnomAD) is a resource developed by an international coalition of investigators, with the goal of aggregating and harmonizing both exome and genome sequencing data from a wide variety of large-scale sequencing projects, and making summary data available for the wider scientific community.


## bam文件检查
VerifyBamID2：https://github.com/Griffan/VerifyBamID 检查BAM文件中的读数是否与特定样品的先前基因型匹配

## 结构变异
LUMPY：https://github.com/arq5x/lumpy-sv<br>
delly：https://github.com/dellytools/delly<br>
### nanopore SV
### PacBio SV
### both
## 单倍型推断
Shapeit2: https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html<br>
Shapeit4: https://github.com/odelaneau/shapeit4<br>
hapcut2: https://github.com/vibansal/hapcut2<br>
nPhase: https://github.com/OmarOakheart/nPhase<br>

## 基因型填充
Impute2: https://mathgen.stats.ox.ac.uk/impute/impute_v2.html<br>
Beagle: https://faculty.washington.edu/browning/beagle/beagle.html<br>
eagleimp: https://github.com/ikmb/eagleimp<br>
## CNV分析
### 靶向测序（PCR或者芯片捕获）
### 全外显子
### 全基因组

### 分析软件
#### FREEC：可用于全基因组和全外显子，make不好安装
https://github.com/BoevaLab/FREEC

## snp注释
SnpEff: https://pcingola.github.io/SnpEff<br>
VEP: https://github.com/Ensembl/ensembl-vep<br>
### ANNOVAR
人、小鼠、蠕虫、果蝇、酵母等基因组结构变异功能注释，主要可以做如下3类注释
1 基于基因的注释<br>
确认SNPs和CNVs造成的编码蛋白氨基酸的变化和影响<br>
2 基于区域的注释<br>
确定基因组特定区域的变异，如转录因子结合位点、片段重复区域、GWAS hits区域、组蛋白修饰区域等<br>
3 Filter-based annotation，SNPs和indel<br>
确认变异是否在特定的数据库中被记录描述，如dbDSNP，1000 Genome Project, gnomAD等，计算有害突变等得分<br>
#### 查看可用数据库
```
#https://www.biostars.org/p/196985/
perl annotate_variation.pl -downdb -webfrom annovar avdblist humandb/ -buildver hg38
```
#### 下载数据
```
perl annotate_variation.pl -downdb clinvar_20220320 -webfrom annovar humandb/ -buildver hg19
```

#### VCF格式转换为annovar输入格式
```
perl convert2annovar.pl -format vcf4 case_22BY12800_filter.vcf > case_22BY12800_filter.vcf_variant.avinput
```

## 问题
### 修改染色体名字，关于为什么只需要修改BAM文件的header？
https://www.biostars.org/p/386231/
