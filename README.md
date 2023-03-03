# Resequencing
## GATK4教程
https://hpc.nih.gov/training/gatk_tutorial/preproc.html#preproc-tools<br>
## 群体遗传学
https://speciationgenomics.github.io/
## 统计遗传学
https://genome.sph.umich.edu/wiki/Abecasis_Lab<br>
An Introduction to Statistical Genetic Data Analysis<br>
Handbook of statistical genomics<br>
http://faculty.washington.edu/tathornt/SISG2019.html<br>

### lab
https://yanglab.westlake.edu.cn/ <br>

## variant calling pipeline
### VarDict:
https://github.com/AstraZeneca-NGS/VarDict<br>
### bcftools:
http://samtools.github.io/bcftools/howtos/cnv-calling.html

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

truvari: https://github.com/acenglish/truvari
survivor: https://github.com/fritzsedlazeck/SURVIVOR

## 基因组重排
gridss：https://github.com/PapenfussLab/gridss

## 单倍型推断
Shapeit2: https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html<br>
Shapeit4: https://github.com/odelaneau/shapeit4<br>
hapcut2: https://github.com/vibansal/hapcut2<br>
nPhase: https://github.com/OmarOakheart/nPhase<br>

## 基因型填充
Impute2: https://mathgen.stats.ox.ac.uk/impute/impute_v2.html<br>
Beagle: https://faculty.washington.edu/browning/beagle/beagle.html<br>
Minimac4: https://github.com/statgen/Minimac4<br>
eagleimp: https://github.com/ikmb/eagleimp<br>
eagle2: https://alkesgroup.broadinstitute.org/Eagle/<br>

## 连锁不平衡（Linkage disequilibrium）
### LDSC
https://github.com/bulik/ldsc<br>



## CNV分析
### 靶向测序（PCR或者芯片捕获）
### 全外显子
### 全基因组

### 分析软件
#### FREEC：适用全基因组和全外显子，make不好安装
https://github.com/BoevaLab/FREEC
下载release版本，解压安装即可
#### CNVpytor: 适用全基因组
https://github.com/abyzovlab/CNVpytor
#### cn.mops: 适合全外显子和panel
https://bioconductor.org/packages/release/bioc/vignettes/cn.mops/inst/doc/cn.mops.pdf
#### xhmm: 适合全外显子
https://zzz.bwh.harvard.edu/xhmm/index.shtml
#### ONCOCNV
https://github.com/BoevaLab/ONCOCNV
#### CODEX2
https://github.com/yuchaojiang/CODEX2
#### pureCN 需要100X数据量
https://github.com/lima1/PureCN
#### scarHRD
https://github.com/sztup/scarHRD
#### facet
https://github.com/mskcc/facets

## snp注释
SnpEff: https://pcingola.github.io/SnpEff<br>
VEP: https://github.com/Ensembl/ensembl-vep<br>
### vcfanno
https://github.com/brentp/vcfanno
### ANNOVAR
https://annovar.openbioinformatics.org/en/latest/<br>
人、小鼠、蠕虫、果蝇、酵母等基因组结构变异功能注释，主要可以做如下3类注释<br>
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

#### 单个数据库注释
```
perl annotate_variation.pl -filter -dbtype clinvar_20220320 --thread 30 -buildver hg19 case_22BY12800_filter.vcf_variant.avinput humandb/
```

#### 联合注释
```
# annotation
perl table_annovar.pl example/ex1.avinput humandb/ \
-buildver hg19 \
-out myanno \
-remove \
-protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a \
-operation g,r,f,f,f \
-nastring . \
-csvout

# -buildver hg19 表示使用的参考基因组版本为hg19
# -out myanno 指定输出文件前缀为myanno
# -remove 表示删除中间文件
# -protocol 后跟注释来源数据库名称，每个protocal名称或注释类型之间只有一个逗号，并且没有空白
# -operation 后跟指定的注释类型，和protocol指定的数据库顺序是一致的，g代表gene-based、r代表region-based、f代表filter-based
# -nastring . 表示用.替代缺省值
# -csvout 表示最后输出.csv文件

# 可用-vcfinput选项直接对vcf文件进行注释
perl table_annovar.pl example/ex2.vcf humandb/ -buildver hg19 -out myanno -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp138,ljb23_all -operation g,r,r,f,f,f,f -nastring . -vcfinput

```
## runs of homozygosity 和杂合性缺失
https://cran.r-project.org/web/packages/detectRUNS/vignettes/detectRUNS.vignette.html<br>



## 亲缘关系分析
akt：https://github.com/Illumina/akt<br>
king：https://www.kingrelatedness.com<br>

## PRS 多基因风险分数
https://github.com/getian107/PRScs<br>

## GWAS
https://gwaslab.com/<br>
https://pbreheny.github.io/adv-gwas-tutorial/index.html<br>
### EMMAX
https://genome.sph.umich.edu/wiki/EMMAX

## 合并vcf
### 相当于cat，例如把分染色体call的vcf文件合并到一块
```
java -jar ~/software/snpEff/SnpSift.jar split -j  *bed.vcf > merged.vcf
java -jar ~/software/snpEff/SnpSift.jar sort *bed.vcf > merged.vcf
#不同在于注释行的差异，其他均一致
```
https://github.com/reneshbedre/bioinfokit<br>

##  loss of heterozygosity(杂合性缺失)
https://www.bioconductor.org/packages/release/bioc/vignettes/TitanCNA/inst/doc/TitanCNA.pdf

## 问题
### 修改染色体名字，关于为什么只需要修改BAM文件的header？
https://www.biostars.org/p/386231/
### deepvariant variant call 结果中PASS和RefCall的解释？
https://github.com/google/deepvariant/issues/278
### bcftools 给VCF文件增加Tags
https://github.com/samtools/bcftools/issues/1731
### b allele frequency
www.biostars.org/p/254848/
### VCF文件中的星号（*）/<*>的意思和区别
https://gatk.broadinstitute.org/hc/en-us/articles/360035531912?id=11029<br>
http://samtools.github.io/hts-specs/VCFv4.3.pdf 页数8<br>
### SNP的LD剪枝与聚集 LD pruning & clumping
https://www.biostars.org/p/343818/
https://zhuanlan.zhihu.com/p/373217037
