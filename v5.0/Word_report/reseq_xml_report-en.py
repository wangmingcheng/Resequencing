# -*- coding: utf-8 -*- 
import sys, os, argparse, glob, os.path,time,glob
reload(sys)
sys.setdefaultencoding('utf8')
import numpy as np
import math
import re
import xml.dom.minidom 
from xml.dom.minidom import  parseString, getDOMImplementation
from collections import OrderedDict
import urllib
import shutil
import codecs
import subprocess
Bin=os.path.split(os.path.realpath(__file__))[0]
t1=0
t2=0
t3=0    
class parseXML(object): 
    def __init__(self, file): 
        self.f = file
    def __iter__(self): 
        return self 
    def parse(self):
        file=open(self.f,'r')
        
    def next(self): 
        r = [self.s,self.e] 
        self.s, self.e = self.bin+self.s, self.e + self.bin 
        return r 
def cp_template(o,s=None):
    if os.path.exists(o+'/src'): shutil.rmtree(o+'/src')
    if s:
        shutil.copytree(s,o+'/src')
    else:
        shutil.copytree(Bin+'/src',o+'/src')
    
def convert_pic():
    shutil.copyfile("oldfile","newfile")
def addTable(tmp,isfirst):
    result=""
    if  isfirst:
        for i in tmp:
            result+="<th>%s</th>"%i
    else:
        for i in tmp:
            result+="<td>%s</td>"%i
    return result
    
def timenow():
    """return current time as a string
    """
    return time.strftime('[%d/%m/%Y %H:%M:%S]', time.localtime(time.time()))
def Indent(dom, node, indent = 0):
    children = node.childNodes[:]
    if indent:
        text = dom.createTextNode('\n' + '\t' * indent)
        node.parentNode.insertBefore(text, node)
    if children:
        if children[-1].nodeType == node.ELEMENT_NODE:
            text = dom.createTextNode('\n' + '\t' * indent)
            node.appendChild(text)
        for n in children:
            if n.nodeType == node.ELEMENT_NODE:
                Indent(dom, n, indent + 1)
def addElement(d,r,l,**args):
    global t1
    global t2
    global t3
    if(re.match("^h1$", l)):
        t1+=1
        t2=0
        t3=0
    if(re.match("^h2$", l)):
        t2+=1
        t3=0
    if(re.match("^h3$", l)):
        t3+=1
    if t2==0:
        myleble="%s"%(t1)
    elif t3==0:
        myleble="%s.%s"%(t1,t2)
    else:
        myleble="%s.%s.%s"%(t1,t2,t3)
        
    e=d.createElement(l)
    for k,v in args.iteritems():
        if(re.match("^h[1-3]$", l) and re.match("^name$", k)) :
            e.setAttribute(k,"%s %s"%(myleble,v))
        else:
            e.setAttribute(k,v)
    r.appendChild(e)
    return e
def addElementListImage(d,r,l,p,analysis_path,desc=None):
    if desc==None:
        desc=r.getAttribute('desc')
    for i in sorted(p):
        dir,suffix=os.path.splitext(i)
        name=os.path.basename(dir)
        addElement(d, r,l ,name=name ,type="type1", path=i.replace(analysis_path+'/',""),desc=desc)

def get_reflen(refinfo):
    print refinfo
    r=re.split('\s+', refinfo[-1])[1]
    print r
    m=re.match('(\d*\.?\d*?)([kKmMGg]?)', r)
    if m:
        reflen=m.group(1)
        k=1
        if not m.group(2)=='':
            
            if re.match('[kK]',m.group(2)):k=1000
            if re.match('[mM]',m.group(2)):k=1000000
            if re.match('[gG]',m.group(2)):k=1000000000
        return float(reflen)*k
    else:
        sys.stderr.write("Can't find reference length please check your config file\n")
        exit(1)

parser = argparse.ArgumentParser(description='This script was used to generate a xml report from reseq analysis.')
parser.add_argument('-detail', '--detail', dest='detail', required=True, help='input detail config  file')
parser.add_argument('-data', '--data', dest='data', required=True, help='input data config  file')
parser.add_argument('-l', '--tableLineNum', dest='tableLineNum', required=False,type=int,default=100, help='Input max table line numbers to show, default:100')
parser.add_argument('-n','--name',dest='name',required=False,default='configtest',help='specify the output file prefix,default is "configtest"')
args = parser.parse_args()
###################################################################
args.detail=os.path.abspath(args.detail)
args.data = os.path.abspath(args.data)
lineNum=args.tableLineNum
sys.stdout.write("\n%s:Start to generate xml file...\n\n"%timenow())
###########################################################################
f=codecs.open(args.detail,'rU','utf-8')
ref_info=[]
for i in f:
    i=i.strip()
    
    if i.startswith('analysis_dir'):analysis_dir=os.path.abspath(re.split(u'\s+', i)[1])
    if i.startswith('Sample_Project'):Sample_Project=re.split(u'\s+', i)[1]
    if i.startswith('ref_url'):ref_url=re.split(u'\s+', i)[1]
    if i.startswith('Project'):Project=re.split(u'\s+', i)[1]
    if i.startswith('title_Project'):title_Project=re.split(u'\s+', i)[1]
    if i.startswith('%'):ref_info.append(re.split(u'\s+', i,1)[1])
    if i.startswith('platform'):Sequencing_platform=re.split(u'\t', i)[1]
    else:
        Sequencing_platform='Illumina'
    if i.startswith('read_length'):Read_length=re.split(u'\s+', i)[1]
    if i.startswith('delete'):delete=re.split(u'\s+',i)[1]
f.close()

fdata = codecs.open(args.data,"rU","utf-8")
sampleInfo=[]
for i in fdata:
	i=i.strip()
	if i.startswith("Sample"):sampleInfo.append(re.split(u'\s+',i,1)[1])
fdata.close()

web_dir = analysis_dir;
if not os.path.exists(web_dir):
    os.mkdir(web_dir)
f=codecs.open(web_dir+'/sampleinfo.txt','w','utf-8')
f.write("BMK ID\t"+"Client ID\n")
for i in sampleInfo:
    f.write(i.encode('utf-8')+'\n')
f.close()
f=codecs.open(web_dir+'/refinfo.txt','w','utf-8',)
for i in ref_info:
    f.write(i.encode('utf-8')+'\n')
reflen=get_reflen(ref_info)
f.close()
sampleNum=len(sampleInfo)
print sampleNum
#####################################get project guideLine#############################################
cleanData=0.
averageQ30=0.
averagemapRate=0.
averageDepth=0.
coverageRate=0.
averageSNP=0
averageINDEL=0
averageSV=0

f=open(analysis_dir+'/Up_load/Dataassess/sample_data_assess.list','r')
Clean_Base=0.
Q30=0.
for i in f:
    i=i.strip()
    if "BMK" in i:continue
    tmp=re.split('\t', i)
    Clean_Base+=float(tmp[2])
    Q30+=float(tmp[4])
    
averageQ30=Q30/sampleNum
cleanData=Clean_Base/1000000000
f.close()

f=open(analysis_dir+"/Up_load/Mapping/%s.map_stat.xls"%Project,'r')
Mapped=0.
for i in f:
    i=i.strip()
    if "BMK" in i:continue
    tmp=re.split('\t', i)
    Mapped+=float(tmp[2])
averagemapRate=Mapped/sampleNum
f.close()

f=open(analysis_dir+"/Up_load/Mapping/%s.depth_stat.xls"%Project,'r')
Ave_depth=0.
Cov_ratio_1X=0.
for i in f:
    i=i.strip()
    if "BMK" in i:continue
    tmp=re.split('\t', i)
    Ave_depth+=float(tmp[1])
    Cov_ratio_1X+=float(tmp[2])
    
averageDepth=Ave_depth/sampleNum
coverageRate=Cov_ratio_1X/sampleNum
f.close()

#####################################################################################
impl = xml.dom.minidom.getDOMImplementation()
dom = impl.createDocument(None, 'report', None)
root = dom.documentElement
addElement(dom,root,'report_version',value="v1.3")
addElement(dom,root,'report_name',value=title_Project)
addElement(dom,root,'report_code',value=Project)
addElement(dom,root,'report_user',value="XXXX")
addElement(dom,root,'report_user_addr',value="XXXX")
addElement(dom,root,'report_time',value="XXXX/XX/XX;XXXX/XX/XX;XXXX/XX/XX;%s" % time.strftime('%Y/%m/%d', time.localtime(time.time())))

if os.path.exists(analysis_dir+'/Up_load/SV') and os.path.exists(analysis_dir + '/Up_load/CNV') and os.path.exists(analysis_dir + '/Up_load/Diff_analysis'):
    addElement(dom,root,'report_abstract',value="<p class= \"p-abstract\" >Analysis content:</p><p class=\" p-abstract\" >Completed sequencing of %s samples, the specific analysis contents are as follows, data evaluation: statistics of sequencing data amount, sequencing data quality and GC content. Alignment with genome: statistics of alignment rate, genome coverage and the depth of genome coverage. Variation detection and annotation: detection and annotation of SNP, InDel, SV and CNV. DNA level mutation gene analysis: detection of genes with SNP non-synonymous mutation, InDel mutation, SV mutation and CNV mutation in the coding region. Gene annotation: DNA level mutation genes were annotated by KEGG, GO, COG, NR, SwissProt database. </p><p class=\" p-abstract\" > Summary of results：</p><p class=\" p-abstract\" > Amount of data used in this analysis is %.2fGbp Clean Data, it's Q30 reaches %.2f%%. The average alignment rate between sample and reference genome is %.2f%%, the mean coverage depth is %dX, the genome coverage is %.2f%% (at least covered by one base). The detailed mutation detection results (SNP, Indel) are shown in the report below. SNP non-synonymous mutations, InDel mutations, SV and CNV genes were detected between the sample and the reference genome, and the genes at DNA level variations were annotated by databases such as KEGG, GO, COG, NR, SwissProt, etc. (for details, the contract shall prevail). </p>"%(sampleNum,cleanData,averageQ30,averagemapRate,averageDepth,coverageRate) )
else:
    addElement(dom,root,'report_abstract',value="<p class= \"p-abstract\" > Analysis content:</p><p class=\" p-abstract\" > Complete resequencing for %s samples, the detailed analysis content is as follows, data evaluation: statistics of sequencing data amount, sequencing data quality and GC content. Alignment with genome: statistics of alignment rate, genome coverage and the depth of genome coverage. Variation detection and annotation: detection and annotation of SNP, InDel. </p><p class=\" p-abstract\" > Summary of results：</p><p class=\" p-abstract\" > Amount of data used in this analysis is %.2fGbp Clean Data,  it's Q30 reaches %.2f%%. The average alignment rate between sample and reference genome is %.2f%%，the mean coverage depth is %dX，the genome coverage is %.2f%% (at least covered by one base). The detailed mutation detection results (SNP, Indel) are shown in the report below. </p>"%(sampleNum,cleanData,averageQ30,averagemapRate,averageDepth,coverageRate) )
    
addElement(dom,root,'h1',name="Project basic information" ,type="type1" ,desc="Project basic information")
addElement(dom,root,'p',type="type1" ,desc="(1)the corresponding ID of sample information:")
addElement(dom,root,'table',name="sample information" ,type="full" ,desc="Note：BMK ID：BMK will set a unique ID for a sample, both the library construction and subsequent bioinformatics will use this ID." ,path="sampleinfo.txt")
addElement(dom,root,'p',type="type1", desc="(2) Information of reference genome：")
addElement(dom,root,'p',type="type1" ,desc="Sequenced species：%s；Download URL：<a href=\"%s\" title=\"click\" class=\"mylink\" target=\"_blank\" >%s</a>"%(Project,ref_url,ref_url) )
addElement(dom, root,'table',name="Genome information" ,type="full", desc="", path="refinfo.txt")
addElement(dom, root,'h1',  name="Project process", type="type1", desc="Project process")
addElement(dom, root,'h2' ,name="Experiment workflow" ,type="type1" ,desc="Experiment workflow")
if Sequencing_platform == 'Illumina':
	addElement(dom, root,"p" ,type="type1" ,desc="Experiment workflow follows standard protocol provided by Illumina Corporation, including sample quality control, library construction, library quality control and library sequencing, etc. The specific flow chart is as follows." )
	addElement(dom, root,"pic" ,name="Experiment workflow" ,type="type1", desc=" ", path="src/images/experiment_process_en.png" )
	addElement(dom, root,"p" ,type="type1" ,desc="Library construction after DNA sample qualification includes:DNA fragmentation, End repair, Add “A” bases to the end, Ligate adapter, Production size selection, PCR reaction, PCR production fragment size selection and Library quality control. qualified libraries were sequenced by using %s."%Sequencing_platform)
else:
	addElement(dom, root,"p" ,type="type1" ,desc="Experiment workflow follows standard protocol provided by BGI Corporation, including sample quality control, library construction, library quality control and library sequencing, etc. The specific flow chart is as follows." )
	addElement(dom, root,"pic" ,name="Experiment workflow" ,type="type1", desc=" ", path="src/images/experiment_process.BGI_en.png" )
	addElement(dom, root,"p" ,type="type1" ,desc="Library construction after DNA sample qualification includes:DNA fragmentation, End repair, Add “A” bases to the end, Ligate adapter , PCR reaction, the DNA fragments amplified by PCR single chain cycling to form circular DNA, DNA nanoball (DNB) was prepared by using circular DNA as a template through a unique linear amplification mode (roll ring amplification). The prepared DNB library was combined with the array sites on the sequencing chip for library sequencing." )
addElement(dom, root,"h2" ,name="Bioinformatics workflow" ,type="type1", desc="Bioinformatics workflow" )
addElement(dom, root,"p", type="type1", desc="Make quality evaluation and filtering for the raw reads (Paired ends) obtained from sequencing to get Clean Reads, which can be used for subsequent bioinformatics. Align Clean Reads with reference genome, perform variation detection and annotation for SNP, InDel, etc. of aligned result to perform DNA-level DEGs mining and DEGs function annotation." )
addElement(dom, root,"p", type="type1", desc="Bioinformatics workflow for resequencing is shown as below：" )
if os.path.exists(analysis_dir+'/Up_load/SV') and os.path.exists(analysis_dir + '/Up_load/CNV') and os.path.exists(analysis_dir + '/Up_load/Diff_analysis'):
    addElement(dom, root,"pic", name="Resequencing bioinformatics flow chart", type="type1", desc=" ", path="src/images/Analysis_process_en.png" )
else:
    addElement(dom, root,"pic", name="Resequencing bioinformatics flow chart", type="img-width-normal", desc=" ", path="src/images/Analysis_process2_en.png" )
addElement(dom, root,"h1" ,name="Bioinformatics analysis method and result" ,type="type1", desc="Bioinformatics analysis method and result" )
addElement(dom, root,"h2", name="Sequencing data quality control" ,type="type1" ,desc="Sequencing data quality control" )
addElement(dom, root,"h3", name="Sequencing data introduction",type="type1", desc="Sequencing data introduction" )
addElement(dom, root,"p" ,type="type1" ,desc="Raw image data files obtained from high-throughput sequencing （%s sequencing platform）were converted into raw sequenced reads through base calling, which was called as Raw Data or Raw Reads and stored as FASTQ（fq）file, including sequence information and the corresponding sequence quality information of Reads. The randomly selected real data of sequenced sample are as follows：" %Sequencing_platform)
if Sequencing_platform == 'Illumina':
	addElement(dom, root,"pic" ,name="Figure Fastq format", type="type1" ,desc="Note：Each Read in the FASTQ format file is described by four lines, the first of which begins with “@” followed by Illumina Sequence Identifiers and description (optional part); The second line is the sequence of bases; The third line begins with “+”, followed by the Illumina sequencing identifier (optional part); The fourth line shows the sequencing quality of the corresponding sequence." ,path="src/images/fastq_format.png" ) 
	addElement(dom, root,"p" ,type="type1" ,desc="Details of the Illumina Sequence Identifiers are as follows：" )
	addElement(dom, root,"table", name="Detailed information table of Illumina Sequence Identifiers", type="full", desc="", path="src/images/table3-en" ) 
	addElement(dom, root,"p", type="type1" ,desc="By using the ASCII value corresponding to each character in the fourth line, the sequencing quality value corresponding to bases in the second line is obtained. If the sequencing error rate is represented by e and the base mass quality value of %s is represented by Qphred, then the following relationship exists："%Sequencing_platform)
	addElement(dom, root,"pic", name="Figure Formula 1", type="img-width-normal" ,desc=" " ,path="src/images/formula1.png" )
else:
	addElement(dom, root,"pic" ,name="Figure Fastq format", type="type1" ,desc="Note：Each Read in the FASTQ format file is described by four lines, the first of which begins with “@” followed by BGI Sequence Identifiers and description (optional part); The second line is the sequence of bases; The third line begins with “+”, followed by the BGI sequencing identifier (optional part); The fourth line shows the sequencing quality of the corresponding sequence." ,path="src/images/fastq_format_BGI.png" ) 
	addElement(dom, root,"p" ,type="type1" ,desc="Details of the BGI Sequence Identifiers are as follows：" )
	addElement(dom, root,"table", name="Detailed information table of BGI Sequence Identifiers", type="full", desc="", path="src/images/table3_BGI-en" )
	addElement(dom, root,"p", type="type1" ,desc="By using the ASCII value corresponding to each character in the fourth line, the sequencing quality value corresponding to bases in the second line is obtained. If the sequencing error rate is represented by e and the base mass quality value of %s is represented by Qphred, then the following relationship exists："%Sequencing_platform)
	addElement(dom, root,"pic", name="Figure Formula 1", type="img-width-normal" ,desc=" " ,path="src/images/formula1.png" )
addElement(dom, root,"p" ,type="type1" ,desc="The concise correlation between the sequencing error rate from Illumina Casava version 1.8 and the sequencing quality value is shown in the following table：" )
addElement(dom, root,"table" ,name="The correlation between sequencing error rate and sequencing quality value", type="full" ,desc="Note：Analysis software for Base Calling：Illumina Casava 1.8 version；Sequencing parameters：Paired end；Sequencing reads length：%sbp（or in cycles)"%Read_length, path="src/images/table4-en" ) 
addElement(dom, root,"h3", name="data quality statistics" ,type="type1", desc="data quality statistics" )
addElement(dom, root,"p", type="type1" ,desc="Raw Sequenced Reads or Raw Reads were obtained from sequencing, which also have low quality Reads with adaptors. To ensure bioinformatics quality, Raw Reads were filtered to get Clean Reads for subsequent bioinformatics. The main steps of data filtering are as follows :(1) remove the reads with adapter; (2) filter reads with N content over 10%; (3) remove reads whose base value (with quality value less than 10) is more than 50%." )
addElement(dom, root,"p", type="type1", desc="The evaluation results of sequencing output data of each sample are shown in the following table：" )
addElement(dom, root,"table", name="Statistics of evaluation of sample sequencing data" ,type="full" ,path="Up_load/Dataassess/sample_data_assess.list",desc="Note：BMK_ID：The uniform sample ID set by BMK；Clean_Reads：Number of filtered reads； Clean_Base：number of filtered bases, which equals to the number of Clean_Reads times the sequence length；Q20(%)：The percentage of bases with a mass value greater than or equal to 20 in the total number of bases；Q30(%)：The percentage of bases with a mass value greater than or equal to 30 in the total number of bases；GC(%)：Sample GC content, the percentage of G and C type bases in the total bases."  ) 
addElement(dom, root,"h3", name="Base sequencing quality distribution", type="type1", desc="Base sequencing quality distribution" )
addElement(dom, root,"p" ,type="type1", desc="The sequencing error rate of each base was obtained by conversion of sequencing Phred score (Qphred) by formula 1, while the Phred value was calculated by a probabilistic model for predicting the occurrence of error in base discrimination during Base Calling. The correlation is shown in the following table：" )
addElement(dom, root,"table", name="Prediction of error probability for base discrimination", type="full", desc="", path="src/images/table31-en" ) 
addElement(dom, root,"p" ,type="type1", desc="When %s sequencing system was used in sequencing, firstly chip preparation was performed for library, aiming at fixing library DNA template on chip. During fixation of DNA template, each DNA molecule will form a cluster, one cluster is one sequencing site, very few clusters will have physical overlap during fixation, during sequencing, software will use the first four bases to analyse and distinguish these overlapping sites, then they can be separated to ensure that one DNA molecule will be sequenced in one site, thus the error rate of first few bases at 5' end of sequences is relatively higher. In addition, the sequencing error rate will increase as the length of sequenced reads increases, which is caused by the consumption of chemical reagents in the sequencing process. Therefore, during analysis of quality distribution of base sequencing, the quality value of the sample's base quality distribution in the first four bases and the last a dozen or so bases would be lower than that of the intermediate sequencing bases, but the quality value would be higher than Q30%%. According to the relationship between the quality value and the error rate, we converted the quality value into the error rate and plotted the error rate distribution diagram as follows："%(Sequencing_platform))
pic_list=addElement(dom, root,"pic_list" ,name="Sample base error rate distribution" ,type="type1", desc="Note：X-coordinate shows base position of reads, Y-coordinate shows the error rate of single base. The first %sbp shows the error rate distribution of the first end of paired-ends sequences, the last %sbp shows the error rate distribution of another end."%(Read_length,Read_length))
qualityImages=glob.glob(analysis_dir+'/Up_load/Dataassess/*.quality.png')
addElementListImage(dom,pic_list,'pic',qualityImages,web_dir)
addElement(dom, root,"h3", name="Base type distribution" ,type="type1", desc="Base type distribution" )
addElement(dom, root,"p", type="type1", desc="Base type distribution check is used to detect whether there is AT or GC separation, which may be caused by sequencing or library construction, and may affect subsequent analysis. Sequences used in high-throughput sequencing are DNA fragments after random breaking of genome. Due to that the distribution of sites on the genome is approximately uniform, and the contents of G/C and A/T are also approximately uniform, so according to the law of large numbers, in each sequencing cycle, the GC and AT contents should be equal to each other and also equal to the GC and AT contents of the genome respectively. Similarly, due to the relationship of overlapping clusters, the first several bases of samples have higher fluctuant AT and GC, which is higher than other sequencing regions, while the contents of GC and AT in other regions are equal, and the distribution is uniform without separation, as shown in the figure below：" )
pic_list=addElement(dom, root,"pic_list" ,name="Proportion distribution of each base of the sample" ,type="type1", desc="Note：X-coordinate shows base location of reads, Y-coordinate shows the proportion of bases; different colors represent different base types, green for base G, blue for base C, red for base A, purple for base T, and gray for base N that is not recognized by sequencing. The first %sbp shows the base distribution of the first end of paired-ends sequences, the last %sbp shows the base distribution of another end. Each cycle represents each base sequenced. For example, the first cycle represents the distribution of A, T, G, C and N of all sequencing reads at the first base. This figure shows that the AT and CG bases are basically not separated, the curve is relatively flat, so the sequencing result is normal."%(Read_length,Read_length))
acgtnImages=glob.glob(analysis_dir+'/Up_load/Dataassess/*.acgtn.png')
addElementListImage(dom,pic_list,'pic',acgtnImages,web_dir)

addElement(dom, root,"h2", name="Statistics of alignment to reference genome", type="type1", desc="Statistics of alignment to reference genome" )
addElement(dom, root,"p", type="type1", desc="Sequencing reads obtained by resequencing need to be relocated to the reference genome for subsequent variation analysis. bwa<a href=\"#ref1\">[1]</a> software is mainly used to align the short sequences obtained by second-generation high-throughput sequencing （Sequencing platforms such as %s） to reference genome. Clean Reads were located on the reference genome by alignment, then the sequencing depth, genome coverage and other information of each sample were collected, and the variation was detected."%(Sequencing_platform))

addElement(dom, root,"h3",name="Alignment result statistics",type="type1",desc="Alignment result statistics" )
addElement(dom, root,"p", type="type1", desc="Reads were aligned to reference genome, then information such as samples alignment rate was collected. Alignment rate: the proportion of Clean Reads located on the reference genome to the total Clean Reads. If the appropriate reference genome is selected and there is no contamination in experiment, the alignment rate of sequencing reads will be higher than 70%. In addition, the alignment efficiency is related to the relationship between sequenced species and the reference genome, the assembly quality of reference genome and the reads sequencing quality, the closer the relationship is, the more complete the assembly of the reference genome is, the higher quality the reads sequencing is, the more reads that can be located to the reference genome, the higher the alignment rate is." )
addElement(dom, root,"p", type="type1", desc="The statistical method of alignment is as follows：" )
addElement(dom, root,"p",type="type1", desc="Clean Reads：Collect data file of Clean Reads, every four lines are as one unit, paired-ends data are collected separately, i.e. reads1 and reads2 will be counted as two reads.")
addElement(dom, root,"p",type="type1", desc="Mapped(%)：Percentage of number of Clean Reads located to the reference genome to the number of all Clean Reads, it was realized by samtools flagstat command.")
addElement(dom, root,"p",type="type1", desc="Properly mapped(%)：length distribution of paired-ends sequences that are located to the reference genome and their distance matches sequencing fragments, it was realized by samtools flagstat command.")

addElement(dom, root,"p", type="type1", desc="Sample alignment result is shown in the following table:" )
addElement(dom, root,"table" ,name="Alignment result statistics", type="full",  path="Up_load/Mapping/%s.map_stat.xls"%Project,desc="Note：BMK ID：The uniform sample ID set by BMK；Total_reads：number of Clean Reads；Mapped(%)：Percentage of number of Clean Reads located to the reference genome to the number of all Clean Reads；Properly_mapped(%)：length distribution of paired-ends sequences that are located to the reference genome and their distance matches sequencing fragments." ) 
addElement(dom, root,"h3", name="Insersion distribution statistics", type="type1", desc="Insersion distribution statistics" )
addElement(dom, root,"p", type="type1", desc="By detecting the start and end positions of the paired-ends sequences on the reference genome, the actual size of sequenced fragments from sample DNA can be obtained, i.e. the insert size, which is an important parameter in bioinformatics. The distribution of insertion size generally conforms to the normal distribution, it has only one single peak, the Insertion Size distribution map shows the length distribution of insertions of each sample. The analysis of insersion size of sequencing data of each sample is realized by CollectInsertSizeMetric.jar from picard software tool pack." )
pic_list=addElement(dom, root,"pic_list", name="Insertion distribution map", type="type1", desc="Note：X-coordinate shows the length of insertions, Y-coordinate shows number of its corresponding reads.")
insertImages=glob.glob(analysis_dir+'/Up_load/Mapping/*insert.cut.png')
addElementListImage(dom,pic_list,'pic',insertImages,web_dir)
addElement(dom, root,"p", type="type1", desc="As can be seen from the above figure, the insertion length distribution conforms to the normal distribution, indicating that there is no anomaly in the library construction of sequencing data." )
addElement(dom, root,"h3", name="Depth distribution statistics", type="type1", desc="Depth distribution statistics" )
addElement(dom, root,"p", type="type1", desc="After locating Reads to the reference genome, the base coverage situation on the reference genome can be calculated. The percentage of the number of bases covered by reads on the reference genome in the genome is called genome coverage; the number of reads covering on bases is called the coverage depth. Genome coverage can reflect the integrity of variation detection on the reference genome. The more regions covered, the more mutation sites that can be detected. The coverage was mainly affected by the sequencing depth and the genetic relationship between the sample and the reference genome. The coverage depth of the genome will affect the accuracy of variation detection, the accuracy of variation detection will be higher in the region with higher coverage depth (non-repetitive sequence region). In addition, if the distribution of base coverage depth on the genome is relatively uniform, it also indicates that the randomness of sequencing is good. The distribution curves of base coverage depth and coverage degree are shown in the following figure：" )
pic_list=addElement(dom, root,"pic_list", name="Sample depth distribution", type="type1", desc="Note：The above figure reflects the basic situation of sequencing depth distribution, the X-coordinate shows sequencing depth, the left Y-coordinate shows the percentage of bases corresponding to this depth, which corresponds to the red curve, the right Y-coordinate shows the percentage of bases corresponding to the below depth, which corresponds to the blue curve.")
depthImages=glob.glob(analysis_dir+'/Up_load/Mapping/*depth.dis.cut.png')
addElementListImage(dom,pic_list,'pic',depthImages,web_dir)
addElement(dom, root,"p" ,type="type1", desc="The average coverage depth of each sample and the corresponding genome coverage ratio of each depth are shown in the following table：" )
addElement(dom, root,"table" ,name="Statistics of sample coverage depth and the coverage ratio", type="full", desc="Note：BMK ID：The uniform sample ID set by BMK；Ave_depth：Average coverage depth of sample; The following three columns represent the ratio of the number of bases at or above a given coverage depth to the total number of bases of the reference genome, which are above 1X, 5X and 10X respectively." ,path="Up_load/Mapping/%s.depth_stat.xls"%Project ) 
addElement(dom, root,"p", type="type1", desc="Plot according to the coverage depth of each site of chromosome, if the coverage depth is evenly distributed on the chromosome, it can be considered that the randomness of sequencing is good. The sample chromosome coverage depth distribution map is shown below：" )
pic_list=addElement(dom, root,"pic_list", name="Distribution map of chromosome coverage depth of samples", type="type1", desc="Note：X-coordinate shows position of chromosome, Y-coordinate shows the value of the logarithm (log2) of the coverage depth of the corresponding position on the chromosome.")
covImages=glob.glob(analysis_dir+'/Up_load/Mapping/*depth.cov.txt.png')
addElementListImage(dom,pic_list,'pic',covImages,web_dir)
addElement(dom, root,"p", type="type1", desc="It can be seen from the above figure that the genome is covered  evenly, indicating good randomness of sequencing. The uneven depth may be caused by repetitive sequences and PCR bias." )

addElement(dom, root,"h2", name="Variation detection and annotation", type="type1", desc="Variation detection and annotation" )
addElement(dom, root,"h3", name="Tools and methods for variation detection", type="type1", desc="Tools and methods for variation detection" )
addElement(dom, root,"p", type="type1", desc="Detection of SNP（Single Nucleotide Polymorphism）and small InDel（small Insertion and Deletion）mainly uses GATK<a href=\"#ref2\">[2]</a> software tools package. According to the localization result of Clean Reads in the reference genome, Picard<a href=\"#ref3\">[3]</a> is used to filter abundant reads（MarkDuplicates）to ensure the accuracy of testing result. Then HaplotypeCaller（Local haplotype assembly）algorithm from GATK is used to detect SNP and InDel variations. gVCF was generated from each sample firstly, then population joint-genotype was performed. Finally, mutation sites set was collected after filtering." )
addElement(dom, root,"h3", name="Quality control of variation detection result", type="type1", desc="Quality control of variation detection result" )
addElement(dom, root,"p", type="type1", desc="Variation result was filtered strictly to ensure its reliability, the main filtering parameters are as follows：" )
addElement(dom, root,"p", type="type1", desc="（1）Use subprogram vcfutils.pl（varFilter -w 5 -W 10）from bcftools to filter SNP in the vicinity of INDEL within 5bp and SNP in the adjacent INDEL within 10bp；" )
addElement(dom, root,"p", type="type1", desc="（2）clusterSize 2  clusterWindowSize 5，Indicating that the number of variations in the 5bp window should not exceed 2；" )
addElement(dom, root,"p", type="type1", desc="（3）QUAL < 30，quality value in Phred format, indicating the possibility of variant variation in this site. Values with quality value less than 30 should be filtered out；" )
addElement(dom, root,"p", type="type1", desc="（4）QD < 2.0，The ratio of the Quality of the variation divided by the coverage depth, the coverage depth is the sum of the coverage depth of all samples at this site containing variant bases. Values with OD less than 2.0 are filtered out；" )
addElement(dom, root,"p", type="type1", desc="（5）MQ < 40，the mean square root of alignment quality value of alll reads aligned to this site. Values with MQ less than 40 are filtered out；" )
addElement(dom, root,"p", type="type1", desc="（6）FS > 60.0，value converted from the p-value of Fisher's test, it describes whether there is an obvious positive and negative chain specificity in the read that only contains variations or read that only contains reference sequence bases during sequencing or alignment. In other words, FS should be close to zero if there is no chain-specific alignment. Values with FS above 60 are filtered out；" )
addElement(dom, root,"p", type="type1", desc="（7）Other variation filtering parameters are treated by following the default values officially specified by GATK." )
addElement(dom, root,"h3", name="High quality variation result showing", type="type1", desc="High quality variation result showing" )
addElement(dom, root,"p", type="type1", desc="Variation results are shown in vcf format files, vcf file consists of comment line, title line, and data line. The comment line contains the meaning interpretation of the various identifiers used in the INFO and FORMAT columns of the data line, while the title line and data line contain the variation detection result information of each sample, the format is shown as follows：" )
addElement(dom, root,"table", name="Demo list of SNP variation result information", type="full", desc="", path="src/images/table9" )
addElement(dom, root,"p", type="type1", desc="VCF file of SNP variation result of this project can be seen from：Up_load/SNP/*.vcf." )
addElement(dom, root,"p", type="type1", desc="The meaning of each column is as follows：" )
addElement(dom, root,"table", name="Meaning of each column", type="full", desc="" ,path="src/images/table10-en" ) 
addElement(dom, root,"p", type="type1", desc="Its FORMAT（GT:AD:DP:GQ:PL）：Keywords are separated by colon." )
addElement(dom, root,"p", type="type1", desc="GT：genotype，indicating the genotype of this sample, for a diploid organism, the GT value represents types of two alleles that the sample carries at this locus, 0 means that it's the same as REF；1 means that it's the same as ALT，0/0 means homozygous and it is consistent with REF；0/1 means heterozygosity，two alleles are ALT and REF；1/1 means homozygous and both are ALT；" )
addElement(dom, root,"p", type="type1", desc="AD：allele depth, corresponds to two comma-separated values, these two values mean numbers of reads that are covered to REF and ALT, equivalent to the sequencing depth that supports REF and ALT；" )
addElement(dom, root,"p", type="type1", desc="DP：depth of coverage, the total number of reads covered to this site, equals to the depth of this site；" )
addElement(dom, root,"p", type="type1", desc="GQ：Quality of the assigned genotype, represents the mass value of the most likely genotype；" )
addElement(dom, root,"p", type="type1", desc="PL：Normalized Phred-scaled likelihoods of the possible genotypes, they correspond to three comma-separated values, these three values represent standardized phred-scaled likelihood value (L) without prior test of this site when its genotype is 0/0, 0/1 and 1/1. If they are converted to possibility (P) that supports this genotype, since that L=-10lgP，so P=10^（-L/10）, when L equals to 0, P=10^0=1. So the smaller the value is, the greater the probability of support is, the greater the likelihood of this kind of genotype." )
addElement(dom, root,"p", type="type1", desc="Please find the web for the detailed description of vcf file：<a href=\"http://gatkforums.broadinstitute.org/discussion/1268/how-should-i-interpret-vcf-files-produced-by-the-gatk\" target=\"_blank、\">http://gatkforums.broadinstitute.org/discussion/1268/how-should-i-interpret-vcf-files-produced-by-the-gatk</a>" )
addElement(dom, root,"p", type="type1" ,desc="To ensure the reliability of SNP, collect the accumulative distribution of support number of reads of detected SNP and the distance between adjacent SNPs." )
addElement(dom, root,"pic", name="SNP quality distribution map", type="type1" ,desc="Note：The left is the cumulative map of support number of SNP reads, the right is the cumulative map of the distance between adjacent SNPs. " ,path="Up_load/SNP/SNP.quality.distribution.png" )
addElement(dom, root,"p", type="type1", desc="SNP mutation can be divided into conversion and transference, mutations between the same type of bases are called transitions, such as those between purines and purines and between pyrimidines and pyrimidines, while mutations between different types of bases are called transversions, such as those between purines and pyrimidines. Generally conversion is more likely to occur than transition, so the ratio of transition/transition (Ti/Tv) is generally greater than 1, the specific value is related to the species tested. For diploid or polyploid species, if a certain SNP site on homologous chromosomes is the same base, the SNP site is called homozygous SNP site; If the SNP site on homologous chromosomes contains different types of bases, the SNP site is called heterozygous SNP site. The higher the number of homozygous SNP sites, the greater the difference between the sample and the reference genome, the higher the number of heterozygous SNP sites, the higher the heterozygosity of the sample. The specific results are related to the material selection of the sample. The SNP detection results between the sample and the reference genome are shown below." )
addElement(dom, root,"table", name="Statistics of detected SNP", type="full", desc="", path="Up_load/SNP/%s.raw.filter.snp.snp.stat"%Project ) 
addElement(dom, root,"p", type="type1", desc="Description of each column is shown as the follow table：" )
addElement(dom, root,"table", name="Meaning of each column", type="full", desc="", path="src/images/table12-en" ) 
if sampleNum ==1:
    addElement(dom, root,"h3", name="Sample SNP detection (single sample)", type="type1", desc="Sample SNP detection (single sample)" )
    addElement(dom, root,"p", type="type1" ,desc="According to the alignment result between sample and the reference genome, the variation sites of the sample and the reference genome are summarized. The format of the SNP list file of the sample is as follows：" )
    addElement(dom, root,"table", name="Sample SNP list" ,type="type1", desc="", path="Up_load/SNP/%s.raw.filter.snp.snp"%Project )
    addElement(dom, root,"p", type="type1" ,desc="Please find the detailed list data here：Up_load/SNP/%s.raw.filter.snp.snp."%Project )
else:
    addElement(dom, root,"h3", name="SNP detection between samples (multiple samples)", type="type1", desc="SNP detection between samples (multiple samples)" )
    addElement(dom, root,"p", type="type1" ,desc="According to the alignment result between sample and the reference genome, the variation sites of all samples are summarized. The format of the SNP list file of all samples is as follows：" )
    addElement(dom, root,"table", name="Samples SNP list" ,type="type1", desc="", path="Up_load/SNP/%s.raw.filter.snp.DEG.snp"%Project )
    addElement(dom, root,"p", type="type1" ,desc="Please find the detailed list data here：Up_load/SNP/%s.raw.filter.snp.DEG.snp." %Project)
addElement(dom, root,"p", type="type1" ,desc="SNP genotypes are encoded using standard nucleotide symbols, as shown in the symbol table below：" )
addElement(dom, root,"table", name="SNP genotypes code", type="full", desc="", path="src/images/table14-en" )
addElement(dom, root,"p", type="type1" ,desc="Genome-wide SNP mutations can be divided into six classes. Use T:A>C:G as an example, this type of SNP mutation has T>C and A>G. Since the sequencing data can be aligned to the positive chain or the negative chain of the reference genome, when T>C type mutation occurs at the positive chain of the reference genome, A>G type mutation will occur at the same position of the genome negative chain, so T>C and A>G are classified as one class." )
addElement(dom, root,"pic", name="SNP", type="type1" ,desc="Note：Y-coordinate shows number of SNP, X-coordinate shows SNP mutation type. " ,path="Up_load/SNP/SNP.mutation.distribution.png" )
if sampleNum == 1:
	pass
elif os.path.exists(analysis_dir+"/Up_load/SNP/SNP.venn.png"):
    addElement(dom, root,"p", type="type1", desc="The statistical results of SNP among samples are shown in the following figure：" )
    addElement(dom, root,"pic", name="Statistical Venn diagam of SNP among samples", type="type1" ,desc="Note：Venn statistics of the number of mutation sites only consider whether the locations are the same, not whether the genotypes are the same." ,path="Up_load/SNP/SNP.venn.png" )
else:
    addElement(dom, root,"p", type="type1", desc="The statistical results of SNP among samples are shown in the following table：" )
    addElement(dom, root,"table", name="Statistical table of SNP among samples", type="full" ,desc="BMK ID：The uniform sample ID set by BMK；The values in the table are the number of SNPS between the corresponding horizontal and vertical samples." ,path="Up_load/SNP/DEG_snp.stat" )

addElement(dom, root,"h3", name="SNP result annotation", type="type1" ,desc="SNP result annotation" )
addElement(dom, root,"p", type="type1", desc="SnpEff<a href=\"#ref4\">[4]</a> is a software used to annotate the variation (SNP, Small InDel) and predict the effect of the variation. According to the location of the mutation locus on the reference genome and the gene location information on the reference genome, the region where the mutation locus occurs in the genome (intergene region, gene region or CDS region, etc.) and the impact of the mutation (synonymous mutation, non-synonymous mutation, etc.) can be obtained. The software can use VCF format files as input and output. The output result will add the following fields to the INFO column of the VCF file：EFF= Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_Length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype_Number [ | ERRORS | WARNINGS ] )" )

addElement(dom, root,"p" ,type="type1", desc="Each identifier is described as follows：" )
addElement(dom, root,"table",name="Identifier description", type="full",desc="Note：If the above results cannot be obtained, the corresponding column is null. Find F file of SnpEff for the details：<a href=\"http://snpeff.sourceforge.net/SnpEff_manual.html#output\" target=\"_blank\">http://snpeff.sourceforge.net/SnpEff_manual.html#output.</a> ", path="src/images/table16-en" ) 
pic_list=addElement(dom, root,"pic_list", name="SNP annotation result statistics", type="type1", desc="")
snpAnnImages=glob.glob(analysis_dir+'/Up_load/SNP/*pie.png')
addElementListImage(dom,pic_list,'pic',snpAnnImages,web_dir)
addElement(dom, root,"p", type="type1", desc="The meanings of the above symbols are shown in the following table：" )
addElement(dom, root,"table",name="Description of meaning of functional region",type="full" ,desc="",path="src/images/table18-en" )
addElement(dom, root,"h3", name="Small InDel detection between sample and the reference genome", type="type1", desc="mall InDel detection between sample and the reference genome" )
addElement(dom, root,"p",type="type1", desc="According to the localization results of sample Clean Reads on the reference genome, the insertion and deletion of Small fragments (Small InDel) between the sample and the reference genome were detected. GATK was used to detect InDel of the sample. Small InDel variation is generally less than SNP variation, which also reflects the difference between the sample and the reference genome, InDel in the coding region will cause code shift mutation, leading to changes in gene function. Please find the table below" )

addElement(dom, root,"table", name="InDel statistics of whole genome and the coding region",type="full", path="Up_load/INDEL/%s.combine.Indel.stat"%Project,desc="Note：CDS：coding region InDel statistics；Genome：Genome-wide InDel statistics；Insertion：Number of detected insertions；Deletion：Number of detected deletions；Het：Number of heterozygous InDel；Homo：Number of homozygous InDel；Total：Total number of detected InDel (Duplicates removed)."  ) 

addElement(dom, root,"p", type="type1", desc="Collect statistics on basis of the InDel length of the sample in the CDS region and the genome-wide region, the length distribution is shown in the figure below：" )
addElement(dom, root,"pic", name="Distribution map of InDel of whole genome and coding region" ,type="type1", desc="Note：Y-coordinate shows the length of InDel (within 10bp)，ordinate larger than 0 means Insertion，ordinate less than 0 means Deletion，X-coordinate shows the corresponding number.", path="Up_load/INDEL/all.sample.indel.length.distribution.png" ) 

addElement(dom, root,"h3", name="Small InDel detection between samples", type="type1", desc="Small InDel detection between samples" )
addElement(dom, root,'p', type="type1" ,desc="According to the detection results of Small InDel between sample and reference genome, the alignment results of sample sequencing data are shown in the following table." )
addElement(dom, root,"table", name="Statistics of Small InDel sequencing data of samples", type="type1", desc="", path="Up_load/INDEL/%s.raw.filter.indel.anno.gatk.list"%Project )
addElement(dom, root,"p", type="type1", desc="Please find the details of sample Small indel sequencing data here：Up_load/INDEL/%s.raw.filter.indel.anno.gatk.list."%Project)
addElement(dom, root,"p", type="type1", desc="Description of meaning of each column is as follows：" )
addElement(dom, root,"table", name="Description of meaning of each column", type="full", desc="", path="src/images/table21-en" ) 

if sampleNum == 1:
	pass
elif os.path.exists(analysis_dir+"/Up_load/INDEL/INDEL.venn.png"):
    addElement(dom, root,"p", type="type1", desc="The statistical result of INDEL between samples is shown in the following figure：" )
    addElement(dom, root,"pic", name="Statistical Venn diagam of INDEL among samples", type="type1" ,desc="Note：Venn statistics of the number of mutation sites only consider whether the locations are the same (start position of indel), not whether the genotypes are the same." ,path="Up_load/INDEL/INDEL.venn.png" )
else:
    addElement(dom, root,"p", type="type1", desc="The statistical results of Small InDel detection between two samples are shown in the following table：" )
    addElement(dom, root,"table", name="Statistical result of InDel between samples", type="full" ,desc="Note：BMK ID：The uniform sample ID set by BMK；Each value in the table is the number of Small InDel between the corresponding horizontal and vertical samples." ,path="Up_load/INDEL/DEG_indel.stat" )

addElement(dom, root,"h3", name="Small InDel annotation", type="type1", desc="Small InDel annotation" )
addElement(dom, root,"p", type="type1", desc="According to the location information of detected Small InDel locus on the reference genome, make alignment of information such as genes in the reference genome and CDS location, etc. (generally can be found in gff files) to annotate whether the InDel locus occurs in the intergene region, gene region or CDS region, and whether it is a code shift mutation. The annotation of Small InDel is realized by SnpEff software. InDel with code shift mutation may result in change of gene function, the detailed annotation result is shown in the figure below：" )

pic_list=addElement(dom, root,"pic_list", name="InDel annotation result statistics", type="type1", desc="")
snpAnnImages=glob.glob(analysis_dir+'/Up_load/INDEL/*pie.png')
addElementListImage(dom,pic_list,'pic',snpAnnImages,web_dir)
addElement(dom, root,"p", type="type1", desc="Description of labels above is shown in the following table：" )
addElement(dom, root,"table", name="Description of each functional region", type="full", desc="", path="src/images/table24-en" ) 

reference = 5
#SV detection and annotation
if(os.path.exists(analysis_dir+'/Up_load/SV')):
    addElement(dom, root,"h2", name="SV detection and annotation", type="type1", desc="SV detection and annotation" )
    addElement(dom, root,"h3", name="SV detection", type="type1", desc="SV detection" )
    addElement(dom, root,"p", type="type1", desc="Genomic structural variation (SV) refers to the variation of chromosome structure of species, such as large fragment insertion, deletion, inversion and translocation. breakDancer is commonly used to detect structural variations, it first obtained the insertion fragment size and variance of the sequencing library based on the alignment results between the sequence and the reference genome, and then found for possible structural variations by querying abnormal alignment results between the sequence and the reference genome (deviation of the insertion fragment, inconsistent alignment direction, etc.). The SV result file contains the title line and the data line, as shown below：" )
    addElement(dom, root,"table", name="Demo file of SV result (Sample %s)"%Sample_Project,type="type1", desc=" ",path="Up_load/SV/%s.%s.%s.dedup.realn.bam.max"%(Project,Project,Sample_Project))
    addElement(dom, root,"p", type="type1", desc="（Sample R*）SV result file can be found here：Up_load/SV/%s.%s.R*.dedup.realn.bam.max"%(Project,Project) )
    addElement(dom, root,"p", type="type1", desc="Description of meaning of each column in the title column is shown in the following table：" )
    addElement(dom, root,"table", name="Description of meaning of each column", type="full" ,desc="Note：Description of SV result file can be found here：<a href=\"https://github.com/kenchen/breakdancer#readme\" target=\"_blank\">https://github.com/kenchen/breakdancer#readme</a>" ,path="src/images/table26-en" ) 
    addElement(dom, root,"p", type="type1", desc="Use BreakDancer<a href=\"#ref%d\">[%d]</a> software to detect insertion (INS), Deletion (DEL), Inversion (INV), Intra-chromosomal Translocation (ITX) and Inter-chromosomal Translocation (CTX) of  tested samples and the reference genome on the basis of the relationship of Pair-end reads aligned to the reference genome and the actual Insert Size. Statistics of number of each type of SV of %s samples is shown in the following table：" %(reference,reference,sampleNum))
    addElement(dom, root,"table", name="SV number statistical table", type="full", desc="Note：BMK ID：The uniform sample ID set by BMK；SV: Total number of structural variations；INS：number of insertion variations；DEL：number of deletion variations；INV：number of invertion variations；ITX：number of Intra-chromosomal Translocation variations；CTX：number of Inter-chromosomal Translocation variations；UN：Complex structural variations.", path="Up_load/SV/%s.sv.stat.xls"%Project ) 
    addElement(dom, root,"h3", name="SV annotation", type="type1", desc="SV annotation" )
    addElement(dom, root,"p", type="type1", desc="According to the location information of detected SV of sample on the reference genome, make alignment of information such as genes in the reference genome and CDS location, etc. (generally can be found in gff files) to annotate whether the SV occurs in the intergene region, gene region or CDS region, etc. The three types of structural variations of deletion (DEL), insertion (INS) and inversion (INV) were annotated, the statistical results are shown in the table：" )
    addElement(dom, root,"table", name="Statistical table of SV annotation result", type="full", desc="Note：BMK ID：The uniform sample ID set by BMK；Type：SV type；Exon：Variations in exon region；Intron：Variations in intron region；Intergenic：Variations in intergenic region.", path="Up_load/SV/%s.sv.anno.stat"%Project ) 
    reference = reference + 1
    
#CNV detection and annotation
if(os.path.exists(analysis_dir+'/Up_load/CNV')):
    addElement(dom, root,"h2", name="CNV detection", type="type1", desc="CNV detection" )
    addElement(dom, root,"p", type="type1",desc="Use FREEC<a href=\"#ref%d\">[%d]</a> to detect CNV by the depth distribution of sample sequencing reads on the reference genome, then draw distribution of Copy Number Gain and Loss on the reference genome. Part of result file of FREEC of sample %s：" %(reference,reference,Sample_Project))
    addElement(dom, root,"table", name="FREEC result file", type="type1", desc="Note：chr：chromosome；strat：the start position of CNV ；end：the end position of CNV；predicted copy number：predicted copy number；type of alteration：CNV type.",path="Up_load/CNV/%s.%s.dedup.realn.bam_CNVs"%(Project,Sample_Project))
    addElement(dom, root,"p", type="type1", desc="（Sample R*）FREEC result file can be found here：Up_load/CNV/%s.R*.dedup.realn.bam_CNVs." %(Project) )
    reference = reference + 1
#circos plot
if(os.path.exists(analysis_dir+'/Up_load/Circos')):
    addElement(dom, root,"h2", name="Distribution of each type of variation on genome", type="type1", desc="Distribution of each type of variation on genome" )
    addElement(dom, root,"p", type="type1", desc="The distribution of the detected results of each type of variation was shown in the circos plot, which was made by circos software, website of the software：<a href=\"http://circos.ca/\" target=\"_blank\">http://circos.ca/</a>. The distribution of each type of variation on the chromosome of each sample is shown in the figure below：" )
    if os.path.exists(args.od + 'Up_load/SV') and os.path.exists(args.od + 'Up_load/CNV'):
    	pic_list=addElement(dom, root,"pic_list", name="The distribution of each type of variation in chromosome", type="type1", desc="Note：From the outside to the inside are: chromosome coordinates, SNP density distribution, InDel density distribution, CNV, SV (INS, DEL, INV, ITX(red line), CTX(green line)) distribution on the genome (unit M).")
    else:
    	pic_list=addElement(dom, root,"pic_list", name="The distribution of each type of variation in chromosome", type="type1", desc="Note：From the outside to the inside are: chromosome coordinates, SNP density distribution, InDel density distribution.")
    circosImages=glob.glob(analysis_dir+'/Up_load/Circos/*.circos.png')
    addElementListImage(dom,pic_list,'pic',circosImages,web_dir)

#Differentially expressed gene analysis
if(os.path.exists(analysis_dir+'/Up_load/Diff_analysis')):
    addElement(dom, root,"h2", name="mutation genes analysis on DNA level", type="type1", desc="mutation genes analysis on DNA level" )
    addElement(dom, root,"h3", name="mutation genes mining on DNA level", type="type1", desc="mutation genes mining on DNA level" )
    addElement(dom, root,"p", type="type1", desc="Mutations occurring in the CDS region may cause change of gene function, by finding non-synonymous mutation SNP, InDel and SV genes occurring in the CDS region between the reference genome and the sample, genes with possible functional differences between the sample and the reference genome can be found. The variation between each sample and the reference genome is shown in the following table：" )
    addElement(dom, root,"table", name=" Taxonomic statistics of differential genes caused by various mutations", type="full", desc="Note：BMK_ID：The uniform sample ID set by BMK；Genes with Non-synonymous：The number of genes with non-synonymous mutations, the case that one gene with multiple non-synonymous mutations were not counted again；Genes with InDel：Number of genes with Small InDel；Genes with SV：Number of genes with SV." ,path="Up_load/Diff_analysis/Diff_gene.stat" ) 
    addElement(dom, root,"h3", name="DNA level function annotation of mutation genes", type="type1", desc="DNA level function annotation of mutation genes" )
    ref1 = reference
    ref2 = ref1 + 1
    ref3 = ref2 + 1
    ref4 = ref3 + 1
    ref5 = ref4 + 1
    ref6 = ref5 + 1
    addElement(dom, root,"p", type="type1" ,desc="Use BLAST<a href=\"#ref%d\">[%d]</a> to compare mutation genes with function databases such as NR<a href=\"#ref%d\">[%d]</a>，SwissProt<a href=\"#ref%d\">[%d]</a>，GO<a href=\"#ref%d\">[%d]</a>，COG<a href=\"#ref%d\">[%d]</a>，KEGG<a href=\"#ref%d\">[%d]</a>, annotation of these genes can be obtained to analyse gene function." %(ref1,ref1,ref2,ref2,ref3,ref3,ref4,ref4,ref5,ref5,ref6,ref6))
    addElement(dom, root,"p" ,type="type1" ,desc="NR database：It is the non-redundant Protein database from NCBI database, including protein databases like SwissProt、PIR(Protein Information Resource)、PRF(Protein Research Foundation)、PDB(Protein Data Bank) and Protein Data translated from GenBank and RefSeq CDS Data." )
    addElement(dom, root,"p" ,type="type1" ,desc="SwissProt database：It's an annotated protein sequence database maintained by the European institute of bioinformatics (EBI). This database consists of protein sequence entrys, each entry contains protein sequence, citation information, taxonomic information, annotations, etc., the annotation includes protein function, post-transcriptional modification, special sites and regions, secondary structure, quaternary structure, similarity with other sequences, relationship between sequence disability and disease, sequence variants and conflicts, etc. Redundant sequences are minimized in SwissProt, the cross-references are established with more than 30 other data, including nucleic acid sequence library, protein sequence library and protein structure library." )
    addElement(dom, root,"p" ,type="type1" ,desc="GO database：Gene Ontology（GO）is an internationally standardized classification system of Gene function that provides a dynamically updated set of controlled vocabulary to comprehensively describe the properties of genes and Gene products in organisms. The most basic concept in a Gene Ontology is term. Each entry in GO has an unique number tag, like GO:nnnnnnn, and a term ID, like“cell”、“fibroblast growth factor receptor binding”or“signal transduction”. Each term belongs to an ontology, there are three ontologies, which are Molecular function, Cellular component and Biological process." )
    addElement(dom, root,"p" ,type="type1" ,desc="COG database：It is a database for direct homology classification of gene products. Each COG protein is assumed to be from an ancestor protein and is divided into orthologs and paralogs. Orthologs are proteins from different species that have evolved from vertical families and typically retain the same functions as the original proteins. Paralogs are proteins derived from gene replication in a given species that may evolve new functions related to the original." )
    addElement(dom, root,"p" ,type="type1" ,desc="KEGG database：KEGG integrates the current knowledge of chemical compounds, reactions and molecular action networks in biochemistry, it is a database for the systematic analysis of the metabolic pathways of gene products in cells and the functions of these gene products, which is benefitial to the study of gene and expression information as a whole network." )
    addElement(dom, root,"p" ,type="type1" ,desc="The annotation list file of mutation genes is shown in the table below：" )
    addElement(dom, root,"table", name="Annotation list of mutation gene of sample %s"%Sample_Project, type="type1", desc=" ", path="Up_load/Diff_analysis/%s/Integrated_Function.annotation.xls"%Sample_Project)
    addElement(dom, root,"p", type="type1", desc="(Sample R*) Mutation gene annotation list can be found here:Up_load/Diff_analysis/R*/Integrated_Function.annotation.xls" )
    addElement(dom, root,"p", type="type1", desc="The specific statistical result of each sample is shown in the following table:" )
    addElement(dom, root,"table", name="Statistical table of mutation genes annotation", type="full", desc="Note：Sample：sample name；Gene_Number：number of DEGs annotated to the database.", path="Up_load/Diff_analysis/Allsample_Diffgene.list" ) 
    addElement(dom, root,"p", type="type1", desc="Statistical result of GO classification of mutation genes is shown in below figure：" )
    addElement(dom, root,"pic", name="GO annotation clustering of mutation genes of sample %s"%Sample_Project, type="type1",path="Up_load/Diff_analysis/%s/go_enrichment/%s.GO.png" %(Sample_Project,Project),desc="Note：X-coordinate shows each classification of GO, the left Y-coordinate shows the percentage of genes number, the left shows genes number") 
    addElement(dom, root,"p", type="type1" ,desc="statistical result of COG classification of mutation genes is shown in below figure：" )
    addElement(dom, root,"pic" ,name="COG annotation classification map of mutation genes of sample %s"%Sample_Project,type="type1" ,path="Up_load/Diff_analysis/%s/Cog_Anno/%s.Cog.classfy.png"%(Sample_Project,Project) ,desc="Note：X-coordinate shows each classification of COG, Y-coordinate shows number of genes. In different functional classes, the proportion of genes reflects the metabolic or physiological bias in the corresponding period and environment, which can be explained scientifically in combination with the distribution of research objects in each functional class." ) 
    snpAnnImages=glob.glob(analysis_dir+'/Up_load/Diff_analysis/%s/pathway/kegg_map/*.png'%Sample_Project)
    if snpAnnImages != []:
	    addElement(dom, root,"p", type="type1", desc="Result of metabolic pathways of mutation genes is shown in the figure below：" )
	    pic_list=addElement(dom, root,"pic_list", name="Metabolic pathways map of mutation genes of sample %s"%Sample_Project, type="type1", desc="Note：Numbers in box stand for enzyme numbers, suggests that the corresponding gene is associated with the enzyme. The whole pathway is formed by a variety of different enzymes through complex biochemical reactions, all mutation genes related to this pathway are shown in red boxes. Researchers can focus on the genes of related metabolic pathways according to their research objects, and explain the origin of the corresponding genes through metabolism.")
	    addElementListImage(dom,pic_list,'pic',snpAnnImages,web_dir)
	    
addElement(dom, root,"h1", name="Description to view result files", type="type1", desc="Description to view result files" )
addElement(dom, root,"p", type="type1", desc="(1) There is a readme.txt description in the upload directory that detailedly introduces what each file represents. The uploaded result data files are mainly in text format (fa file, TXT file, detail file, XLS file, etc.). For viewing files on Windows, it is recommended to use Editplus or UltraEdit as text browsers, otherwise the files will be too large and cause a crash, or to open the large text files using PilotEdit Lite（<a href=\"http://www.pilotedit.com/\" target=\"_blank\">http://www.pilotedit.com/</a>）.On a Unix or Linux system, large text files can be viewed with operating commands such as Less, etc." )
addElement(dom, root,"p" ,type="type1" ,desc="(2) The report file contains an image file in SVG format, a vectorized image file that can be enlarged without distortion. To view files in SVG format, please install the SVG plug-in first." )

#References
reference = 5
ref_list=addElement(dom, root,"ref_list", name="References", type="type1", desc="")
addElement(dom, ref_list,"ref", id="1" ,name="Li H,Durbin R.Fast and accurate short read alignment with Burrows-Wheeler Transform.Bioinformatics, 2009 25:1754-60", link="http://bioinformatics.oxfordjournals.org/content/25/14/1754.long" )
addElement(dom, ref_list,"ref" ,id="2" ,name="McKenna A, Hanna M, Banks E, Sivachenko A, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010 20:1297-303" ,link="http://genome.cshlp.org/content/20/9/1297.long" )
addElement(dom, ref_list,"ref" ,id="3" ,name=" Picard: http://sourceforge.net/projects/picard/ .(Picard)" ,link="http://sourceforge.net/projects/picard/" )
addElement(dom, ref_list,"ref" ,id="4" ,name="Cingolani P, Platts A, Wang le L, et al. A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3., Fly (Austin). 2012 Apr-Jun;6(2):80-92.",link="http://www.tandfonline.com/doi/full/10.4161/fly.19695" )

if(os.path.exists(analysis_dir+'/Up_load/SV')):
    addElement(dom, ref_list,"ref" ,id="%d" %(reference),name="Ken Chen, John W Wallis, Michael D McLellan, et al. BreakDancer: an algorithm for high-resolution mapping of genomic structural variation. Nature Methods, 2009 6:677-681" ,link="http://www.nature.com/doifinder/10.1038/nmeth.1363" )
    reference = reference + 1
    
if(os.path.exists(analysis_dir+'/Up_load/CNV')):
    addElement(dom, ref_list,"ref" ,id="%d" %(reference) ,name="Boeva V, Popova T, Bleakley K, Chiche P, Cappo J, Schleiermacher G, Janoueix-Lerosey I, Delattre O, Barillot E. (2012) Control-FREEC: a tool for assessing copy number and allelic content using next generation sequencing data. Bioinformatics. Bioinformatics, 2012, 28(3):423-5. PubMed PMID: 22155870. " ,link="http://bioinformatics.oxfordjournals.org/content/28/3/423.long")
    reference = reference + 1
    
if(os.path.exists(analysis_dir+'/Up_load/Diff_analysis')):
    addElement(dom, ref_list,"ref" ,id="%d" %(reference),name="Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch?ffer, et al. Gapped BLAST and PSI-BLAST: A New Generation of Protein Database Search Programs. Nucleic Acids Res, 1997 25(17): 3389" ,link="http://nar.oxfordjournals.org/content/25/17/3389.long" )
    reference = reference + 1
    addElement(dom, ref_list,"ref",id="%d" %(reference),name="Yangyang DENG, Jianqi LI, Songfeng WU, et al. Integrated nr Database in Protein Annotation System and Its Localization. Computer Engineering, 2006 32(5):71_74" ,link="http://scholar.google.com.hk/schhp?hl=zh-CN" )
    reference = reference + 1
    addElement(dom, ref_list,"ref", id="%d" %(reference),name="Yangyang DENG, Jianqi LI, Songfeng WU, et al. Integrated nr Database in Protein Annotation System and Its Localization. Computer Engineering, 2006 32(5):71_74" ,link="http://scholar.google.com.hk/schhp?hl=zh-CN" )
    reference = reference + 1
    addElement(dom, ref_list,"ref" ,id="%d"%(reference) ,name="Michael Ashburner, Catherine A. Ball, Judith A. Blake, et al. Gene ontology:  tool for the unification of biology. The Gene Ontology Consortium. Nat Genet, (25): 25_29",link="http://www.nature.com/doifinder/10.1038/75556" )
    reference = reference + 1
    addElement(dom, ref_list,"ref" ,id="%d"%(reference) ,name="Roman L.Tatusov, Michael Y.Galperin, Darren A.Natale, et al. The COG database: a tool for genome_scale analysis of protein functions and evolution. Nucleic Acids Res, 2000_7_1 28(1):33_6" ,link="http://nar.oxfordjournals.org/content/28/1/33.long" )
    reference = reference + 1
    addElement(dom, ref_list,"ref" ,id="%d" %(reference),name="Minoru Kanehisa, Susumu Goto, Shuichi Kawashima, et al. The KEGG resource for deciphering the genome. Nucleic Acids Res 2004, (32):D277_D280" ,link="http://nar.oxfordjournals.org/content/32/suppl_1/D277.long" )
    reference = reference + 1


#################################################output my xml object####################################################
domcopy = dom.cloneNode(True)
Indent(domcopy, domcopy.documentElement)
f = open(web_dir+'/'+args.name+'.xml', 'wb')
writer = codecs.lookup('utf-8')[3](f)
domcopy.writexml(writer, encoding = 'utf-8')
domcopy.unlink()
f.close()
sys.stdout.write("%s:complete xml file: %s \n"%(timenow(),web_dir+'/'+args.name+'.xml'))
sys.stdout.write('[cmd]: /share/nas2/genome/biosoft/Python/2.7.8/bin/python ' + Bin+'/htmlConvert/xml2HtmlConverter.py  -i ' +web_dir+'/'+args.name+'.xml' +' -o '+web_dir+" -l 5 -t "+Bin+'/src'+'\n')
os.system('/share/nas2/genome/biosoft/Python/2.7.8/bin/python  '+Bin+'/htmlConvert/xml2HtmlConverter.py  -i ' +web_dir+'/'+args.name+'.xml' +' -o '+web_dir+" -l 5 -t "+Bin+'/src ')

if not os.path.exists(analysis_dir + "/Web_Report"):
	os.mkdir(analysis_dir + "/Web_Report")
else:
	os.system("rm -r " + analysis_dir + "/Web_Report")
	os.mkdir(analysis_dir + "/Web_Report")
		
os.system('mv '+analysis_dir+'/index.html '+analysis_dir+'/Web_Report')
os.system('mv '+analysis_dir+'/src '+analysis_dir+'/Web_Report')
os.system('mv '+analysis_dir+'/refinfo.txt '+analysis_dir+'/Web_Report')
os.system('mv '+analysis_dir+'/sampleinfo.txt '+analysis_dir+'/Web_Report')
os.system('mv '+analysis_dir+'/configtest.xml '+analysis_dir + '/Web_Report')
os.system('mv '+analysis_dir+'/Up_load '+analysis_dir + '/Web_Report')

if os.path.exists(analysis_dir + '/Web_Report/index.html'):
        if(delete == 'Y'):
                os.system('rm -r ' + analysis_dir + '/Analysis')
                os.system('rm -r ' + analysis_dir + '/DataAssess')
                os.system('rm -r ' + analysis_dir + '/Diff_analysis')




