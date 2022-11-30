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
    if i.startswith('platform'):Sequencing_platform=re.split(u'\s+', i)[1]
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
f.write("BMK编号\t"+"客户编号\n")
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
    addElement(dom,root,'report_abstract',value="<p class= \"p-abstract\" >分析内容:</p><p class=\" p-abstract\" >完成%s个样品的重测序，具体分析内容如下，数据评估：测序数据量，测序数据质量和GC含量的统计。与基因组比对：比对率，基因组覆盖度，基因组覆盖深度统计。变异检测和注释：SNP、InDel、SV、CNV的检测和注释。DNA水平变异基因分析：检测编码区发生SNP非同义突变、InDel突变、SV突变、CNV突变的基因。基因注释：对DNA水平变异基因进行KEGG、GO、COG、NR、SwissProt数据库注释。</p><p class=\" p-abstract\" >结果概述：</p><p class=\" p-abstract\" >本次分析数据量为%.2fGbp的Clean Data，Q30达到%.2f%%。样品与参考基因组平均比对率为%.2f%%，平均覆盖深度为%dX，基因组覆盖度为%.2f%%（至少一个碱基覆盖）。具体变异检测结果(SNP，Indel)见下面结题报告。样品与参考基因组之间检测发生SNP非同义突变、InDel突变、SV、CNV基因，并对DNA水平变异的基因进行了KEGG、GO、COG、NR、SwissProt等数据库注释(具体以合同签订为准)。</p>"%(sampleNum,cleanData,averageQ30,averagemapRate,averageDepth,coverageRate) )
else:
    addElement(dom,root,'report_abstract',value="<p class= \"p-abstract\" >分析内容:</p><p class=\" p-abstract\" >完成%s个样品的重测序，具体分析内容如下，数据评估：测序数据量，测序数据质量和GC含量的统计。与基因组比对：比对率，基因组覆盖度，基因组覆盖深度统计。变异检测和注释：SNP、InDel的检测和注释。</p><p class=\" p-abstract\" >结果概述：</p><p class=\" p-abstract\" >本次分析数据量为%.2fGbp的Clean Data，Q30达到%.2f%%。样品与参考基因组平均比对率为%.2f%%，平均覆盖深度为%dX，基因组覆盖度为%.2f%%（至少一个碱基覆盖）。具体变异检测结果(SNP，Indel)见下面结题报告。</p>"%(sampleNum,cleanData,averageQ30,averagemapRate,averageDepth,coverageRate) )
    
addElement(dom,root,'h1',name="项目基本信息" ,type="type1" ,desc="项目基本信息")
addElement(dom,root,'p',type="type1" ,desc="(1)样品信息对应编号:")
addElement(dom,root,'table',name="样品信息" ,type="full" ,desc="注：BMK编号：百迈客对样品的统一编号，实验建库和后续信息分析均使用该编号。" ,path="sampleinfo.txt")
addElement(dom,root,'p',type="type1", desc="(2) 参考基因组信息：")
addElement(dom,root,'p',type="type1" ,desc="测序物种名称：%s；下载网址：<a href=\"%s\" title=\"click\" class=\"mylink\" target=\"_blank\" >%s</a>"%(Project,ref_url,ref_url) )
addElement(dom, root,'table',name="基因组信息" ,type="full", desc="", path="refinfo.txt")
addElement(dom, root,'h1',  name="项目流程", type="type1", desc="项目流程")
addElement(dom, root,'h2' ,name="实验流程" ,type="type1" ,desc="实验流程")
if Sequencing_platform == 'Illumina':
	addElement(dom, root,"p" ,type="type1" ,desc="实验流程按照Illumina公司提供的标准protocol执行，包括样品质量检测、文库构建、文库质量检测和文库测序等流程，具体流程图如下。" )
	addElement(dom, root,"pic" ,name="实验流程图" ,type="type1", desc=" ", path="src/images/experiment_process.png" )
	addElement(dom, root,"p" ,type="type1" ,desc="样品基因组DNA检测合格后，用机械打断的方法（超声波）将DNA片段化，然后对片段化的DNA进行片段纯化、末端修复、3′端加A、连接测序接头，再用琼脂糖凝胶电泳进行片段大小选择，进行PCR扩增形成测序文库，建好的文库先进行文库质检，质检合格的文库用%s进行测序。"%Sequencing_platform)
else:
	addElement(dom, root,"p" ,type="type1" ,desc="实验流程按照华大公司提供的标准protocol执行，包括样品质量检测、文库构建、文库质量检测和文库测序等流程，具体流程图如下。" )
	addElement(dom, root,"pic" ,name="实验流程图" ,type="type1", desc=" ", path="src/images/experiment_process.BGI.png" )
	addElement(dom, root,"p" ,type="type1" ,desc="样品基因组DNA检测合格后，用机械打断的方法（超声波）将DNA片段化，然后进行片段大小选择，再对片段化的DNA末端修复、3′端加A、连接测序接头，再进行PCR扩增并纯化，将PCR扩增的DNA片段进行单链环化，形成环状DNA，以环化DNA为模板，通过独特的线性扩增模式（滚环扩增）将文库制备成DNA纳米球（DNA nanoball，DNB），制备好的DNB文库与测序芯片上的阵列式位点进行结合，从而进行文库测序。" )
addElement(dom, root,"h2" ,name="信息分析流程" ,type="type1", desc="信息分析流程" )
addElement(dom, root,"p", type="type1", desc="对测序得到的原始reads(双端序列)进行质量评估并过滤得到Clean Reads，用于后续生物信息学的分析。将Clean Reads与参考基因组序列进行比对，基于比对结果进行SNP、InDel等变异检测和注释，并实现DNA水平差异基因挖掘和差异基因功能注释等。" )
addElement(dom, root,"p", type="type1", desc="重测序生物信息分析流程见下图：" )
if os.path.exists(analysis_dir+'/Up_load/SV') and os.path.exists(analysis_dir + '/Up_load/CNV') and os.path.exists(analysis_dir + '/Up_load/Diff_analysis'):
    addElement(dom, root,"pic", name="重测序生物信息分析流程图", type="type1", desc=" ", path="src/images/Analysis_process.png" )
else:
    addElement(dom, root,"pic", name="重测序生物信息分析流程图", type="img-width-normal", desc=" ", path="src/images/Analysis_process2.png" )
addElement(dom, root,"h1" ,name="生物信息学分析方法和结果" ,type="type1", desc="生物信息学分析方法和结果" )
addElement(dom, root,"h2", name="测序数据质控" ,type="type1" ,desc="测序数据质控" )
addElement(dom, root,"h3", name="测序数据介绍",type="type1", desc="测序数据介绍" )
addElement(dom, root,"p" ,type="type1" ,desc="高通量测序（%s等测序平台）得到的原始图像数据文件，经碱基识别（Base Calling）分析转化为原始测序序列（Sequenced Reads），我们称之为Raw Data或Raw Reads，结果以FASTQ（简称为fq）文件格式存储，其中包含测序序列（Reads）的序列信息以及其对应的测序质量信息。测序样品中真实数据随机截取结果如下：" %Sequencing_platform)
if Sequencing_platform == 'Illumina':
	addElement(dom, root,"pic" ,name="图 Fastq格式", type="type1" ,desc="注：FASTQ格式文件中每个Read由四行描述，其中第一行以“@”开头，随后为Illumina测序识别符（Sequence Identifiers）和描述文字（选择性部分）；第二行是碱基序列；第三行以“+”开头，随后为Illumina测序识别符（选择性部分）；第四行是对应序列的测序质量。" ,path="src/images/fastq_format.png" ) 
	addElement(dom, root,"p" ,type="type1" ,desc="Illumina测序识别符（Sequence Identifiers）详细信息如下：" )
	addElement(dom, root,"table", name="Illumina 测序标识详细信息表", type="full", desc="", path="src/images/table3" ) 
	addElement(dom, root,"p", type="type1" ,desc="通过使用第四行中每个字符对应的ASCII值进行计算，即得到对应第二行碱基的测序质量值。如果测序错误率用e表示，%s的碱基质量值用Qphred表示，则有下列关系："%Sequencing_platform)
	addElement(dom, root,"pic", name="图 公式1", type="img-width-normal" ,desc=" " ,path="src/images/formula1.png" )
else:
	addElement(dom, root,"pic" ,name="图 Fastq格式", type="type1" ,desc="注：FASTQ格式文件中每个Read由四行描述，其中第一行以“@”开头，随后为BGI测序识别符（Sequence Identifiers）和描述文字（选择性部分）；第二行是碱基序列；第三行以“+”开头，随后为BGI测序识别符（选择性部分）；第四行是对应序列的测序质量。" ,path="src/images/fastq_format_BGI.png" ) 
	addElement(dom, root,"p" ,type="type1" ,desc="BGI测序识别符（Sequence Identifiers）详细信息如下：" )
	addElement(dom, root,"table", name="BGI 测序标识详细信息表", type="full", desc="", path="src/images/table3_BGI" )
	addElement(dom, root,"p", type="type1" ,desc="通过使用第四行中每个字符对应的ASCII值进行计算，即得到对应第二行碱基的测序质量值。如果测序错误率用e表示，%s的碱基质量值用Qphred表示，则有下列关系："%Sequencing_platform)
	addElement(dom, root,"pic", name="图 公式1", type="img-width-normal" ,desc=" " ,path="src/images/formula1.png" )
addElement(dom, root,"p" ,type="type1" ,desc="Illumina Casava 1.8版本测序错误率与测序质量值简明对应关系如下表所示：" )
addElement(dom, root,"table" ,name="测序错误率与测序质量值对应关系", type="full" ,desc="注：碱基识别（Base Calling）分析软件：Illumina Casava 1.8版本；测序参数：双端测序(Paired end)；测序序列读长：%sbp（或者单位为循环数(cycle))"%Read_length, path="src/images/table4" ) 
addElement(dom, root,"h3", name="数据质量统计" ,type="type1", desc="数据质量统计" )
addElement(dom, root,"p", type="type1" ,desc="测序得到的原始测序序列（Sequenced Reads）或者Raw Reads，里面含有带接头的、低质量的Reads，为了保证信息分析质量，对Raw Reads进行过滤，得到Clean Reads，用于后续信息分析。数据过滤的主要步骤如下：(1) 去除带接头（adapter）的reads；(2) 过滤N含量超过10%的reads；(3) 去除质量值低于10的碱基超过50%的reads。" )
addElement(dom, root,"p", type="type1", desc="各样品测序产出数据评估结果见下表：" )
addElement(dom, root,"table", name="样品测序数据评估统计" ,type="full" ,path="Up_load/Dataassess/sample_data_assess.list",desc="注：BMK_ID：百迈客对项目样品的统一编号；Clean_Reads：过滤后的reads数； Clean_Base：过滤后的碱基数，Clean_Reads数乘以序列长度；Q20(%)：质量值大于等于20的碱基占总碱基数的百分比；Q30(%)：质量值大于等于30的碱基占总碱基数的百分比；GC(%)：样品GC含量，即G和C类型的碱基占总碱基的百分比。"  ) 
addElement(dom, root,"h3", name="碱基测序质量分布", type="type1", desc="碱基测序质量分布" )
addElement(dom, root,"p" ,type="type1", desc="每个碱基测序错误率是通过测序 Phred数值（Phred score，Qphred）通过公式一转化得到，而Phred数值是在碱基识别（Base Calling）过程通过一种预测碱基判别发生错误概率模型计算得到的，对应关系如下表所显示：" )
addElement(dom, root,"table", name="预测碱基判别发生错误概率", type="full", desc="", path="src/images/table31" ) 
addElement(dom, root,"p" ,type="type1", desc="在%s测序系统测序时，首先会对文库进行芯片制备，目的是将文库DNA模板固定到芯片上，在固定DNA模板的过程中，每个DNA分子会形成一个簇，一个簇就是一个测序位点，在进行固定过程中极少量的簇与簇之间物理位置会发生重叠，在测序时，测序软件通过前4个碱基对这些重叠的点进行分析和识别，将这些重叠点位置分开，保证每个点测到的是一个DNA分子，因此测序序列5′端前几个碱基的错误率相对较高。另外测序错误率会随着测序序列（Sequenced Reads）的长度的增加而升高，这是由于测序过程中化学试剂的消耗而导致的。因此在进行碱基测序质量分布分析时，样品的碱基质量分布在前4个碱基和后十几个碱基的质量值会低于中间测序碱基，但其质量值都高于Q30%%，根据质量值和错误率的关系，我们将质量值转换成错误率，绘制错误率分布图如下："%(Sequencing_platform))
pic_list=addElement(dom, root,"pic_list" ,name="样品碱基错误率分布" ,type="type1", desc="注：横坐标为reads的碱基位置，纵坐标为单碱基错误率。前%sbp为双端测序序列的第一端测序Reads的错误率分布情况，后%sbp为另一端测序reads的错误率分布情况。"%(Read_length,Read_length))
qualityImages=glob.glob(analysis_dir+'/Up_load/Dataassess/*.quality.png')
addElementListImage(dom,pic_list,'pic',qualityImages,web_dir)
addElement(dom, root,"h3", name="碱基类型分布" ,type="type1", desc="碱基类型分布" )
addElement(dom, root,"p", type="type1", desc="碱基类型分布检查用于检测有无AT、GC分离现象，而这种现象可能是测序或者建库所带来的，并且会影响后续分析。高通量所测序列为基因组随机打断后的DNA片段，由于位点在基因组上的分布是近似均匀的，同时，G/C、A/T含量也是近似均匀的。因此，根据大数定理，在每个测序循环上，GC、AT含量应当分别相等，且等于基因组的GC、AT含量。同样因为重叠簇的关系会导致样品前几个碱基AT、GC不等波动较大，高于其他测序区段，而其它区段的GC、AT的含量相等，且分布均匀无分离现象,下图所示：" )
pic_list=addElement(dom, root,"pic_list" ,name="样品各碱基比例分布" ,type="type1", desc="注：横坐标为reads的碱基位置，纵坐标为碱基所占的比例；不同颜色代表不同的碱基类型，绿色代表碱基G，蓝色代表碱基C，红色代表碱基A，紫色代表碱基T，灰色代表测序中识别不出的碱基N。前%sbp为双端测序序列的第一端测序Reads的碱基分布，后%sbp为另一端测序reads的碱基分布。每个cycle代表测序的每个碱基，如第一cycle即表示该项目所有测序reads在第一个碱基的A、T、G、C、N的分布情况。该图的结果显示AT、CG碱基基本未发生分离，曲线较平缓，测序结果正常。"%(Read_length,Read_length))
acgtnImages=glob.glob(analysis_dir+'/Up_load/Dataassess/*.acgtn.png')
addElementListImage(dom,pic_list,'pic',acgtnImages,web_dir)

addElement(dom, root,"h2", name="与参考基因组比对统计", type="type1", desc="与参考基因组比对统计" )
addElement(dom, root,"p", type="type1", desc="重测序获得的测序reads需要重新定位到参考基因组上，才可以进行后续变异分析。bwa<a href=\"#ref1\">[1]</a>软件主要用于二代高通量测序（如%s等测序平台）得到的短序列与参考基因组的比对。通过比对定位Clean Reads在参考基因组上的位置，统计各样品的测序深度、基因组覆盖度等信息，并进行变异的检测。"%(Sequencing_platform))

addElement(dom, root,"h3",name="比对结果统计",type="type1",desc="比对结果统计" )
addElement(dom, root,"p", type="type1", desc="与参考基因组比对，统计样品比对率等信息。比对率：即可以定位到参考基因组上的Clean Reads占总的Clean Reads数比例，如果参考基因组选择合适，且相关实验过程不存在污染，测序reads的比对率会高于70%。另外，比对率的高低受测序物种与参考基因组亲缘关系远近、参考基因组组装质量高低及reads测序质量有关，物种越近缘、参考基因组组装越完整、测序reads质量越高，则可以定位到参考基因组的reads也越多，比对率越高。" )
addElement(dom, root,"p", type="type1", desc="比对统计方法如下：" )
addElement(dom, root,"p",type="type1", desc="Clean Reads：统计Clean Reads数据文件，每四行为一个单位，双端分别统计，即read1和read2记为2条reads。")
addElement(dom, root,"p",type="type1", desc="Mapped(%)：定位到参考基因组的Clean Reads数占所有Clean Reads数的百分比，使用samtools flagstat命令实现。")
addElement(dom, root,"p",type="type1", desc="Properly mapped(%)：双端测序序列均定位到参考基因组上且距离符合测序片段的长度分布，使用samtools flagstat命令实现。")

addElement(dom, root,"p", type="type1", desc="样品的比对结果如下表所示:" )
addElement(dom, root,"table" ,name="比对结果统计", type="full",  path="Up_load/Mapping/%s.map_stat.xls"%Project,desc="注：BMK ID：百迈客对项目样品的统一编号；Total_reads：Clean Reads数；Mapped(%)：定位到参考基因组的Clean Reads数占所有Clean Reads数的百分比；Properly_mapped(%)：双端测序序列均定位到参考基因组上且距离符合测序片段的长度分布。" ) 
addElement(dom, root,"h3", name="插入片段分布统计", type="type1", desc="插入片段分布统计" )
addElement(dom, root,"p", type="type1", desc="通过检测双端序列在参考基因组上的起止位置，可以得到样品DNA打断后得到的测序片段的实际大小，即插入片段大小(Insert Size)，它是信息分析时的一个重要参数。插入片段大小的分布一般符合正态分布，且只有一个单峰，Insert Size分布图可以展示各个样品的插入片段的长度分布情况。每个样品测序数据插入片段大小的分析使用picard软件工具包中CollectInsertSizeMetric.jar软件实现。" )
pic_list=addElement(dom, root,"pic_list", name="插入片段分布图", type="type1", desc="注：横坐标为插入片段长度，纵坐标为其对应的reads数。")
insertImages=glob.glob(analysis_dir+'/Up_load/Mapping/*insert.cut.png')
addElementListImage(dom,pic_list,'pic',insertImages,web_dir)
addElement(dom, root,"p", type="type1", desc="由上图可知，插入片段长度分布符合正态分布，说明测序数据文库构建无异常。" )
addElement(dom, root,"h3", name="深度分布统计", type="type1", desc="深度分布统计" )
addElement(dom, root,"p", type="type1", desc="Reads定位到参考基因组后，可以统计参考基因组上碱基的覆盖情况。参考基因组上被reads覆盖到的碱基数占基因组的百分比称为基因组覆盖度；碱基上覆盖的reads数为覆盖深度。基因组覆盖度可以反映参考基因组上变异检测的完整性，覆盖到的区域越多，可以检测到的变异位点也越多。覆盖度主要受测序深度以及样品与参考基因组亲缘关系远近的影响。基因组的覆盖深度会影响变异检测的准确性，在覆盖深度较高的区域（非重复序列区），变异检测的准确性也越高。另外，若基因组上碱基的覆盖深度分布较均匀，也说明测序随机性较好。样品的碱基覆盖深度分布曲线和覆盖度分布曲线见下图：" )
pic_list=addElement(dom, root,"pic_list", name="样品深度分布", type="type1", desc="注：上图反映了测序深度分布的基本情况，横坐标为测序深度，左纵坐标为该深度对应的碱基所占百分比，对应红色曲线，右纵坐标为该深度及以下的碱基所占百分比，对应蓝色曲线。")
depthImages=glob.glob(analysis_dir+'/Up_load/Mapping/*depth.dis.cut.png')
addElementListImage(dom,pic_list,'pic',depthImages,web_dir)
addElement(dom, root,"p" ,type="type1", desc="各样品的平均覆盖深度和各深度对应的基因组覆盖比例如下表所示：" )
addElement(dom, root,"table" ,name="样品覆盖深度和覆盖度比例统计", type="full", desc="注：BMK ID：百迈客对项目样品的统一编号；Ave_depth：样品平均覆盖深度；后面三列代表覆盖深度在给定深度及以上的碱基数占参考基因组总碱基数的比例，分别是1X、5X、10X以上的比例。" ,path="Up_load/Mapping/%s.depth_stat.xls"%Project ) 
addElement(dom, root,"p", type="type1", desc="根据染色体各位点的覆盖深度情况进行作图，若覆盖深度在染色体上的分布比较均匀，则可以认为测序随机性比较好。样品的染色体覆盖深度分布图如下所示：" )
pic_list=addElement(dom, root,"pic_list", name="样品染色体覆盖深度分布图", type="type1", desc="注：横坐标为染色体位置，纵坐标为染色体上对应位置的覆盖深度取对数（log2）得到的值。")
covImages=glob.glob(analysis_dir+'/Up_load/Mapping/*depth.cov.txt.png')
addElementListImage(dom,pic_list,'pic',covImages,web_dir)
addElement(dom, root,"p", type="type1", desc="由上图可以看出基因组被覆盖的较均匀，说明测序随机性较好。图上深度不均一的地方可能是由于重复序列、PCR偏好性引起的。" )

addElement(dom, root,"h2", name="变异检测及注释", type="type1", desc="变异检测及注释" )
addElement(dom, root,"h3", name="变异检测工具及方法", type="type1", desc="变异检测工具及方法" )
addElement(dom, root,"p", type="type1", desc="SNP（Single Nucleotide Polymorphism，单核苷酸多态性）和small InDel（small Insertion and Deletion，小片段的插入与缺失）的检测主要使用GATK<a href=\"#ref2\">[2]</a>软件工具包实现。根据Clean Reads在参考基因组的定位结果，使用Picard<a href=\"#ref3\">[3]</a>过滤冗余reads（MarkDuplicates），以保证检测结果的准确性。然后使用GATK的HaplotypeCaller（局部单体型组装）算法进行SNP和InDel的变异检测，每个样本先各自生成gVCF，再进行群体joint-genotype。最后过滤，并得到最终的变异位点集。" )
addElement(dom, root,"h3", name="变异检测结果质量控制", type="type1", desc="变异检测结果质量控制" )
addElement(dom, root,"p", type="type1", desc="变异结果经过严格的过滤，保证变异结果的可靠性，主要过滤参数如下：" )
addElement(dom, root,"p", type="type1", desc="（1）基于bcftools中的子程序vcfutils.pl（varFilter -w 5 -W 10）对INDEL附近5bp内的SNP以及相邻INDEL在10bp内的过滤掉；" )
addElement(dom, root,"p", type="type1", desc="（2）clusterSize 2  clusterWindowSize 5，表示5bp窗口内的变异数量不应该超过2个；" )
addElement(dom, root,"p", type="type1", desc="（3）QUAL < 30，Phred格式的质量值，表示该位点存在variant变异的可能性。质量值低于30的则过滤掉；" )
addElement(dom, root,"p", type="type1", desc="（4）QD < 2.0，变异质量值（Quality）除以覆盖深度（Depth）得到的比值，覆盖深度是这个位点上所有含有变异碱基的样本的覆盖深度之和。QD低于2.0的则过滤掉；" )
addElement(dom, root,"p", type="type1", desc="（5）MQ < 40，所有比对至该位点上的read的比对质量值的均方根。MQ低于40的则过滤掉；" )
addElement(dom, root,"p", type="type1", desc="（6）FS > 60.0，通过Fisher检验的p-value转换而来的值，描述的是测序或者比对时对于只含有变异的read以及只含有参考序列碱基的read是否存在着明显的正负链特异性。也就是说，不会出现链特异的比对结果，FS应该接近于零。FS高于60的则过滤掉；" )
addElement(dom, root,"p", type="type1", desc="（7）其它变异过滤参数采用GATK官方指定的默认值处理。" )
addElement(dom, root,"h3", name="高质量变异结果展示", type="type1", desc="高质量变异结果展示" )
addElement(dom, root,"p", type="type1", desc="变异结果使用vcf文件格式展示。vcf文件包括注释行、标题行和数据行三部分。其中注释行包含文件数据行的INFO和FORMAT列中使用的各种标识符的意义解释，而标题行和数据行包含各样品的变异检测结果信息，格式如下所示：" )
addElement(dom, root,"table", name="SNP变异结果信息示范列表", type="full", desc="", path="src/images/table9" )
addElement(dom, root,"p", type="type1", desc="本项目SNP变异结果VCF文件详见：Up_load/SNP/*.vcf。" )
addElement(dom, root,"p", type="type1", desc="各列意义说明如下：" )
addElement(dom, root,"table", name="各列意义说明", type="full", desc="" ,path="src/images/table10" ) 
addElement(dom, root,"p", type="type1", desc="其中FORMAT（GT:AD:DP:GQ:PL）：关键字之间用冒号隔开。" )
addElement(dom, root,"p", type="type1", desc="GT：genotype，表示这个样本的基因型，对于一个二倍体生物来说，GT值表示的是这个样本在这个位点所携带的两个等位基因的类型，0表示跟REF一样；1表示跟ALT一样，0/0表示纯合且跟REF一致；0/1表示杂合，两个allele一个是ALT一个是REF；1/1表示纯和且都为ALT；" )
addElement(dom, root,"p", type="type1", desc="AD：allele depth，对应两个以逗号隔开的值，这两个值分别表示覆盖到REF和ALT碱基的reads数，相当于支持REF和支持ALT的测序深度；" )
addElement(dom, root,"p", type="type1", desc="DP：depth of coverage，覆盖到这个位点的总的reads数量，相当于这个位点的深度；" )
addElement(dom, root,"p", type="type1", desc="GQ：Quality of the assigned genotype，表示最可能的基因型的质量值；" )
addElement(dom, root,"p", type="type1", desc="PL：Normalized Phred-scaled likelihoods of the possible genotypes，对应3个以逗号隔开的值，这三个值分别表示该位点基因型是0/0，0/1，1/1的没经过先验的标准化Phred-scaled似然值（L）。如果转换成支持该基因型概率（P）的话，由于L=-10lgP，那么P=10^（-L/10），当L值为0时，P=10^0=1。因此，这个值越小，支持概率就越大，也就是说是这个基因型的可能性越大。" )
addElement(dom, root,"p", type="type1", desc="vcf文件的详细说明信息见网页：<a href=\"http://gatkforums.broadinstitute.org/discussion/1268/how-should-i-interpret-vcf-files-produced-by-the-gatk\" target=\"_blank、\">http://gatkforums.broadinstitute.org/discussion/1268/how-should-i-interpret-vcf-files-produced-by-the-gatk</a>" )
addElement(dom, root,"p", type="type1" ,desc="为了确保SNP的可信性，对检测的SNP的reads支持数，相邻SNP的距离统计累积分布。" )
addElement(dom, root,"pic", name="SNP质量分布图", type="type1" ,desc="注：左边为SNP reads 支持数目累积图，右边为相邻SNP之间的距离累积图。 " ,path="Up_load/SNP/SNP.quality.distribution.png" )
addElement(dom, root,"p", type="type1", desc="SNP类型的变异分为转换和颠换两种，同种类型碱基之间突变称为转换（Transition），如嘌呤与嘌呤之间、嘧啶与嘧啶之间的变异，不同类型碱基之间的突变称为颠换（Transversion），如嘌呤与嘧啶之间的变异。一般来说转换比颠换更容易发生，故转换/颠换（Ti/Tv）的比例一般大于1，具体数值和所测物种有关。对于二倍体或者多倍体物种，若同源染色体上的某一SNP位点均为同一种碱基，则该SNP位点称为纯合SNP位点；若同源染色体上的SNP位点包含不同类型的碱基，则该SNP位点称为杂合SNP位点。纯合SNP数量越多，则样品与参考基因组之间差异越大，杂合SNP数量越多，则样品的杂合程度越高，具体结果和样品的材料选择有关。样品与参考基因组之间的SNP检测结果如下所示。" )
addElement(dom, root,"table", name="检测得到的SNP统计", type="full", desc="", path="Up_load/SNP/%s.raw.filter.snp.snp.stat"%Project ) 
addElement(dom, root,"p", type="type1", desc="各列说明如下表所示：" )
addElement(dom, root,"table", name="各列意义说明", type="full", desc="", path="src/images/table12" ) 
if sampleNum ==1:
    addElement(dom, root,"h3", name="样品SNP的检测(单个样品)", type="type1", desc="样品SNP的检测(单个样品)" )
    addElement(dom, root,"p", type="type1" ,desc="根据样品与参考基因组的比对结果，汇总样品与参考基因组的变异位点，样品的SNP列表文件格式如下所示：" )
    addElement(dom, root,"table", name="样品SNP列表示意" ,type="type1", desc="", path="Up_load/SNP/%s.raw.filter.snp.snp"%Project )
    addElement(dom, root,"p", type="type1" ,desc="详细列表数据请见：Up_load/SNP/%s.raw.filter.snp.snp。"%Project )
else:
    addElement(dom, root,"h3", name="样品之间SNP的检测(多个样品)", type="type1", desc="样品之间SNP的检测(多个样品)" )
    addElement(dom, root,"p", type="type1" ,desc="根据样品与参考基因组的比对结果，汇总样品之间所有有差异的变异位点，样品间的SNP列表文件格式如下所示：" )
    addElement(dom, root,"table", name="样品间SNP列表示意" ,type="type1", desc="", path="Up_load/SNP/%s.raw.filter.snp.DEG.snp"%Project )
    addElement(dom, root,"p", type="type1" ,desc="详细列表数据请见：Up_load/SNP/%s.raw.filter.snp.DEG.snp。" %Project)
addElement(dom, root,"p", type="type1" ,desc="SNP基因型的编码采用标准核苷酸符号，符号表如下所示：" )
addElement(dom, root,"table", name="SNP基因型的编码", type="full", desc="", path="src/images/table14" )
addElement(dom, root,"p", type="type1" ,desc="全基因组SNP突变可以分成6类。以T:A>C:G为例，此种类型SNP突变包括T>C和A>G。由于测序数据既可比对到参考基因组的正链，也可比对到参考基因组的负链，当T>C类型突变出现在参考基因组正链上，A>G类型突变即在参考基因组负链的相同位置，所以将T>C和A>G划分成一类。" )
addElement(dom, root,"pic", name="SNP", type="type1" ,desc="注：纵坐标表示SNP的数目，横坐标表示SNP突变类型。 " ,path="Up_load/SNP/SNP.mutation.distribution.png" )
if sampleNum == 1:
	pass
elif os.path.exists(analysis_dir+"/Up_load/SNP/SNP.venn.png"):
    addElement(dom, root,"p", type="type1", desc="样品间SNP的统计结果如下图所示：" )
    addElement(dom, root,"pic", name="样品间SNP统计Venn图", type="type1" ,desc="注：变异位点数量venn统计只考虑位置是否相同，不考虑基因型是否相同。" ,path="Up_load/SNP/SNP.venn.png" )
else:
    addElement(dom, root,"p", type="type1", desc="样品间SNP的统计结果如下表所示：" )
    addElement(dom, root,"table", name="样品间SNP统计表", type="full" ,desc="BMK ID：百迈客对项目样品的统一编号；表中各数值为对应的横纵两样品之间的SNP数。" ,path="Up_load/SNP/DEG_snp.stat" )

addElement(dom, root,"h3", name="SNP结果注释", type="type1" ,desc="SNP结果注释" )
addElement(dom, root,"p", type="type1", desc="SnpEff<a href=\"#ref4\">[4]</a>是一款用于注释变异（SNP、Small InDel）和预测变异影响的软件。根据变异位点在参考基因组上的位置以及参考基因组上的基因位置信息，可以得到变异位点在基因组发生的区域（基因间区、基因区或CDS区等），以及变异产生的影响（同义、非同义突变等）。软件可以使用vcf格式文件作为输入和输出。输出结果会在vcf文件的INFO列添加以下字段：EFF= Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_Length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype_Number [ | ERRORS | WARNINGS ] )" )

addElement(dom, root,"p" ,type="type1", desc="各标识符说明如下：" )
addElement(dom, root,"table",name="标识符说明", type="full",desc="注：以上的结果若无法得到，则其对应列为空。具体说明可参见SnpEff的说明文档：<a href=\"http://snpeff.sourceforge.net/SnpEff_manual.html#output\" target=\"_blank\">http://snpeff.sourceforge.net/SnpEff_manual.html#output。</a> ", path="src/images/table16" ) 
pic_list=addElement(dom, root,"pic_list", name="SNP注释结果统计", type="type1", desc="")
snpAnnImages=glob.glob(analysis_dir+'/Up_load/SNP/*pie.png')
addElementListImage(dom,pic_list,'pic',snpAnnImages,web_dir)
addElement(dom, root,"p", type="type1", desc="上图各标示意义说明如下表所示：" )
addElement(dom, root,"table",name="功能区意义说明",type="full" ,desc="",path="src/images/table18" )
addElement(dom, root,"h3", name="样品与参考基因组间Small InDel的检测", type="type1", desc="样品与参考基因组间Small InDel的检测" )
addElement(dom, root,"p",type="type1", desc="根据样品的Clean Reads在参考基因组上的定位结果，检测样品与参考基因组之间是否存在小片段的插入与缺失（Small InDel）。样品的插入缺失使用GATK检测。Small InDel变异一般比SNP变异少，同样反映了样品与参考基因组之间的差异，并且编码区的InDel会引起移码突变，导致基因功能上的变化。见下表" )

addElement(dom, root,"table", name="全基因组和编码区InDel统计",type="full", path="Up_load/INDEL/%s.combine.Indel.stat"%Project,desc="注：CDS：编码区的InDel统计；Genome：全基因组范围的InDel统计；Insertion：检测到的插入数量；Deletion：检测到的删除数量；Het：杂合InDel的数量；Homo：纯合InDel的数量；Total：检测到的InDel总数（去除重复）。"  ) 

addElement(dom, root,"p", type="type1", desc="根据样品在CDS区和全基因组范围的InDel长度进行统计，其长度分布图见下图：" )
addElement(dom, root,"pic", name="全基因组和编码区InDel长度分布图" ,type="type1", desc="注：纵坐标为InDel长度(10bp以内)，大于0为Insertion，小于0为Deletion，横坐标为对应的数量。", path="Up_load/INDEL/all.sample.indel.length.distribution.png" ) 

addElement(dom, root,"h3", name="样品之间Small InDel的检测", type="type1", desc="样品之间Small InDel的检测" )
addElement(dom, root,'p', type="type1" ,desc="根据样品与参考基因组的Small InDel检测结果，样品测序数据部分比较结果见下表。" )
addElement(dom, root,"table", name="样品Small InDel测序数据统计", type="type1", desc="", path="Up_load/INDEL/%s.raw.filter.indel.anno.gatk.list"%Project )
addElement(dom, root,"p", type="type1", desc="样品Small indel测序数据详见：Up_load/INDEL/%s.raw.filter.indel.anno.gatk.list。"%Project)
addElement(dom, root,"p", type="type1", desc="各列意义说明如下：" )
addElement(dom, root,"table", name="各列意义说明", type="full", desc="", path="src/images/table21" ) 

if sampleNum == 1:
	pass
elif os.path.exists(analysis_dir+"/Up_load/INDEL/INDEL.venn.png"):
    addElement(dom, root,"p", type="type1", desc="样品间INDEL的统计结果如下图所示：" )
    addElement(dom, root,"pic", name="样品间INDEL统计Venn图", type="type1" ,desc="注：变异位点数量venn统计只考虑位置是否相同（indel的起始位置），不考虑基因型是否相同。" ,path="Up_load/INDEL/INDEL.venn.png" )
else:
    addElement(dom, root,"p", type="type1", desc="两两样品间的Small InDel检测的统计结果如下表所示：" )
    addElement(dom, root,"table", name="样品间InDel统计结果", type="full" ,desc="注：BMK ID：百迈客对项目样品的统一编号；表中各数值为对应的横纵两样品之间的Small InDel数。" ,path="Up_load/INDEL/DEG_indel.stat" )

addElement(dom, root,"h3", name="Small InDel的注释", type="type1", desc="Small InDel的注释" )
addElement(dom, root,"p", type="type1", desc="根据样品检测得到的Small InDel位点在参考基因组上的位置信息，对比参考基因组的基因、CDS位置等信息(一般在gff文件中)，可以注释InDel位点是否发生在基因间区、基因区或CDS区、是否为移码突变等。Small InDel的注释通过SnpEff软件实现。发生移码突变的InDel可能会导致基因功能的改变，具体注释结果见下图：" )

pic_list=addElement(dom, root,"pic_list", name="InDel注释结果统计", type="type1", desc="")
snpAnnImages=glob.glob(analysis_dir+'/Up_load/INDEL/*pie.png')
addElementListImage(dom,pic_list,'pic',snpAnnImages,web_dir)
addElement(dom, root,"p", type="type1", desc="上图各标示说明如下表所示：" )
addElement(dom, root,"table", name="各功能区意义说明", type="full", desc="", path="src/images/table24" ) 

reference = 5
#SV检测与注释
if(os.path.exists(analysis_dir+'/Up_load/SV')):
    addElement(dom, root,"h2", name="SV检测与注释", type="type1", desc="SV检测与注释" )
    addElement(dom, root,"h3", name="SV的检测", type="type1", desc="SV的检测" )
    addElement(dom, root,"p", type="type1", desc="基因组结构变异（SV）是指物种个体在染色体结构上的变异，比如大片段插入、缺失、倒位、易位等。结构变异一般使用breakDancer检测，breakDancer首先基于序列与参考基因组的比对结果，得到测序数据文库的插入片段大小和方差，然后通过查找序列和参考基因组之间的异常比对结果（插入片段发生偏离、比对方向不一致等），寻找可能的结构变异。SV结果文件包含标题行和数据行，如下所示：" )
    addElement(dom, root,"table", name="SV结果文件展示(样品%s)"%Sample_Project,type="type1", desc=" ",path="Up_load/SV/%s.%s.%s.dedup.realn.bam.max"%(Project,Project,Sample_Project))
    addElement(dom, root,"p", type="type1", desc="（样品R*）SV结果文件详见：Up_load/SV/%s.%s.R*.dedup.realn.bam.max"%(Project,Project) )
    addElement(dom, root,"p", type="type1", desc="标题行每列意义说明如下表所示：" )
    addElement(dom, root,"table", name="各列意义说明", type="full" ,desc="注：关于SV的结果文件的相关说明可以参见：<a href=\"https://github.com/kenchen/breakdancer#readme\" target=\"_blank\">https://github.com/kenchen/breakdancer#readme</a>" ,path="src/images/table26" ) 
    addElement(dom, root,"p", type="type1", desc="利用BreakDancer<a href=\"#ref%d\">[%d]</a>软件，基于Pair-end reads比对到参考基因组上面的关系及实际Insert Size大小检测样品与参考基因组间的插入（Insertion, INS）、缺失（Deletion, DEL）、倒位（Inversion, INV）、染色体内部易位（Intra-chromosomal Translocation, ITX）、染色体间易位（Inter-chromosomal Translocation, CTX)。%s个样品检测得到的各类型SV数量统计见表：" %(reference,reference,sampleNum))
    addElement(dom, root,"table", name="SV数量统计表", type="full", desc="注：BMK ID：百迈客对项目样品的统一编号；SV：结构变异总数目；INS：插入类型变异数量；DEL：缺失类型变异数量；INV：反转类型变异数量；ITX：染色体内易位类型变异数量；CTX：染色体间易位类型变异数量；UN：复杂的结构变异。", path="Up_load/SV/%s.sv.stat.xls"%Project ) 
    addElement(dom, root,"h3", name="SV的注释", type="type1", desc="SV的注释" )
    addElement(dom, root,"p", type="type1", desc="根据样品检测得到的SV变异在参考基因组上的位置信息，对比参考基因组的基因、CDS位置等信息（一般在gff文件中），可以注释SV变异是否发生在基因间区、基因区或CDS区等。对缺失（DEL）、插入（INS）、反转（INV）3种类型的结构变异注释进行注释，统计结果如表所示：" )
    addElement(dom, root,"table", name="SV注释结果统计表", type="full", desc="注：BMK ID：百迈客对项目样品的统一编号；Type：SV类型；Exon：处在外显子区的变异；Intron：处在内含子区的变异；Intergenic：处在基因间区的变异。", path="Up_load/SV/%s.sv.anno.stat"%Project ) 
    reference = reference + 1
    
#CNV检测与注释
if(os.path.exists(analysis_dir+'/Up_load/CNV')):
    addElement(dom, root,"h2", name="CNV检测", type="type1", desc="CNV检测" )
    addElement(dom, root,"p", type="type1",desc="利用FREEC<a href=\"#ref%d\">[%d]</a>，通过样品测序reads在参考基因组上的深度分布来检测CNV，并绘制参考基因组上的Copy Number Gain 及Loss的分布。样品%s的FREEC的部分结果文件：" %(reference,reference,Sample_Project))
    addElement(dom, root,"table", name="FREEC结果文件", type="type1", desc="注：chr：表示染色体；strat：表示CNV起始位置；end：表示CNV的结束位置；predicted copy number：表示预测到的拷贝数；type of alteration：表示CNV类型。",path="Up_load/CNV/%s.%s.dedup.realn.bam_CNVs"%(Project,Sample_Project))
    addElement(dom, root,"p", type="type1", desc="（样品R*）FREEC结果文件详见：Up_load/CNV/%s.R*.dedup.realn.bam_CNVs。" %(Project) )
    reference = reference + 1
#circos图
if(os.path.exists(analysis_dir+'/Up_load/Circos')):
    addElement(dom, root,"h2", name="各变异在基因组上的分布", type="type1", desc="各变异在基因组上的分布" )
    addElement(dom, root,"p", type="type1", desc="检测得到的各类型变异结果分布使用circos图展示，使用circos软件作图，软件网址：<a href=\"http://circos.ca/\" target=\"_blank\">http://circos.ca/</a>。各样品的各类型变异在染色体上的分布见下图：" )
    pic_list=addElement(dom, root,"pic_list", name="样品各类型变异在染色体的分布", type="type1", desc="注：从外到里依次为：染色体坐标、SNP密度分布、InDel密度分布、CNV、SV（INS、DEL、INV、ITX(red line)、CTX(green line)）在基因组上的分布(单位为M）。")
    circosImages=glob.glob(analysis_dir+'/Up_load/Circos/*.circos.png')
    addElementListImage(dom,pic_list,'pic',circosImages,web_dir)

#差异基因分析
if(os.path.exists(analysis_dir+'/Up_load/Diff_analysis')):
    addElement(dom, root,"h2", name="DNA水平变异基因分析", type="type1", desc="DNA水平变异基因分析" )
    addElement(dom, root,"h3", name="DNA水平的变异基因挖掘", type="type1", desc="DNA水平的变异基因挖掘" )
    addElement(dom, root,"p", type="type1", desc="发生在CDS区的变异可能会引起基因功能的改变，通过寻找参考基因组与样品间发生非同义突变SNP、CDS区发生的InDel与SV的基因，寻找样品与参考基因组之间可能存在功能差异的基因。各样品与参考基因组之间的变异如下表所示：" )
    addElement(dom, root,"table", name=" 各种变异产生的差异基因的分类统计", type="full", desc="注：BMK_ID：百迈客对项目样品的统一编号；Genes with Non-synonymous：发生非同义突变基因数量，同一基因发生多个非同义突变的不做重复统计；Genes with InDel：发生小的插入缺失(Small InDel)的基因数量；Genes with SV：发生SV变异基因数量。" ,path="Up_load/Diff_analysis/Diff_gene.stat" ) 
    addElement(dom, root,"h3", name="DNA水平的变异基因功能注释", type="type1", desc="DNA水平的变异基因功能注释" )
    ref1 = reference
    ref2 = ref1 + 1
    ref3 = ref2 + 1
    ref4 = ref3 + 1
    ref5 = ref4 + 1
    ref6 = ref5 + 1
    addElement(dom, root,"p", type="type1" ,desc="通过BLAST<a href=\"#ref%d\">[%d]</a>将变异基因与NR<a href=\"#ref%d\">[%d]</a>，SwissProt<a href=\"#ref%d\">[%d]</a>，GO<a href=\"#ref%d\">[%d]</a>，COG<a href=\"#ref%d\">[%d]</a>，KEGG<a href=\"#ref%d\">[%d]</a>等功能数据库比对，得到这些基因的注释，用以分析基因功能。" %(ref1,ref1,ref2,ref2,ref3,ref3,ref4,ref4,ref5,ref5,ref6,ref6))
    addElement(dom, root,"p" ,type="type1" ,desc="NR数据库：是NCBI数据库的非冗余蛋白质数据库，包含了SwissProt、PIR(Protein Information Resource)、PRF(Protein Research Foundation)、PDB(Protein Data Bank)蛋白质数据库及从GenBank和RefSeq的CDS数据翻译过来的蛋白质数据。" )
    addElement(dom, root,"p" ,type="type1" ,desc="SwissProt数据库：是经过注释的蛋白质序列数据库，由欧洲生物信息学研究所(EBI)维护。数据库由蛋白质序列条目构成，每个条目包含蛋白质序列、引用文献信息、分类学信息、注释等，注释中包括蛋白质的功能、转录后修饰、特殊位点和区域、二级结构、四级结构、与其它序列的相似性、序列残缺与疾病的关系、序列变异体和冲突等信息。SwissProt中尽可能减少了冗余序列，并与其它30多个数据建立了交叉引用，其中包括核酸序列库、蛋白质序列库和蛋白质结构库等。" )
    addElement(dom, root,"p" ,type="type1" ,desc="GO数据库：Gene Ontology（简称GO）是一个国际标准化的基因功能分类体系，提供了一套动态更新的标准词汇表（controlled vocabulary）来全面描述生物体中基因和基因产物的属性。Gene Ontology中最基本的概念是term。GO里面的每一个entry都有一个唯一的数字标记，形如GO:nnnnnnn，还有一个term名，比如“cell”、“fibroblast growth factor receptor binding”或者“signal transduction”。每个term都属于一个ontology，总共有三个ontology，它们分别是Molecular function（分子功能），Cellular component（细胞组件）和Biological process（生物过程）。" )
    addElement(dom, root,"p" ,type="type1" ,desc="COG数据库：是对基因产物进行直系同源分类的数据库。每个COG的蛋白都被假定来自祖先蛋白，分为orthologs和paralogs。Orthologs是指来自于不同物种的由垂直家系进化而来的蛋白，并且典型的保留与原始蛋白相同的功能。Paralogs是指那些在一定物种中来源于基因复制的蛋白，可能会进化出新的与原来有关的功能。" )
    addElement(dom, root,"p" ,type="type1" ,desc="KEGG数据库：KEGG整合了当前生物化学中关于化合物、反应和分子作用网络的知识，是系统分析基因产物在细胞中的代谢途径以及这些基因产物功能的数据库，有助于把基因及表达信息作为一个整体的网络进行研究。" )
    addElement(dom, root,"p" ,type="type1" ,desc="变异基因的注释列表文件如下表所示：" )
    addElement(dom, root,"table", name="样品%s变异基因的注释列表"%Sample_Project, type="type1", desc=" ", path="Up_load/Diff_analysis/%s/Integrated_Function.annotation.xls"%Sample_Project)
    addElement(dom, root,"p", type="type1", desc="(样品R*)变异基因注释列表详见:Up_load/Diff_analysis/R*/Integrated_Function.annotation.xls" )
    addElement(dom, root,"p", type="type1", desc="各样品的具体统计结果见下表:" )
    addElement(dom, root,"table", name="变异基因注释统计表", type="full", desc="注：Sample：样品名；Gene_Number：注释到数据库的差异基因数目。", path="Up_load/Diff_analysis/Allsample_Diffgene.list" ) 
    addElement(dom, root,"p", type="type1", desc="变异基因GO分类统计结果见下图：" )
    addElement(dom, root,"pic", name="样品%s变异基因GO注释聚类"%Sample_Project, type="type1",path="Up_load/Diff_analysis/%s/go_enrichment/%s.GO.png" %(Sample_Project,Project),desc="注：横坐标为GO各分类内容，纵坐标左边为基因数目所占百分比，右边为基因数目。") 
    addElement(dom, root,"p", type="type1" ,desc="变异基因COG分类统计结果见下图：" )
    addElement(dom, root,"pic" ,name="样品%s变异基因COG注释分类图"%Sample_Project,type="type1" ,path="Up_load/Diff_analysis/%s/Cog_Anno/%s.Cog.classfy.png"%(Sample_Project,Project) ,desc="注：横坐标为COG各分类内容，纵坐标为基因数目。在不同的功能类中，基因所占多少反映对应时期和环境下代谢或者生理偏向等内容，可以结合研究对象在各个功能类的分布作出科学的解释。" ) 
    snpAnnImages=glob.glob(analysis_dir+'/Up_load/Diff_analysis/%s/pathway/kegg_map/*.png'%Sample_Project)
    if snpAnnImages != []:
	    addElement(dom, root,"p", type="type1", desc="变异基因的代谢通路结果示意见下图：" )
	    pic_list=addElement(dom, root,"pic_list", name="样品%s变异基因通路代谢图"%Sample_Project, type="type1", desc="注：框内的数字代表enzyme的号码，说明对应基因与此酶相关，而整个通路是有很多种不同的酶经过复杂的生化反应形成的，变异基因中与此通路相关的均用红色的框标出，研究人员可以根据自己的研究对象，重点研究相关代谢通路的基因情况，通过代谢解释对应基因的根源。")
	    addElementListImage(dom,pic_list,'pic',snpAnnImages,web_dir)
	    
addElement(dom, root,"h1", name="结果文件查看说明", type="type1", desc="结果文件查看说明" )
addElement(dom, root,"p", type="type1", desc="(1) 上传目录中有Readme.txt说明，详细介绍了每个文件所代表的内容。上传的结果数据文件多以文本格式为主(fa文件、txt文件，detail文件，xls文件等)。在Windows系统下查看文件，推荐使用Editplus或 UltraEdit作为文本浏览程序，否则会因文件过大造成死机,或者使用PilotEdit Lite打开超大文本文件（<a href=\"http://www.pilotedit.com/\" target=\"_blank\">http://www.pilotedit.com/</a>）。在Unix或Linux系统下可以浏览较大的文本文件，用Less等操作命令可以顺利地查看。" )
addElement(dom, root,"p" ,type="type1" ,desc="(2) 报告文件含有SVG格式的图片文件，SVG是矢量化的图片文件，可以随意放大而不失真。要查看SVG格式的文件，请先安装SVG插件。" )

#参考文献
reference = 5
ref_list=addElement(dom, root,"ref_list", name="参考文献", type="type1", desc="")
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


'''
if not os.path.exists(analysis_dir + "/Needed_Data"):
	os.mkdir(analysis_dir + "/Needed_Data")
else:
	os.system("rm -r " + analysis_dir + "/Needed_Data")
	os.mkdir(analysis_dir + "/Needed_Data")
	
os.chdir(analysis_dir)	
os.system("zip -r ./Needed_Data/biomarker_htmlReport.zip Web_Report")
os.system("zip -r ./Needed_Data/biomarker_Web_Report.zip configtest.xml Up_load ") 
'''

os.system('mv '+analysis_dir+'/configtest.xml '+analysis_dir + '/Web_Report')
os.system('mv '+analysis_dir+'/Up_load '+analysis_dir + '/Web_Report')

if os.path.exists(analysis_dir + '/Web_Report/index.html'):
	if(delete == 'Y'):
		os.system('rm -r ' + analysis_dir + '/Analysis')
		os.system('rm -r ' + analysis_dir + '/DataAssess')
		os.system('rm -r ' + analysis_dir + '/Diff_analysis')


