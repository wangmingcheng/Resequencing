'''
This script was used to plot a connected zoom pie in one picture.

Author:
       haungls
Version:
        1.0;2015-7-17
'''





import sys, os, argparse, glob, os.path
#from pylab import *
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import math
import re
import operator
#sys.path.append('/share/bioCloud/huangls/02.package/python/svgwrite-1.1.6')
#import svgwrite
mpl.rcParams['interactive']=False


import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
#from svgwrite import cm, mm 
#from cloud import join
#from casuarius import required
#from enstaller.config import default









class myColor():
    '''define my own color object ,color was base on colorbrewer (http://colorbrewer2.org/)'''
    def __init__(self,name,*args,**kwargs):
        self.color={'Paired':['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'],
                    'Set1':['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999'],
                    'Set3':["#8DD3C7", "#FFFFB3" ,"#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69" ,"#FCCDE5" ,"#D9D9D9" ,"#BC80BD", "#CCEBC5" ,"#FFED6F",'#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494','#b3b3b3'],
                    'Set2':['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494','#b3b3b3'],
                    'Pastel1':["#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC" ,"#E5D8BD", "#FDDAEC" ,"#F2F2F2"],
                    'Accent':["#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17","#666666","#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC" ,"#E5D8BD", "#FDDAEC" ,"#F2F2F2"]}
        self.colorName=name
        if len(args)!=0:
            self.colorNumber=args[0]
        else:
            self.colorNumber=len(self.color[name])
        self.instanseColor=self.color[name][:self.colorNumber]
    def __getitem__(self,key):
        return self.instanseColor[key]
    def getCol(self,N):
        return self.instanseColor[0:N]
    def getAllColor(self,N):
        col=[]
        for k,v in self.color.items():
            for i in v:
                col.append(i)
        return col[0:N]

def get_key(info):
    d={}
    item=re.split(';',info);
    for i in item:
        if '=' in i:
            k,v=re.split('=',i)
            d[k]=v
        else:
            d[i]=True
    return d
  



def filter_snp(snp_file,outfile):
    f=open(snp_file,'r')
    fout=open(outfile,'w')
    snp_stat={'unknown':0}
    exonic_stat={'unknown':0}
    
    for i in f:
        i=i.strip()
        if i[0]=='#':
            fout.write(i+'\n')
            continue
        line=re.findall('\./\.',i)
        if line:continue
        
        tmp=re.split("\t", i)
        d=get_key(tmp[7])
        
        #\x3b
        
        a=re.sub('\\\\x3b.*$',"",d['Func.refGene'])
        d['Func.refGene']=a
        if d['Func.refGene']=='.':
            snp_stat['unknown']+=1
        else:
            if snp_stat.has_key(d['Func.refGene']):
                snp_stat[d['Func.refGene']]+=1
            else:
                snp_stat[d['Func.refGene']]=1
        
        if d['Func.refGene']=='exonic':
            if d['ExonicFunc.refGene']=='.':
                exonic_stat['unknown']+=1
            else:
                if exonic_stat.has_key(d['ExonicFunc.refGene']):
                    exonic_stat[d['ExonicFunc.refGene']]+=1
                else:
                    exonic_stat[d['ExonicFunc.refGene']]=1
        
        #filter SNP 
        if not re.match(r'exonic|splicing|\./',d['Func.refGene']):continue
        if not re.match(r'nonsynonymous_SNV|\.',d['ExonicFunc.refGene']):continue
        
        if not d['esp6500siv2_all'] == "." :
            if float(d['esp6500siv2_all']) >0.05 and (1.-float(d['esp6500siv2_all'])) >0.05 :continue
        if not d['1000g2014oct_all'] == "." :
            if float(d['1000g2014oct_all']) >0.05 and (1.-float(d['1000g2014oct_all'])) >0.05 :continue
        fout.write(i+'\n')
    f.close()
    fout.close()
    
    if snp_stat['unknown']==0:
        del snp_stat['unknown']
    if exonic_stat['unknown']==0:
        del exonic_stat['unknown']
        
    fout=open(outfile+'.variation.stat','w')
    for k,v in snp_stat.items():
        fout.write('%s\t%s\n'%(k,v))
    fout.close()
    
    fout=open(outfile+'.exonic.stat','w')
    for k,v in exonic_stat.items():
        fout.write('%s\t%s\n'%(k,v))
    fout.close()
    
    return snp_stat,exonic_stat

def label(xy, text):
    y = xy[1] - 0.15 # shift y-value for label so that it's below the artist
    plt.text(xy[0], y, text, ha="center", family='sans-serif', size=14)
def get_percent(snp_stat,exonic_stat):
    snp={}
    snp_all=0.
    exonic={}
    exonic_all=0.
    for k,v in snp_stat.items():
        snp_all+=v
    for k,v in exonic_stat.items():
        exonic_all+=v
    
    for k,v in snp_stat.items():
        r=v/snp_all
        snp[k]=360.*r
    for k,v in exonic_stat.items():
        r=v/exonic_all
        exonic[k]=360. * r
    return snp,exonic
def sort_dict(d):
    sorted_x = sorted(d.iteritems(), key=operator.itemgetter(1))
    r=[]
    i=0
    while not len(sorted_x)==0:
        e=sorted_x.pop()
        if not len(sorted_x)==0:
            s=sorted_x[i]
            del sorted_x[i]
            r.append(s)
        r.append(e)
    return r
        

def plot_connect_pie(snp_stat,exonic_stat,width,height,outfile,item):
    col1=myColor("Set3")
    col2=myColor("Accent")
    assert snp_stat.has_key(item)
    snp,exonic=get_percent(snp_stat,exonic_stat)
    #plt.figure(1, figsize=(width,height))
    #fig, ax = plt.subplots()
    
    fig = plt.figure(1,figsize=(width,height))
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(0,1), ylim=(0,1))
    
    pi=np.math.pi
    r1=0.16
    r2=0.1
    patches = []
    cy=0.65
    #matplotlib.patches.Wedge(center, r, theta1, theta2, width=None, **kwargs)
    grid=np.array([[0.15,cy],[0.6,cy]])
    
    wedge = mpatches.Wedge(grid[0], r1, -snp[item]/2, snp[item]/2, ec="none",fc=col1[0],lw=1.5)
    patches.append(wedge)
    r0=snp[item]/2
    colorOrderSNP=[]
    colorOrderExonic=[]
    colorOrderSNP.append(['%s:%s(%.2f%%)' %(item,str(int(snp_stat[item])),snp[item]*100/360),col1[0]])
    i=1
    fs1=True
    fs2=True
    fs3=True
    step=15
    if(re.match("SNP", os.path.basename(outfile))):
        step=20
    offstart=10
    offsetX=offstart
    offsetY=offstart
    for k,v in sort_dict(snp):
        if k==item:
            aa=ax.annotate(k, xy=grid[0]+[r1,0],  xycoords='data',
                xytext=(20, 0), textcoords='offset points',
                arrowprops=dict(arrowstyle="->",
                connectionstyle="angle3,angleA=%s,angleB=%s"%(0,90)))
            print '%s:%s,%s\n'%(k,aa.get_position()[0],aa.get_position()[1])
            
            continue
        wedge = mpatches.Wedge(grid[0], r1, r0, r0+v, ec="none",fc=col1[i])
        patches.append(wedge)

        
        x=grid[0][0]+r1*np.math.cos(2*pi*(r0+v/2)/360)
        y=grid[0][1]+r1*np.math.sin(2*pi*(r0+v/2)/360)
        if(r0+v/2) >90 and fs1:
            offsetY=offstart
            offsetX=-offstart
            fs1=False
        if(r0+v/2) >180 and fs2:
            offsetY=-offstart
            offsetX=-offstart
            fs2=False
        if(r0+v/2) >270 and fs3:
            offsetY=-offstart
            offsetX=offstart
            fs3=False

        if (r0+v/2) >0 and (r0+v/2) <= 180:
            #offsetX=offsetX*np.math.cos(2*pi*(r0+v/2)/360)
            #offsetY=offsetY*np.math.sin(2*pi*(r0+v/2)/360)
            if (r0+v/2) >0 and (r0+v/2) <= 90:
                offsetX+=step
                offsetY+=step
                aa=ax.annotate(k, xy=(x,y),  xycoords='data',
                xytext=(offsetX, offsetY), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),horizontalalignment='left', verticalalignment='bottom'
                )
            else:
                offsetX-=step
                offsetY+=step
                aa=ax.annotate(k, xy=(x,y),  xycoords='data',
                xytext=(offsetX, offsetY), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),horizontalalignment='right', verticalalignment='bottom'
                )

        if (r0+v/2) >180 and (r0+v/2) <= 360:
            #offsetX=offsetX*np.math.cos(2*pi*(r0+v/2)/360)
            #offsetY=offsetY*np.math.sin(2*pi*(r0+v/2)/360)
            if (r0+v/2) >180 and (r0+v/2) <= 270:
                offsetX-=step
                offsetY-=step
                aa=ax.annotate(k, xy=(x,y),  xycoords='data',
                xytext=(offsetX, offsetY), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),horizontalalignment='right', verticalalignment='top'
                )
            else:
                offsetX+=step
                offsetY-=step
                aa=ax.annotate(k, xy=(x,y),  xycoords='data',
                xytext=(offsetX, offsetY), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),horizontalalignment='left', verticalalignment='top'
                )
        colorOrderSNP.append(["%s:%s(%.2f%%)"%(k,str(int(snp_stat[k])),v*100/360),col1[i]])
        i+=1
        r0+=v
        
        
    i=0
    r0=0
    offsetX=offstart
    offsetY=offstart
    fs1=True
    fs2=True
    fs3=True
    stepx=10
    if(re.match("SNP", os.path.basename(outfile))):
        step=18
    for k,v in sort_dict(exonic):
        wedge = mpatches.Wedge(grid[1], r2, r0, r0+v, ec="none",fc=col2[i])
        colorOrderExonic.append(["%s:%s(%.2f%%)"%(k,str(int(exonic_stat[k])),v*100/360),col2[i]])
        patches.append(wedge)
        x=grid[1][0]+r2*np.math.cos(2*pi*(r0+v/2)/360)
        y=grid[1][1]+r2*np.math.sin(2*pi*(r0+v/2)/360)
        if(r0+v/2) >90 and fs1:
            offsetY=offstart
            offsetX=-offstart
            fs1=False
        if(r0+v/2) >180 and fs2:
            offsetY=-offstart
            offsetX=-offstart
            fs2=False
        if(r0+v/2) >270 and fs3:
            offsetY=-offstart
            offsetX=offstart
            fs3=False

        
            
        
        if (r0+v/2) >0 and (r0+v/2) <= 180:
            #offsetX=offsetX*np.math.cos(2*pi*(r0+v/2)/360)
            #offsetY=offsetY*np.math.sin(2*pi*(r0+v/2)/360)
            if (r0+v/2) >0 and (r0+v/2) <= 90:
                offsetX+=stepx
                offsetY+=step
                aa=ax.annotate(k, xy=(x,y),  xycoords='data',
                xytext=(offsetX, offsetY), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),horizontalalignment='left', verticalalignment='bottom'
                )
            else:
                offsetX-=stepx
                offsetY+=step
                aa=ax.annotate(k, xy=(x,y),  xycoords='data',
                xytext=(offsetX, offsetY), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),horizontalalignment='right', verticalalignment='bottom'
                )
                print '%s:%s,%s\n'%(k,aa.get_position()[0],aa.get_position()[1])
#  
        if (r0+v/2) >180 and (r0+v/2) <= 360:
            #offsetX=offsetX*np.math.cos(2*pi*(r0+v/2)/360)
            #offsetY=offsetY*np.math.sin(2*pi*(r0+v/2)/360)
            if (r0+v/2) >180 and (r0+v/2) <= 270:
                offsetX-=stepx
                offsetY-=step
                aa=ax.annotate(k, xy=(x,y),  xycoords='data',
                xytext=(offsetX, offsetY), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),horizontalalignment='right', verticalalignment='top'
                )
            else:
                offsetX+=stepx
                offsetY-=step
                aa=ax.annotate(k, xy=(x,y),  xycoords='data',
                xytext=(offsetX, offsetY), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),horizontalalignment='left', verticalalignment='top'
                )

        r0+=v
        i+=1
        
        
    if snp[item]/2<90:
        x1=grid[0][0]+r1*np.math.cos(2*pi*snp[item]/2/360)
        y1=grid[0][1]+r1*np.math.sin(2*pi*snp[item]/2/360)
        l=np.math.fabs(grid[0][0]-grid[1][0])
        a=np.math.pow(l**2-r2**2,0.5)
        
        x2=grid[1][0]-r2*r2/l+0.005
        y2=r2*a/l+grid[1][1]
        line1 = mlines.Line2D([x1,x2], [y1,y2], lw=1.5, alpha=1,c='black')
        ax.add_line(line1)
        line2 = mlines.Line2D([x1,x2], [2*cy-y1,2*cy-y2], lw=1.5, alpha=1,c='black')
        ax.add_line(line2)
    else:
        l=np.math.fabs(grid[0][0]-grid[1][0])
        a=np.math.pow(l**2-(r1-r2)**2,0.5)
        
        x1=grid[0][0]+r1*(r1-r2)/l
        y1=grid[0][1]+r1*a/l

        
        x2=grid[1][0]-r2*r2/l
        y2=r2*a/l+grid[1][1]
        line1 = mlines.Line2D([x1,x2], [y1,y2], lw=1.5, alpha=1,c='black')
        ax.add_line(line1)
        line2 = mlines.Line2D([x1,x2], [2*cy-y1,2*cy-y2], lw=1.5, alpha=1,c='black')
        ax.add_line(line2)
#add legend: aaaaaaaaaaaaaaaaa
    snpRectangles=[]
    snpLegendLebel=[]
    for i in colorOrderSNP:
        m=mpatches.Rectangle([1,1], 0.05, 0.1, ec="none",fc=i[1])
        snpRectangles.append(m)
        snpLegendLebel.append(i[0])
    for i in colorOrderExonic:
        m=mpatches.Rectangle([1,1], 0.05, 0.1, ec="none",fc=i[1])
        snpRectangles.append(m)
        snpLegendLebel.append(i[0]) 
    #ax.legend(exonicRectangles,exonicLegendLebel,loc='lower right', shadow=False, fontsize='x-small',framealpha=0.0,ncol=2)
    ax.legend(snpRectangles,snpLegendLebel,loc='lower center', shadow=False, fontsize='large',framealpha=0.0,ncol=3)
    #colors = np.linspace(0, 1, len(patches))
    #collection = PatchCollection(patches, cmap=plt.cm.hsv, alpha=0.3)
    #collection.set_array(np.array(colors))
    for i in patches:
        ax.add_patch(i)
    #ax.set_xlim(0,1)
    #ax.set_ylim(0,1)
    
#connect two circle


    
    
    
    plt.subplots_adjust(left=0, right=0.8, bottom=0, top=1)
    #plt.xlim(0,1)
    #plt.ylim(0,1)
    plt.axis('equal')
    plt.axis('off')
    plt.savefig(outfile+'.pdf')
    plt.savefig(outfile+'.png',dpi=300)
    plt.close()

parser = argparse.ArgumentParser(description='This script was used to plot a connected zoom pie in one picture.')
parser.add_argument('-i', '--input', dest='input', required=True, help='input main SNPEFF ann  file')
parser.add_argument('-c', '--connectLebel', dest='connectLebel', required=False,default="CDS", help='specify which item to connect,default CDS')
parser.add_argument('-o','--outDir',dest='outDir',required=False,default=os.getcwd(),help='specify the output file dir,default is current dir')
parser.add_argument('-p','--outFilePrefix',dest='outFilePrefix',required=False,default='demo',help='specify the output file prefix,demo')
parser.add_argument('-W','--width',required=False,default=20,type=int,help='specify the width column index,default is 15')
parser.add_argument('-H','--height',required=False,default=15,type=int,help='specify the height column index,default is 15')

args = parser.parse_args()




if not os.path.exists(args.outDir): os.mkdir(args.outDir)
args.outDir=os.path.abspath(args.outDir)


#snp_stat,exonic_stat=filter_snp(args.input,args.outDir+'/'+args.outFilePrefix+'.vcf')

all={}
sample={}
f=open(args.input,'r')
for i in f:
    i=i.strip()
    tmp=re.split('\s+',i)
    if '#' in i:
        for n,v in enumerate(tmp):
             if n>1:
                 sample[n]=v
        continue
    if 'CDS' in tmp[0]:
        for n,v in enumerate(tmp):
            if n>1:
                if all.has_key(sample[n]):
                    all[sample[n]]['cds'][tmp[1]]=int(tmp[n])
                    all[sample[n]]['snp']['CDS']+=int(tmp[n])
                else:
                    all[sample[n]]={'snp':{'CDS':0},'cds':{}}
                    all[sample[n]]['cds'][tmp[1]]=int(tmp[n])
                    all[sample[n]]['snp']['CDS']+=int(tmp[n])
    else:
        for n,v in enumerate(tmp):
            if n>1:
                if all.has_key(sample[n]):
        
                    all[sample[n]]['snp'][tmp[1]]=int(tmp[n])
                else:
                    all[sample[n]]={'snp':{'CDS':0},'cds':{}}
                    all[sample[n]]['snp'][tmp[1]]=int(tmp[n])

    
f.close()

for i in all.keys():
    plot_connect_pie(all[i]['snp'],all[i]['cds'],args.width,args.height,args.outDir+'/'+args.outFilePrefix+"."+i+'.pie',args.connectLebel)

