'''


Author:
       huangls
Version:
        1.0;2015-09-22
'''
import numpy as np
import sys, os, argparse, os.path,re,math
import matplotlib as mpl
import threading
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import matplotlib.transforms as transform
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerPatch
from pylab import *
mpl.rcParams['interactive']=False
from scipy import interpolate
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline
import multiprocessing
from scipy.interpolate import spline
import matplotlib.transforms as transforms
################################################################################
parser = argparse.ArgumentParser(description='This script is used to ')

parser.add_argument('-i','--input',help='Please input indel stat file',required=True)
parser.add_argument('-o','--out_dir',help='Please input complete out_put directory path',default = os.getcwd(),required=False)
#parser.add_argument('-u','--up_cut',type=int,default=10,help='cut to stat',required=False)
parser.add_argument('-l','--maxline',type=int,default=50,help='number of gene to draw,default 50',required=False)
parser.add_argument('-W','--width',required=False,default=20,type=int,help='specify the width of pic ,default is 20')
parser.add_argument('-H','--height',required=False,default=15,type=int,help='specify the height of pic ,default is 15')

parser.add_argument('-n','--name',default ='demo',required=False,help='Please specify the output demo')
#################################################################################
#Below scripts serve to parse the command line parameters
#################################################################################
args = parser.parse_args()
try:
    assert args.input
except AssertionError:
    print 'Please specify the parameter "--input"'
    sys.exit(1)
try:
    assert args.name
except AssertionError:
    print 'Please specify the parameter "--name"'
    sys.exit(1)
################################################################################
dout=''
if os.path.exists(args.out_dir):
    dout=os.path.abspath(args.out_dir)
else:
    os.mkdir(args.out_dir)
    dout=os.path.abspath(args.out_dir)
args.input=os.path.abspath(args.input)

cmaps = [('Sequential', ['Blues', 'BuGn', 'BuPu',
'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd',
'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu',
'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd']),
('Sequential (2)', ['afmhot', 'autumn', 'bone', 'cool', 'copper',
'gist_heat', 'gray', 'hot', 'pink',
'spring', 'summer', 'winter']),
('Diverging', ['BrBG', 'bwr', 'coolwarm', 'PiYG', 'PRGn', 'PuOr',
'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral',
'seismic']),
('Qualitative', ['Accent', 'Dark2', 'Paired', 'Pastel1',
'Pastel2', 'Set1', 'Set2', 'Set3']),
('Miscellaneous', ['gist_earth', 'terrain', 'ocean', 'gist_stern',
'brg', 'CMRmap', 'cubehelix',
'gnuplot', 'gnuplot2', 'gist_ncar',
'nipy_spectral', 'jet', 'rainbow',
'gist_rainbow', 'hsv', 'flag', 'prism'])]

class myColor():
    '''define my own color object ,color was base on colorbrewer (http://colorbrewer2.org/)'''
    def __init__(self,name,*args,**kwargs):
        self.color={'Paired':['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'],
                    'Set1':['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999'],
                    'Set3':['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9'],
                    'Set2':['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494','#b3b3b3'],
                    'Set1Paired':['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999','#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']
                    }
        self.colorName=name
        if len(args)!=0:
            self.colorNumber=args[0]
        else:
            self.colorNumber=len(self.color[name])
        self.instanseColor=self.color[name][:self.colorNumber]
    def __getitem__(self,key):
        return self.instanseColor[key]
    

def plot_driver_gene(d,col,gene_num,num_variant,q,w,h,fout):
    nullfmt   = NullFormatter()         # no labels
    fontsize=20
    

    #rect = ax1.patch
    #rect.set_facecolor('lightslategray')

#    if n>21:
#        cm = plt.get_cmap("hsv")
#        mycol=[cm(float(i)/(len(data))) for i in xrange(n)]
    # definitions for the axes
    left, width1, width2,width3,right= 0.045, 0.15,0.55,0.15,0.05
    bottom, height = 0.1, 0.8
    gap1 = 0.05
    gap2 = 0.005
    rect_left = [left, bottom, width1, height]
    rect_middle = [left+width1+gap1, bottom, width2, height]
    rect_right = [left+width1+width2+gap1+gap2, bottom, width3, height]

    fig=plt.figure(1, figsize=(w,h))
    axLeft = plt.axes(rect_left)
    axMiddle = plt.axes(rect_middle)
    axRight=plt.axes(rect_right)
    
    rect = axMiddle.patch
    rect.set_facecolor('grey')
    rect.set_alpha(0.01)
    barWidth=1./gene_num
    #print n,barWidth
    rects=axLeft.barh(np.linspace(0, 1-barWidth, gene_num),num_variant,barWidth,color="grey",edgecolor='none')
    rects=axRight.barh(np.linspace(0, 1-barWidth, gene_num),q,barWidth,color="grey",edgecolor='none')
    
    
    
    legendPatch=[]
    legendLebel=[]
    leg={}
    patches=[]
    sample_num=len(d[0])-2
    x=np.linspace(0,1,sample_num+1)
    y=np.linspace(0,1,gene_num+1)
    rw=x[1]-x[0]
    rh=y[1]-y[0]
    label_gene=[]
    print col
    for i,v in enumerate(d):
        label_gene.append(v[0])
        for j,k in enumerate(v[2:]):
            #print k
            if not k ==".":
                k=k.replace("_variant", "")
                #print k
                p=mpatches.Rectangle((x[j], y[i]), rw, rh,lw="none",transform=axMiddle.transAxes, clip_on=False,edgecolor="none",facecolor=col[k])
                patches.append(p)
                if not leg.has_key(k):
                    legendPatch.append(p)
                    legendLebel.append(k)
                    leg[k]=1
                
                #wedge = mpatches.Wedge(grid[1], r2, r0, r0+v, ec="none",fc=col2[i])
                
               
    for i in patches:
        axMiddle.add_patch(i)
        
#################################################
    #transx = transform.blended_transform_factory(axMiddle.transData,axMiddle.transAxes)
    transy = transform.blended_transform_factory(axMiddle.transAxes,axMiddle.transData)
    hl=[]
    for i in y[1:-1]:
        L=mlines.Line2D([0,1], [i,i], color='black', marker=None,markersize=None, linestyle='-',transform=transy)
        hl.append(L)



    for i in hl:
        axMiddle.add_line(i)


    axMiddle.legend(legendPatch,legendLebel,bbox_to_anchor=(0.,1.01 , 1., .05),loc='upper center', shadow=False, fontsize='large',framealpha=0.0,ncol=7)
#    axLeft.xaxis.set_major_formatter(nullfmt)
#    axRight.xaxis.set_major_formatter(nullfmt)
    axLeft.yaxis.set_major_formatter(nullfmt)
    axRight.yaxis.set_major_formatter(nullfmt)
    
    axLeft.xaxis.set_ticks_position('bottom')
    axRight.xaxis.set_ticks_position('bottom')
    axLeft.yaxis.set_ticks_position('right')
    axRight.yaxis.set_ticks_position('left')
    axLeft.spines['left'].set_visible(False)
    axRight.spines['right'].set_visible(False)
    
    axMiddle.xaxis.set_ticks_position("none")
    axMiddle.yaxis.set_ticks_position("none")
    axMiddle.xaxis.set_ticks_position("none")
    axMiddle.yaxis.set_ticks_position("none")
    axMiddle.xaxis.set_major_formatter(nullfmt)
    axMiddle.yaxis.set_major_formatter(nullfmt)
    #axDown.set_ylim(0,yl[1])
    
    
    
    axLeft.set_xlim(axLeft.get_xlim()[::-1])
    
    
    
    axLeft.spines['left'].set_visible(False)
    axLeft.spines['top'].set_visible(False)
    axRight.spines['right'].set_visible(False)
    axRight.spines['top'].set_visible(False)
    
#    axUpLevel.spines['top'].set_visible(False)
#    axUp.spines['right'].set_visible(False)
#    axUp.spines['top'].set_visible(False)
#    axUp.yaxis.set_ticks_position('left')

    
    #WaxDown.spines['left'].set_visible(False)
#    axDown.spines['bottom'].set_visible(False)
    
#format axes
#    axUpLevel.yaxis.set_ticks_position('right')
#    axUpLevel.xaxis.set_ticks_position('bottom')
    for tick in axLeft.yaxis.get_major_ticks():
        tick.label1On = False
        tick.label2On = True
        #tick.label2.set_color('green')



    fig.suptitle('Driver gene mutation distribution',fontsize='xx-large')
    
    axLeft.set_xlabel("variant number",fontsize='x-large')
    axRight.set_xlabel("q-value",fontsize='x-large')
    at=np.linspace(0, 1-barWidth, gene_num)



    #print at,lebel
    axLeft.set_yticks(at+barWidth/2)
    axRight.set_yticks(at+barWidth/2)
    trans = transforms.blended_transform_factory(fig.transFigure,axLeft.transData)
    axLeft.set_yticklabels( label_gene,verticalalignment='center',horizontalalignment='center',x=left+width1+gap1/2,transform=trans )
    
    
#    axLeft.legend(legendPatch,legendLebel,bbox_to_anchor=(left,bottom+height+0.002,width*2+gap,0.049),
#                  loc='upper center', shadow=False, fontsize='large',framealpha=0.0,ncol=15,bbox_transform=fig.transFigure)
    plt.savefig(fout+'.pdf',dpi=600)
    plt.savefig(fout+'.png',dpi=600)
    plt.close()

def read_data(input,geneNum=100):
    f=open(input,'r')
    mycol=myColor('Set1Paired')
    d=[]
    col={}
    gene_num=0
    col_num=0
    q=[]
    num_variant=[]
    for i in f:
        if 'GENE_NAME' in i:
            continue
        i=i.strip()
        tmp=re.split('\s+', i)
        if float(tmp[1])==0:
            q.append(20.)
        else:
            q.append(-math.log10(float(tmp[1])))
        d.append(tmp)
        n=0
        for j in tmp[2:]:
            if j==".":continue
            n+=1
            j=j.replace("_variant", "")
            if col.has_key(j):
                continue
            else:
                col[j]=mycol[col_num]
                col_num+=1
        num_variant.append(n)
        
        gene_num+=1
        if gene_num>geneNum:break
    f.close()
    return d,col,gene_num,num_variant,q

d,col,gene_num,num_variant,q=read_data(args.input,args.maxline)

plot_driver_gene(d,col,gene_num,num_variant,q,args.width,args.height,dout+'/'+args.name)
