'''


Author:
       huangls
Version:
        1.0;2015-08-10
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
parser.add_argument('-u','--up_cut',type=int,default=10,help='cut to stat',required=False)
parser.add_argument('-l','--low_cut',type=int,default=-10,help='cut to stat',required=False)
parser.add_argument('-W','--width',required=False,default=13,type=int,help='specify the width of pic ,default is 13')
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
    

def plot_indel(d,n,w,h,fout,up,down):
    nullfmt   = NullFormatter()         # no labels
    fontsize=20
    

    
    mycol=myColor('Set1Paired')
    if n>21:
        cm = plt.get_cmap("hsv")
        mycol=[cm(float(i)/n) for i in xrange(n)]
    # definitions for the axes
    left, width = 0.07, 0.41
    bottom, height = 0.1, 0.8
    gap = 0.04
    
    rect_left = [left, bottom, width, height]
    rect_right = [left+width+gap, bottom, width, height]

    fig=plt.figure(1, figsize=(w,h))
    axLeft = plt.axes(rect_left)
    axRight=plt.axes(rect_right)
    barWidth=1./(n+1)
    #print n,barWidth
    legendPatch=[]
    legendLebel=[]
    step=0
    c=0
    for i in sorted(d['CDS'].keys()):
        x=[]
        y=[]
        
        
        genome=d['Genome'][i]
        cds=d['CDS'][i]
        for k,v in sorted(sorted(cds.items(), key=lambda cds:cds[0])):
            if k >0:
                y.append(k-1)
                x.append(v)
            else:
                y.append(k)
                x.append(v)
        y=np.array(y)
        rects=axLeft.barh(y+step,x,barWidth,color=mycol[c],edgecolor='none')
        
        legendPatch.append(rects[0])
        legendLebel.append(i)
        x=[]
        y=[]
        for k,v in sorted(sorted(genome.items(), key=lambda genome:genome[0])):
            if k >0:
                y.append(k-1)
                x.append(v)
            else:
                y.append(k)
                x.append(v)
        y=np.array(y)
        rects=axRight.barh(y+step,x,barWidth,color=mycol[c],edgecolor='none')
        
        c+=1
        step+=barWidth
    
#    axLeft.xaxis.set_major_formatter(nullfmt)
#    axRight.xaxis.set_major_formatter(nullfmt)
    axLeft.yaxis.set_major_formatter(nullfmt)
    axRight.yaxis.set_major_formatter(nullfmt)
    
    axLeft.xaxis.set_ticks_position('bottom')
    axRight.xaxis.set_ticks_position('bottom')
    axLeft.yaxis.set_ticks_position('right')
    axRight.yaxis.set_ticks_position('left')

    

    

    
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



    fig.suptitle('Indel Length Distribution',fontsize='xx-large')
    
    axLeft.set_xlabel("CDS",fontsize='x-large')
    axRight.set_xlabel("Genome",fontsize='x-large')
    at=np.arange(down, up)

    lebel=[]
    #print at
    for i,v in enumerate(at):
        if i==(len(at)-1):
            lebel.append(">=%s"%(v+1))
            continue
        if v >=0:
            lebel.append("%s"%(v+1))
            continue
        if i==0:
            lebel.append("<=%s"%(v))
            continue
        lebel.append("%s"%(v))

    #print at,lebel
    axLeft.set_yticks(at+barWidth*n/2)
    axRight.set_yticks(at+barWidth*n/2)
    trans = transforms.blended_transform_factory(fig.transFigure,axLeft.transData)
    axLeft.set_yticklabels( lebel,verticalalignment='center',horizontalalignment='center',x=left+width+gap/2,transform=trans )
    
    if(len(legendLebel)<20):
        axLeft.legend(legendPatch,legendLebel,bbox_to_anchor=(left,bottom+height+0.01,width*2+gap,0.049),
                  loc='upper center', shadow=False, fontsize='large',framealpha=0.0,ncol=10,bbox_transform=fig.transFigure)
    elif(len(legendLebel)<40):
        axLeft.legend(legendPatch,legendLebel,bbox_to_anchor=(left,bottom+height+0.01,width*2+gap,0.049),
                  loc='upper center', shadow=False, fontsize='small',framealpha=0.0,ncol=13,bbox_transform=fig.transFigure)
    else:
        axLeft.legend(legendPatch,legendLebel,bbox_to_anchor=(left,bottom+height+0.01,width*2+gap,0.049),
                  loc='upper center', shadow=False, fontsize='xx-small',framealpha=0.0,ncol=20,bbox_transform=fig.transFigure)
    plt.savefig(fout+'.pdf',dpi=600)
    plt.savefig(fout+'.png',dpi=600)
    plt.close()

def read_data(input,up,down):
    f=open(input,'r')
    d={'CDS':{},'Genome':{}}
    sample_num=0
    for i in f:
        if '#' in i:
            continue
        i=i.strip()
        tmp=re.split('\s+', i)
        if d[tmp[1]].has_key(tmp[3]):
            if int(tmp[0]) >= up:
                d[tmp[1]][tmp[3]][up]+=int(tmp[2])
                continue
            if int(tmp[0]) <= down:
                d[tmp[1]][tmp[3]][down]+=int(tmp[2])
                continue
            d[tmp[1]][tmp[3]][int(tmp[0])]=int(tmp[2])
        else:
            sample_num+=1
            d[tmp[1]][tmp[3]]={}
            d[tmp[1]][tmp[3]][up]=0
            d[tmp[1]][tmp[3]][down]=0
            if tmp[0] >= up:
                d[tmp[1]][tmp[3]][up]+=int(tmp[2])
                continue
            if int(tmp[0]) <= down:
                d[tmp[1]][tmp[3]][down]+=int(tmp[2])
                continue
            d[tmp[1]][tmp[3]][int(tmp[0])]=int(tmp[2])
    f.close()
    #print d
    sample_num=len(d['CDS'].keys())
    return d,sample_num
d={}
d,num=read_data(args.input,args.up_cut,args.low_cut)

plot_indel(d,num,args.width,args.height,dout+'/'+args.name,args.up_cut,args.low_cut)
