#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript

library('getopt');

#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'help' , 'hl', 0, "logical",
  'infile' , 'i', 1, "character",
  'outfile' , 'o', 0, "character",
  'height' , 'h', 0, "integer",
  'width' , 'w', 0, "integer",
  'col','c',0,"character",
  'sample','s',1,"character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);


# define usage function
print_usage=function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("discription:
        plot reads average error rate from quality file.see test file.
      Usage example: 
      Rscript /share/nas1/wangm/yf/plot/quality/quality_bar.R --infile /share/nas1/wangm/yf/plot/quality/T08.quality --outfile /share/nas1/wangm/yf/plot/quality/bar.png 
      Options: 
      --help  	-hl 	NULL 		get this help
      --infile 	-i 	character 	the input file [forced]
      --outfile -o 	character 	the filename for output graph [optional]
      --height 	-h 	integer 	the height of graph,the unit is mm [optional, default: 120]
      --width 	-w 	integer 	the width of graph, the unit is mm [optional, default: 180]
      --col 		-c 	character [optional, default: blue]
      \n")
  q(status=1);
}


# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }


# check non-null args
if ( is.null(opt$infile) )	{ print_usage(spec) }else{opt$infile<-gsub("\\\\",replacement = "/",x = opt$infile)}


if ( is.null(opt$outfile) )  { opt$outfile =paste(c(getwd(),"/","quality_bar.png"),collapse = "")}
if ( is.null(opt$height ) )  	{ opt$height = 100 }
if ( is.null(opt$width ) )		{ opt$width = 100 }
if ( is.null(opt$col ) )  { opt$col ="blue" }
if ( is.null(opt$sample ) )  { opt$sample ="Demo" }
#-----------------------------------------------------------------
# reading data
#-----------------------------------------------------------------
# reading data

data=read.table(opt$infile,comment.char ="",header=T,check=F,row.names=1)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

#convert Q value to error p value
p=10^(-(0:(ncol(data)-1))/10)

#calculate average error p value for every cycle
err=apply(data,1,function(x){
  mean(x*p/100)
})


#y labels
y.at=seq(0,max(err),by=0.0005)
if(max(err)>max(y.at)){
  y.at=c(y.at,max(y.at)+0.0005)
}
labels=y.at*100
labels[2]="0.05"

## png(filename=paste(opt$outfile,".png",sep=""),width=opt$width,height=opt$height,units="mm",res=600)
## #windowsFonts(Arial=windowsFont("TT Arial")) 
## par(mar=c(4,3.5,3,1),mgp=c(1.5,0.1,0),tck=-0.005,
##     #family="Arial",
##     font.lab=2,font.axis=1,cex=0.6)
## out=barplot(err,col=opt$col,border=opt$col,axes =F,axisnames=F,
##             xlab="",ylab="Average Error Rate(%)",
##             main=paste("Reads Average Error Rate(",opt$sample,")",sep=""),ylim=c(0,max(y.at)))
## out=cbind(out,rep(1:(nrow(data)/2),times=2))
## read=max(out[,1])/2
## segments(read,0,read,5,col="darkgray",lwd=0.5,lty=2)
## pos=seq(0,nrow(data)/2,50)
## pos[1]=1
## at=out[,2]%in%pos
## axis(1,labels=rep(pos,times=2),at=out[at,1])
## axis(2,at=y.at,labels=labels)
## mtext("Read1",side=1,line=1.5,at=read/2,font=2,cex=0.6)
## mtext("Read2",side=1,line=1.5,at=read*1.5,font=2,cex=0.6)
## box()
## dev.off()

mycol<-brewer.pal(4, "Set1")

mdf<-data.frame(err=err,pos=1:length(err))
head(mdf)
p <- ggplot(mdf, aes(xmin = pos-1, xmax = pos , ymin = 0, ymax = err*100)) +
		geom_rect(fill=mycol[3])+
		geom_vline(xintercept = (length(mdf$pos) )/2, colour="grey", linetype = "longdash")+
		coord_cartesian(ylim=c(0, max(err*100)+max(err*100)*0.1))+
		xlab("Position along reads") + ylab("Average Error Rate(%)")+
		labs(title=paste("Reads Average Error Rate(",opt$sample,")",sep=""))


p <- p +theme_bw()+ theme(  
				panel.grid=element_blank(), 
				axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),panel.border=element_rect(colour = "black"))+
		#scale_colour_manual(values = c(mycol[1],mycol[2],mycol[3],mycol[4],"grey"))+
		theme(legend.key = element_blank(),legend.title = element_blank())
#
ggsave(filename=paste(opt$outfile,".pdf",sep=""),plot=p)
ggsave(filename=paste(opt$outfile,".png",sep=""),type="cairo-png",plot=p)

