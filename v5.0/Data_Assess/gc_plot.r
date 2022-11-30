
library(getopt)

#+--------------------
# get options
#+--------------------
spec <- matrix(c(
				'help', 'h', 0, "logical", "help",
				'verbose', 'v', 2, "integer", "verbose mode, default [1]",
				'gc', 'i', 1, "character", "input file gc content, forced.",
				'output', 'o', 1, "character", "output png file, forced.",
				'name', 'n', 1, "character", "output sample name, "
				
		
		), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

#+--------------------
# check options
#+--------------------
if ( !is.null(opt$help) | is.null(opt$gc) | is.null(opt$output) ) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

library(ggplot2)
library(reshape2)
library(RColorBrewer)

mycol<-brewer.pal(4, "Set1")
df<-read.table(opt$gc,sep="\t",head=F)
#df<-as.data.frame(cbind(data[,4],data[,7],data[,10],data[,13],data[,16]))
colnames(df)<-c("pos","A", "C"  ,  "G" ,  "T"  ,  "N")
mdf<-melt(df,measure=c("A", "C"  ,  "G" ,  "T"  ,  "N"))
colnames(mdf)<-c("pos","type","percent")
p <- ggplot(mdf, aes(x=pos, y=percent, group=type))+
		geom_line(aes(colour = type),size=0.5)+
		coord_cartesian(ylim=c(-1,50))+
		xlab("Position along reads") + ylab("Percent of bases")+
		labs(title=paste("Bases content along reads(",opt$name,")",sep=""))


p <- p +theme_bw()+ theme(  
		panel.grid=element_blank(), 
		axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),panel.border=element_rect(colour = "black"))+
		scale_colour_manual(values = c(mycol[1],mycol[2],mycol[3],mycol[4],"grey"))+
		theme(legend.key = element_blank(),legend.title = element_blank())
#
ggsave(filename=paste(opt$output,"/",opt$name,".acgtn.pdf",sep=""),plot=p)
ggsave(filename=paste(opt$output,"/",opt$name,".acgtn.png",sep=""),type="cairo-png",plot=p)
