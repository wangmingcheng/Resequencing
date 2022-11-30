#!/share/nas2/genome/biosoft/R/current/bin/Rscript

#####################################################################

#
# Author: huangls 2015-8-20
#
#
#####################################################################

library(getopt)

#+--------------------
# get options
#+--------------------
spec <- matrix(c(
				'help', 'h', 0, "logical", "help",
				'verbose', 'v', 2, "integer", "verbose mode, default [1]",
				'input', 'i', 1, "character", "input file with preprocessed cytosine coverage depth info, forced.",
				'output', 'o', 1, "character", "output png file, forced.",
				'name', 'n', 1, "character", "output file prefix, SNP."
		), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

#+--------------------
# check options
#+--------------------
if ( !is.null(opt$help) | is.null(opt$input)  ) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}





if ( is.null(opt$output)){opt$output=getwd()}
if ( is.null(opt$name)){opt$name="SNP"}


#+--------------------
# some default options
#+--------------------




#+--------------------
# Main
#+--------------------
freqStat<-function(duration,s=1,e=100,b=1){
	
	breaks = seq(s, e+1, by=b)
	duration.cut = cut(duration, breaks, right=FALSE)
	duration.freq = table(duration.cut) 
	duration.cumfreq = cumsum(duration.freq)
	r=as.numeric(duration.cumfreq)
	return(r)
	
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
	library(grid)
	
	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)
	
	numPlots = length(plots)
	
	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) {
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
				ncol = cols, nrow = ceiling(numPlots/cols))
	}
	
	if (numPlots==1) {
		print(plots[[1]])
		
	} else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		
		# Make each plot, in the correct location
		for (i in 1:numPlots) {
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
							layout.pos.col = matchidx$col))
		}
	}
}

get_data_list<-function(str,names,out_file){
	uniq_all_id=c()
	data_list=list()
	for(i in 1:length(str)){
		#cat(infile,"\n")
		mydata<-str[[i]]
		
		uniq_all_id<-union(as.character(mydata),uniq_all_id)
		data_list[[i]]<-mydata
		i<-i+1
	}
	
	aa=matrix(data = "NA", nrow = length(uniq_all_id), ncol = length(data_list)+1, byrow = FALSE,dimnames = NULL)
	#aa<-as.data.frame(aa)
	
	aa[,1]=as.character(uniq_all_id)
	print(head(aa))
	names(data_list)<-as.character(names)
	for(i in 1:length(data_list)){
		for (j in 1:nrow(aa)){
			if(aa[j,1] %in% as.character(data_list[[i]])){
				aa[j,i+1]=aa[j,1]
			}else{
			}
			
		}
	}
	colnames(aa)=c("all_ID",names)
	write.table(aa,file=out_file,quote=F,sep='\t',row.names=F,col.names=T)
	return(data_list)
}

venn9<-function(IDlist,names,pic_name){
	library(Vennerable)
	data_list<-IDlist
	data<-Venn(data_list)
	infile_num=length(data_list)
	isWeight=FALSE
	if(infile_num>=6){
		isWeight=FALSE
	}
	pdf(file=paste(pic_name,".png",sep=""), height=10, width=10)
	plot(data,doWeight=isWeight)
	dev.off()
	png(file==paste(pic_name,".pdf",sep=""), height=1000, width=1000)
	plot(data,doWeight=isWeight)
	dev.off()
}

venn5<-function(data_list,file_name,main=""){
	library(VennDiagram)
	
	
	pdf(paste(file_name,".pdf",sep=""),h=10,w=10)
	venn.plot <- venn.diagram(
			x = data_list,
			filename = NULL,
			#col = c("red","yellow","green","purple","salmon4"),
			fill = c("deeppink", "lightpink", "peachpuff2","khaki1","darkseagreen2"),
			alpha = 0.50,
			cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
					1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
			cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
			cat.cex = 1.5,
			cat.fontface = "bold",
			margin = 0.05
	);
	grid.draw(venn.plot);
	dev.off()
	png(paste(file_name,".png",sep=""),h=3000,w=3000,res=300)
	venn.plot <- venn.diagram(
			x = data_list,
			filename = NULL,
			#col = c("red","yellow","green","purple","salmon4"),
			fill = c("deeppink", "lightpink", "peachpuff2","khaki1","darkseagreen2"),
			alpha = 0.50,
			cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
					1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
			cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
			cat.cex = 1.5,
			cat.fontface = "bold",
			margin = 0.05
	);
	grid.draw(venn.plot);
	dev.off()
	
}

venn4<-function(data_list,file_name,main=""){
	library(VennDiagram)
	
	
	pdf(paste(file_name,".pdf",sep=""),h=10,w=10)
	venn.plot <- venn.diagram(
			x = data_list,
			filename =NULL,
			#col = c("yellow","green","purple","red"),
			#lty = "dotted",
			lwd = 2,
			fill = c("darkseagreen2", "lightpink", "peachpuff2","khaki1"),
			#fill = c("darkorchid1", "darkorchid1", "darkorchid1", "darkorchid1"),
			alpha = 0.80,
			label.col = c("orange", "black", "darkorchid4", "black", "black", "black",
					"black", "black", "darkblue", "black",
					"black", "black", "black", "darkgreen", "black"),
			cex = 1.5,
			fontfamily = "serif",
			fontface = "bold",
			cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
			cat.cex = 2,
			cat.fontfamily = "serif"
	);
	grid.draw(venn.plot);
	dev.off()
	png(paste(file_name,".png",sep=""),h=3000,w=3000,res=300)
	venn.plot <- venn.diagram(
			x = data_list,
			filename =NULL,
			#col = c("yellow","green","purple","red"),
			#lty = "dotted",
			lwd = 2,
			fill = c("darkseagreen2", "lightpink", "peachpuff2","khaki1"),
			alpha = 0.80,
			label.col = c("orange", "black", "darkorchid4", "black", "black", "black",
					"black", "black", "darkblue", "black",
					"black", "black", "black", "darkgreen", "black"),
			cex = 1.5,
			fontfamily = "serif",
			fontface = "bold",
			cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
			cat.cex = 2,
			cat.fontfamily = "serif"
	);
	grid.draw(venn.plot);
	dev.off()
	
}
venn3<-function(data_list,file_name,main=""){
	library(VennDiagram)
	
	
	pdf(paste(file_name,".pdf",sep=""),h=10,w=10)
	venn.plot <- venn.diagram(
			x =data_list,
			#col = c("yellow","green","purple"),
			filename = NULL,
			output = TRUE,
			#height = 3000,
			#width = 3000,
			resolution = 300,
			compression = 'lzw',
			units = 'px',
			lwd = 2,
			#lty = 'blank',
			fill = c("darkseagreen2","peachpuff2","khaki1"),
			alpha = 0.80,
			cex = 1.5,
			#fontface = "bold",
			fontfamily = "sans",
			cat.cex = 2,
			#cat.fontface = "bold",
			cat.default.pos = "outer",
			cat.col = c("darkblue", "darkgreen", "orange"),
			cat.pos = c(-27, 27, 135),
			cat.dist = c(0.055, 0.055, 0.085),
			cat.fontfamily = "sans",
			rotation = 1
	);
	grid.draw(venn.plot);
	dev.off()
	png(paste(file_name,".png",sep=""),h=3000,w=3000,res=300)
	venn.plot <- venn.diagram(
			x =data_list,
			#col = c("yellow","green","purple"),
			filename = NULL,
			output = TRUE,
			# = 3000,
			#width = 3000,
			resolution = 300,
			compression = 'lzw',
			units = 'px',
			lwd = 2,
			#lty = 'blank',
			fill = c("darkseagreen2","peachpuff2","khaki1"),
			alpha = 0.80,
			cex = 1.5,
			#fontface = "bold",
			fontfamily = "sans",
			cat.cex = 2,
			#cat.fontface = "bold",
			cat.default.pos = "outer",
			cat.col = c("darkblue", "darkgreen", "orange"),
			cat.pos = c(-27, 27, 135),
			cat.dist = c(0.055, 0.055, 0.085),
			cat.fontfamily = "sans",
			rotation = 1
	);
	grid.draw(venn.plot);
	dev.off()
	
}
venn2<-function(data_list,file_name,main=""){
	library(VennDiagram)
	
	
	pdf(paste(file_name,".pdf",sep=""),h=10,w=10)
	venn.plot <- venn.diagram(
			x = data_list,
			filename = NULL,
			#col = c("yellow","green"),
			lwd = 2,
			fill = c("darkseagreen2","khaki1"),
			alpha = 0.75,
			label.col = "black",
			cex = 1.5,
			fontfamily = "serif",
			fontface = "bold",
			cat.col = c("cornflowerblue", "darkorchid1"),
			cat.cex = 2,
			cat.fontfamily = "serif",
			cat.fontface = "bold",
			cat.dist = c(0.03, 0.03),
			cat.pos = c(-20, 14)
	);
	grid.draw(venn.plot);
	dev.off()
	png(paste(file_name,".png",sep=""),h=3000,w=3000,res=300)
	venn.plot <- venn.diagram(
			x = data_list,
			#col = c("yellow","green"),
			filename = NULL,
			lwd = 2,
			fill = c("darkseagreen2","khaki1"),
			alpha = 0.75,
			label.col = "black",
			cex = 1.5,
			fontfamily = "serif",
			fontface = "bold",
			cat.col = c("cornflowerblue", "darkorchid1"),
			cat.cex = 2,
			cat.fontfamily = "serif",
			cat.fontface = "bold",
			cat.dist = c(0.03, 0.03),
			cat.pos = c(-20, 14)
	);
	grid.draw(venn.plot);
	dev.off()
	
}

venn1<-function(data_list,file_name,main=""){
	library(VennDiagram)
	
	
	pdf(paste(file_name,".pdf",sep=""),h=10,w=10)
	venn.plot <- venn.diagram(
			x = data_list,
			#col = "yellow",
			filename = NULL,
			col = "black",
			lwd = 2,
			fontface = "bold",
			fill = "darkseagreen2",
			cat.col = "cornflowerblue",
			alpha = 0.75,
			cex = 1.5,
			cat.cex = 2,
			cat.fontface = "bold",
	);
	grid.draw(venn.plot);
	dev.off()
	png(paste(file_name,".png",sep=""),h=3000,w=3000,res=300)
	venn.plot <- venn.diagram(
			x = data_list,
			#col = "yellow",
			filename = NULL,
			col = "black",
			lwd = 2,
			fontface = "bold",
			fill = "darkseagreen2",
			cat.col = "cornflowerblue",
			alpha = 0.75,
			cex = 1.5,
			cat.cex = 2,
			cat.fontface = "bold",
	);
	grid.draw(venn.plot);
	dev.off()
	
}

## load ggplot2 
library(ggplot2)
library(grid)
library(RColorBrewer)


mycol<-brewer.pal(4, "Set1")
## get data
data <- read.table(opt$input, header=FALSE, sep="\t", comment.char = "#")
colnames(data) <- c("id", "value", "sample")

if (ncol(data) != 3) stop("Input file format error: must be <id> <value> <sample>")

## subset data frame if chromsome in scaffold level and specified the top_n option
depth<-subset(data,id=="depth")
distance<-subset(data,id=="distance")
trans<-subset(data,id=="trans")
snpID<-subset(data,id=="id")
samples<-unique(depth$sample)

#########################################################################################
depthData=data.frame()
for (i in samples){
	d<-subset(depth,sample==i)
	#print(summary(d))
	duration =as.numeric(as.character(d$value))
	
	s=length(duration)
	
	duration.cumfreq=freqStat(duration)
	sd<-data.frame(x=seq(1, 100, by=1),y=as.vector(duration.cumfreq)/s,sample=i)
	depthData<-rbind(depthData,sd)
}



p1<-ggplot(depthData,aes(x=x,y=y*100,color=sample)) + geom_line()+xlab("Depth") + ylab("SNP Fraction(%)")+labs(title="Cumulative SNP depth distribution")

p1 <- p1 +theme_bw()+ theme(  
				panel.grid=element_blank(), 
				axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),panel.border=element_rect(colour = "black"))+
		#scale_colour_manual(values = c(mycol[1],mycol[2],mycol[3],mycol[4],"grey"))+
		theme(legend.key = element_blank(),legend.title = element_blank())
	if(nlevels(depthData$sample)>20) {p1<-p1+guides(col = guide_legend(ncol = 4))}

#png(filename=paste(opt$output,"/",opt$name,".depth.distribution.png",sep=""), height = 2000, width = 2000, res = 300, units = "px")
#print(p)
#dev.off()


##########################################################################
distanceData=data.frame()
for (i in samples){
	d<-subset(distance,sample==i)
	#print(summary(d))
	duration =as.numeric(as.character(d$value))
	
	s=length(duration)
	
	duration.cumfreq=freqStat(duration,e=500)
	sd<-data.frame(x=seq(1, 500, by=1),y=as.vector(duration.cumfreq)/s,sample=i)
	distanceData<-rbind(distanceData,sd)
}
p<-ggplot(distanceData,aes(x=x,y=y*100,color=sample)) + geom_line()+xlab("Neighboring SNP distance(bp)") + ylab("SNP Fraction(%)")+labs(title="Cumulative Neighboring SNP distance distribution")
p <- p +theme_bw()+ theme(  
				panel.grid=element_blank(), 
				axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),panel.border=element_rect(colour = "black"))+
		#scale_colour_manual(values = c(mycol[1],mycol[2],mycol[3],mycol[4],"grey"))+
		theme(legend.key = element_blank(),legend.title = element_blank())
#if(nlevels(distanceData$sample)>20) {p<-p+guides(col = guide_legend(ncol = 4))}
#png(filename=paste(opt$output,"/",opt$name,".quality.distribution.png",sep=""), height = 2000, width = 4000, res = 300, units = "px")

if(nlevels(distanceData$sample)<30) {p<-p+guides(col = guide_legend(ncol = 1))}
if(nlevels(distanceData$sample)<30) {png(filename=paste(opt$output,"/",opt$name,".quality.distribution.png",sep=""), height = 2000, width = 4000, res = 300, units = "px")}

if(nlevels(distanceData$sample)>=30 & nlevels(distanceData$sample)<=100) {p<-p+guides(col = guide_legend(nrow=30,ncol=4))}
if(nlevels(distanceData$sample)>=30 & nlevels(distanceData$sample)<=100) {png(filename=paste(opt$output,"/",opt$name,".quality.distribution.png",sep=""), height = 3000, width = 6000, res = 300, units = "px")}

if(nlevels(distanceData$sample)>100 & nlevels(distanceData$sample)<=180) {p<-p+guides(col = guide_legend(nrow=40,ncol = 5))}
if(nlevels(distanceData$sample)>100 & nlevels(distanceData$sample)<=180) {png(filename=paste(opt$output,"/",opt$name,".quality.distribution.png",sep=""), height = 3500, width = 7000, res = 300, units = "px")}

if(nlevels(distanceData$sample)>180) {p<-p+guides(col = guide_legend(nrow=45,ncol = 5))}
if(nlevels(distanceData$sample)>180) {png(filename=paste(opt$output,"/",opt$name,".quality.distribution.png",sep=""), height = 4000, width = 8000, res = 300, units = "px")}

multiplot(p1, p, cols=2)
dev.off()

##############################################################################

mytrans<-matrix(unlist(strsplit(as.character(trans$value),"-")),ncol=2,byrow=T)

mytrans<-cbind(mytrans,as.character(trans$sample))
mytrans<-data.frame(tran=mytrans[,1],value=as.numeric(mytrans[,2]),sample=mytrans[,3])

p<-ggplot(mytrans, aes(x=tran,y=value, fill=sample)) + geom_bar(position="dodge",stat="identity")+coord_cartesian(ylim=c(0, max(mytrans$value)+max(mytrans$value)/10))+
		xlab("Mutation type") + ylab("SNP number")+labs(title="SNP  Mutation type distribution")
p <- p +theme_bw()+ theme(  
				panel.grid=element_blank(), 
				axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),panel.border=element_rect(colour = "black"))+
		theme(legend.key = element_blank(),legend.title = element_blank(),axis.text.x=element_text(angle=45,hjust=1))
#if(nlevels(mytrans$sample)>20) {p<-p+guides(fill = guide_legend(ncol = 4))}
#png(filename=paste(opt$output,"/",opt$name,".mutation.distribution.png",sep=""), height = 2000, width = 3000, res = 300, units = "px")

if(nlevels(mytrans$sample)<30) {p<-p+guides(fill = guide_legend(ncol = 1))}
if(nlevels(mytrans$sample)<30) {png(filename=paste(opt$output,"/",opt$name,".mutation.distribution.png",sep=""), height = 2000, width = 3000, res = 300, units = "px")}

if(nlevels(mytrans$sample)>=30 & nlevels(mytrans$sample)<=100) {p<-p+guides(fill = guide_legend(nrow=30,ncol =4))}
if(nlevels(mytrans$sample)>=30 & nlevels(mytrans$sample)<=100) {png(filename=paste(opt$output,"/",opt$name,".mutation.distribution.png",sep=""), height = 3000, width = 5000, res = 300, units = "px")}

if(nlevels(mytrans$sample)>100 & nlevels(mytrans$sample)<=180) {p<-p+guides(fill = guide_legend(nrow=40,ncol =5))}
if(nlevels(mytrans$sample)>100 & nlevels(mytrans$sample)<=180) {png(filename=paste(opt$output,"/",opt$name,".mutation.distribution.png",sep=""), height = 3500, width = 6000, res = 300, units = "px")}

if(nlevels(mytrans$sample)>180) {p<-p+guides(fill = guide_legend(nrow=45,ncol = 5))}
if(nlevels(mytrans$sample)>180) {png(filename=paste(opt$output,"/",opt$name,".mutation.distribution.png",sep=""), height = 4000, width = 7000, res = 300, units = "px")}

print(p)
dev.off()

##############################################################################


##########################
names=c()
IDlist=list()
venn_pic_name=paste(opt$output,"/",opt$name,".venn",sep="")
for (i in samples){
	
	d<-subset(snpID,sample==i)
	IDlist[[i]]=d$value
	
}



if(length(samples)<6){
	if(length(samples)==1){venn1(IDlist,venn_pic_name)}
	if(length(samples)==2){venn2(IDlist,venn_pic_name)}
	if(length(samples)==3){venn3(IDlist,venn_pic_name)}
	if(length(samples)==4){venn4(IDlist,venn_pic_name)}
	if(length(samples)==5){venn5(IDlist,venn_pic_name)}
	#venn9(IDlist,venn_pic_name)
}else if(length(samples)<10 & length(samples)>5 ){
	#venn9(IDlist,venn_pic_name)
}
#get_data_list(IDlist,samples,paste(venn_pic_name,".stat.xls",sep=""))

