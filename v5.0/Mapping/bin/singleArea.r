#!/share/nas2/genome/biosoft/R/current/bin/Rscript


# load library
library('getopt');
require(ggplot2)
require(RColorBrewer)


#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
	'help' , 'h', 0, "logical",
	'infile' , 'i', 1, "character",
	'outfile' , 'o', 1, "character",
	'x.col' , 'x', 1, "integer",
	'y.col' , 'y', 1, "integer",
	'height' , 'H', 1, "integer",
	'width' , 'W', 1, "integer",
	'color' , 'c', 1, "character",
	'x.lab' , 'X', 1, "character",
	'y.lab' , 'Y', 1, "character",
	'title.lab' , 'T', 1, "character",
	'lab.size' , 'l', 1, "integer",
	'axis.size' , 's', 1, "integer",
	'no.grid' , 'r', 0, "logical",
	'skip' , 'k', 1, "integer"
	), byrow=TRUE, ncol=4);
opt = getopt(spec);


# define usage function
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
1) Rscript singleArea.r --infile in_singleArea.data --outfile out_singleArea.png --x.col 1 --y.col 2 \\
	--x.lab \"x lab\" --y.lab \"y lab\" --skip 1 --no.grid
2) Rscript singleArea.r --infile in_singleArea.data --outfile out_singleArea.png --x.col 1 --y.col 2 \\
	--x.lab \"x lab\" --y.lab \"y lab\" --skip 1 --color \"#FF66EE\" --no.grid

Options: 
--help		-h 	NULL 		get this help
--infile 	-i 	character 	the input file [forced]
--outfile 	-o 	character 	the filename for output graph [forced]
--x.col 	-x 	integer 	the col for x value [forced]
--y.col 	-y 	integer 	the col for y value [forced]
--height 	-H 	integer 	the height of graph [optional, default: 3000]
--width 	-W 	integer 	the width of graph [optional, default: 4000]
--color 	-c 	character 	the color of distribution, RGB value( ","#FF00FF ) [optional, default: gray]
--x.lab 	-X 	character 	the lab for x [forced]
--y.lab 	-Y 	character 	the lab for y [forced]
--title.lab 	-T 	character 	the lab for title [optional, default: NULL]
--lab.size 	-l 	integer 	the font size of lab [optional, default: 14]
--axis.size 	-s 	integer 	the font size of text for axis [optional, default: 14]
--no.grid	-r 	NULL 		Do not drawing grid
--skip 		-k 	integer 	the number of line for skipping [optional, default: 0]
\n")
	q(status=1);
}



# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }


# check non-null args
if ( is.null(opt$infile) )	{ print_usage(spec) }
if ( is.null(opt$outfile) )	{ print_usage(spec) }
if ( is.null(opt$x.col) )	{ print_usage(spec) }
if ( is.null(opt$y.col) )	{ print_usage(spec) }
if ( is.null(opt$x.lab) )	{ print_usage(spec) }
if ( is.null(opt$y.lab) )	{ print_usage(spec) }


#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$skip ) )		{ opt$skip = 0 }
if ( is.null(opt$height ) )		{ opt$height = 3000 }
if ( is.null(opt$width ) )		{ opt$width = 4000 }
if ( is.null(opt$lab.size ) )		{ opt$lab.size = 14 }
if ( is.null(opt$axis.size ) )		{ opt$axis.size = 14 }
if ( is.null(opt$title.lab) )		{ opt$title.lab = NULL }
if ( is.null(opt$color) )		{ opt$color = brewer.pal(9,"Set1")[9] }




#-----------------------------------------------------------------
# reading data
#-----------------------------------------------------------------
# reading data
data <- read.table(opt$infile, skip=opt$skip)
# check dim
data.dim <- dim(data)
if ( is.null(data.dim) ){
	cat("Final Error: the format of infile is error, dim(data) is NULL\n")
	print_usage(spec)
}
# check col size
if ( data.dim[2] < max(opt$x.col, opt$y.col) ){
	cat("Final Error: max(x.col, y.col) > the col of infile\n")
	print_usage(spec)
}
# create df
df <- data.frame(x=data[,opt$x.col], y=data[,opt$y.col])


mycol<-brewer.pal(4, "Paired")

#-----------------------------------------------------------------
# plot
#-----------------------------------------------------------------
# mian plot
#p <- ggplot(df, aes(x=x, y=y),binwidth=0.1) + geom_area(aes(color=group, fill=group))
p <- ggplot(df, aes(x=x, y=y),binwidth=0.1) + geom_area(colour=mycol[2], fill=mycol[1],alpha=1)
#p <- ggplot(df, aes(x=x, y=y),binwidth=0.1) + geom_line(aes(color=group, fill=group))
#p <- p + facet_grid(group~.)
## color
#if( opt$is.color ){
#	p <- ggplot(df, aes(x=x, y=y, group = group, colour=group))
#}else{
#	p <- ggplot(df, aes(x=x, y=y, group = group))
#}
## point
#shape.value <- NULL
#if( opt$is.shape ){
#	shape.value <- c(1:length(levels(df$group)))
#}else{
#	shape.value <- rep(1,length(levels(df$group)))
#}
#if( opt$is.point ){
#	p <- p + geom_point(aes(shape=group)) + scale_shape_manual(values=shape.value)
#}
## line
#linetype.value <- NULL
#if( opt$is.linetype ){
#	linetype.value <- c(1:length(levels(df$group)))
#}else{
#	linetype.value <- rep(1,length(levels(df$group)))
#}
#if( opt$is.line ){
#	p <- p + geom_line(aes(linetype=group))  + scale_linetype_manual(values=linetype.value)
#}
#p <- qplot(x, y, data=df, geom="point", colour=group)
#p <- qplot(x, y, data=df, geom="line", colour=group, linetype=group)
#p <- qplot(x, y, data=df, geom="line", colour=group, linetype=group) + geom_point(aes(shape=group))
#p <- qplot(x, y, data=df, geom="line", colour=group, linetype=group) + geom_point(aes(shape=group)) + scale_shape_manual(values=c(1:15))
#p <- qplot(x, y, data=df, geom="line", colour=group, linetype=group) + geom_point(aes(shape=group)) + scale_shape_manual(values=c(1:15)) + scale_linetype_manual(values=c(1:15))
#p <- ggplot(df, aes(x=x, y=y, group = group, colour=group)) + geom_line(aes(linetype=group)) + geom_point(aes(shape=group)) + scale_shape_manual(values=c(rep(1,15))) + scale_linetype_manual( values=c(1:length(levels(df$group))) )
# faceting
#if( is.null(opt$facet.ncol) ) {
#	p <- p + facet_wrap(~group)
#} else {
#	p <- p + facet_wrap(~group, ncol=opt$facet.ncol)
#}

#-----------------------------------------------------------------
# theme
#-----------------------------------------------------------------
# lab
p <- p + xlab(opt$x.lab) + ylab(opt$y.lab) + labs(title=opt$title.lab)
# set lab and axis test size
# remove legend
p <- p + theme(legend.position = "none")
# grid and background


p<-p+theme_bw()+ theme(  
		panel.grid=element_blank(), 
		axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),panel.border=element_rect(colour = "black"))

#-----------------------------------------------------------------
# output plot
#-----------------------------------------------------------------
pdf(file=paste(opt$outfile,".pdf",sep=""), height=opt$height*2/1000, width=opt$width*2/1000)
print(p)
dev.off()
png(filename=opt$outfile, height=opt$height, width=opt$width, res=500, units="px")
print(p)
dev.off()







