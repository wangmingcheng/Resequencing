#!/share/nas2/genome/biosoft/R/current/bin/Rscript

#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
# load library
library('getopt');


#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
	'help' , 'h', 0, "logical",
	'infile' , 'i', 1, "character",
	'outfile' , 'o', 1, "character",
	'x.col' , 'v', 1, "integer",
	'y1.col' , 'a', 1, "integer",
	'y2.col' , 'b', 1, "integer",
	'height' , 'H', 1, "integer",
	'width' , 'W', 1, "integer",
	'x.lab' , 'V', 1, "character",
	'y1.lab' , 'A', 1, "character",
	'y2.lab' , 'B', 1, "character",
	'title.lab' , 'T', 1, "character",
	'legend.xpos' , 'x', 1, "double",
	'legend.ypos' , 'y', 1, "double",
	'skip' , 'k', 1, "integer"
	), byrow=TRUE, ncol=4);
opt = getopt(spec);


# define usage function
print_usage <- function(spec=NULL){
	cat("\n")
	cat(getopt(spec, usage=TRUE));
	cat("
Usage example: 
1) Rscript dual_axis.r --infile in.data --outfile out.png --height 3000 --width 5000 \\
	--x.col 1 --y1.col 2 --y2.col 3 --x.lab \"x lab\" --y1.lab \"y1 lab\" \\
	--y2.lab \"y2 lab\" --title.lab \"title: hello world\"
2) Rscript dual_axis.r --infile in.data --outfile out.png \\
	--x.col 1 --y1.col 2 --y2.col 3 --x.lab \"x lab\" --y1.lab \"y1 lab\" \\
	--y2.lab \"y2 lab\" --title.lab \"title: hello world\" \\
	--legend.xpos 0.7 --legend.ypos 0.85

Options: 
--help		-h 	NULL 		get this help
--infile 	-i 	character 	the input file [forced]
--outfile 	-o 	character 	the filename for output graph [forced]
--x.col 	-v 	integer 	the col for group factor [forced]
--y1.col 	-a 	integer 	the col for x value [forced]
--y2.col 	-b 	integer 	the col for y value [forced]
--height 	-H 	integer 	the height of graph [optional, default: 3000]
--width 	-W 	integer 	the width of graph [optional, default: 4000]
--x.lab 	-V 	character 	the lab for group factor [forced]
--y1.lab 	-A 	character 	the lab for x [forced]
--y2.lab 	-B 	character 	the lab for y [forced]
--title.lab 	-T 	character 	the lab for title [optional, default: NULL]
--legend.xpos 	-x 	double 		the x relative position for legend, (0.0,1.0) [optional, default: 0.8]
--legend.ypos 	-y 	double 		the y relative position for legend, (0.0,1.0) [optional, default: 0.85]
--skip 		-k 	integer 	the number of line for skipping [optional, default: 1]
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
if ( is.null(opt$y1.col) )	{ print_usage(spec) }
if ( is.null(opt$y2.col) )	{ print_usage(spec) }
if ( is.null(opt$x.lab) )	{ print_usage(spec) }
if ( is.null(opt$y1.lab) )	{ print_usage(spec) }
if ( is.null(opt$y2.lab) )	{ print_usage(spec) }


#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$skip ) )		{ opt$skip = 0 }
if ( is.null(opt$legend.xpos ) )	{ opt$legend.xpos = 0.7 }
if ( is.null(opt$legend.ypos ) )	{ opt$legend.ypos = 0.85 }
if ( is.null(opt$title.lab) )		{ opt$title.lab = NULL }



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
if ( data.dim[2] < max(opt$x.col, opt$y1.col, opt$y2.col) ){
	cat("Final Error: max(x.col, y1.col, y2.col) > the col of infile\n")
	print_usage(spec)
}


library(RColorBrewer)
mycol<-brewer.pal(4, "Set1")

#-----------------------------------------------------------------
# define function
#-----------------------------------------------------------------
plot_dual_axis <- function(x=NULL, y1=NULL, y2=NULL, 
par.mar=c(5, 5.5, 4, 6) + 0.1,
line.col=mycol[1:2], 
line.lwd=c(2.5,2.5),
line.pch=c(15,16),
axis.lwd=2,
axis.cex=1.3,
text.cex=1.3,
text.line=c(2.5,3.5,3.5), # x, y1, y2
legend.cex=1,
legend.pos=c(0.05,0.98),
legend.box.lwd=1,
x.lab="x",
y1.lab="y1",
y2.lab="y2",
main="main",
xpd=TRUE,
limit=rep(0.05, 6), # x.min, x.max, y1.min, y1.max, y2.min, y2.max
... ){
	# get x limit
	x.min <- min(x)
	x.max <- max(x)
	x.min.limit <- x.min - limit[1] * (x.max-x.min)
	x.max.limit <- x.max + limit[2] * (x.max-x.min)

	## add extra space to right margin of plot within frame
	par(mar=par.mar)

	## Plot first set of data and draw its axis
	y1.min <- min(y1)
	y1.max <- max(y1)
	y1.min.limit <- y1.min - limit[3] * (y1.max-y1.min)
	y1.max.limit <- y1.max + limit[4] * (y1.max-y1.min)
	plot(x, y1, pch='', lwd=line.lwd[1], 
		axes=FALSE, xlim=c(x.min.limit,x.max.limit),
		ylim=c(y1.min.limit,y1.max.limit), 
		xlab="", ylab="", type="o",col=line.col[1], 
		main=main)
	axis(2, ylim=c(y1.min.limit,y1.max.limit),
		col="black", col.axis=mycol[1], las=1,col.ticks= mycol[1],
		lwd=axis.lwd, cex.axis=axis.cex) ## las=1 makes horizontal labels
	mtext(y1.lab, side=2, line=text.line[2], col=mycol[1], cex=text.cex)
	box(lwd=axis.lwd)

	## Allow a second plot on the same graph
	par(new=TRUE)
 
	## Plot the second plot and put axis scale on right
	y2.min <- min(y2)
	y2.max <- max(y2)
	y2.min.limit <- y2.min - limit[5] * (y2.max-y2.min)
	y2.max.limit <- y2.max + limit[6] * (y2.max-y2.min)
	plot(x, 1-y2, lwd=line.lwd[2],pch='',
		axes=FALSE, xlim=c(x.min.limit,x.max.limit),
		ylim=c(0,1.0), 
		xlab="", ylab="", type="o",col=line.col[2])
	## a little farther out (line=4) to make room for labels
	axis(4, ylim=c(y2.min.limit,y2.max.limit),
		col="black", col.axis=mycol[2], las=1,col.ticks= mycol[2],
		lwd=axis.lwd, cex.axis=axis.cex) ## las=1 makes horizontal labels
	mtext(y2.lab, side=4, line=text.line[3], col=mycol[2], cex=text.cex)

	## Draw the time axis
	axis(1, xlim=c(x.min.limit,x.max.limit), 
		col="black", col.axis="black",
		lwd=axis.lwd, cex.axis=axis.cex)
	mtext(x.lab, side=1, line=text.line[1], col="black", cex=text.cex)
	box(lwd=axis.lwd)

	## Add Legend
	x.pos <- x.min.limit + legend.pos[1] * (x.max.limit-x.min.limit)
	y.pos <- y2.min.limit + legend.pos[2] * (y2.max.limit-y2.min.limit)
	legend("right", legend=c(y1.lab,y2.lab),
		text.col=line.col, lty=c(1,1), 
		col=line.col, lwd=line.lwd[1], cex=legend.cex, bty="n",ncol=1 )

}




#-----------------------------------------------------------------
# reading data
#-----------------------------------------------------------------
# reading data
data <- read.table(opt$infile)



#-----------------------------------------------------------------
# plot
#-----------------------------------------------------------------
png(filename=opt$outfile, height = 2500, width = 3000, res = 300, units = "px")

plot_dual_axis(
as.numeric(x=data[,opt$x.col]), y1=as.numeric(data[,opt$y1.col]), y2=as.numeric(data[,opt$y2.col]),
x.lab=opt$x.lab,
y1.lab=opt$y1.lab,
y2.lab=opt$y2.lab,
main=opt$title.lab,
legend.pos=c(opt$legend.xpos,opt$legend.ypos),
text.line=c(2.5,3.5,3.5),
par.mar=c(4, 5, 3, 6) + 0.1,
text.cex=1.5,
legend.cex=1,
legend.box.lwd=1)

dev.off()









