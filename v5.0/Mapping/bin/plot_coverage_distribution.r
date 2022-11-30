#!/share/nas2/genome/biosoft/R/current/bin/Rscript

#####################################################################
# Copyright 2014, BMK
#
# Author: macx <macx@biomarker.com.cn>
#
# Function: draw genomewide cytosine coverage distribution map
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
	'window', 'w', 2, "integer", "window size, forced",
	'nchro', 'n', 2, "integer", "Choose the top n chromosome to draw, default [0]",
	'x.title', 'x', 2, "character", "x title, default [Chromosome position]",
	'y.title', 'y', 2, "character", "y title, default [Median reads density(log2)]",
	'title', 't', 2, "character", "graph title, default [Genomewide distribution of base coverage depth]"
	
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

#+--------------------
# check options
#+--------------------
if ( !is.null(opt$help) | is.null(opt$input) | is.null(opt$output) | is.null(opt$window)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

#+--------------------
# some default options
#+--------------------

if ( is.null(opt$nchro) ) opt$nchro <- 0
if ( is.null(opt$x.title) ) opt$x.title <- "Chromosome position"
if ( is.null(opt$y.title) ) opt$y.title <- "Median reads density(log2)"
if ( is.null(opt$title) ) opt$title <- "Genomewide distribution of base coverage depth"

#+--------------------
# Main
#+--------------------

## load ggplot2 
library(ggplot2)
library(grid)

## get data
reads <- read.table(opt$input, header=FALSE, sep="\t", comment.char = "#")
colnames(reads) <- c("chromsome", "position", "coverage")

if (ncol(reads) != 3) stop("Input file format error: must be <chromsome> <position> <coverage>")

## subset data frame if chromsome in scaffold level and specified the top_n option
if (opt$nchro > 0) {
	chro_max <- aggregate(as.numeric(reads$position), by=list(reads$chromsome), FUN=max)
	chro_kept <- as.character(chro_max[order(chro_max$x, decreasing=TRUE),][1:opt$nchro,1])
	reads <- reads[reads$chromsome %in% chro_kept,]
}


## plot 
p <- ggplot(reads, aes(position*as.numeric(opt$window)/1e6, coverage), binwidth = 0.1) + geom_area(aes(color=chromsome, fill=chromsome), position="jitter") ## create plot
p <- p + facet_grid(chromsome~., space="fixed") ## facet by chromsome
p <- p + theme(legend.position="none")
#p <- p + theme(legend.background=element_blank(), legend.key=element_blank(), legend.title=element_blank(), legend.text=element_blank()) ## set legend 

x_max <- round(max(reads$position * as.numeric(opt$window)/1e6), digits=0)
x_step <- round((round(x_max / 5) / 10+0.5)) * 10
x_breaks <- seq(0, x_max, x_step)

#p <- p + scale_y_continuous(limits=c(-0.5, 0.5), breaks=seq(-0.5, 0.5, 1)) + scale_x_continuous(breaks=x_breaks, labels=paste(x_breaks, "M", sep=""))  ## set x and y lables
p <- p + scale_x_continuous(breaks=x_breaks, labels=paste(x_breaks, "M", sep=""))  ## set x and y lables
p <- p + theme(strip.background=element_rect(fill="white", colour="white"), strip.text.y=element_text(angle=360, vjust=0.1), 
               panel.background=element_rect(fill="white", size=10), panel.grid=element_blank(), 
               axis.line=element_line(colour="grey")) ## set other theme opts 
p <- p + labs(title=opt$title, x=opt$x.title, y=opt$y.title) ## set labs


## output to png file 
png(filename=opt$output, height = 3000, width = 3000, res = 300, units = "px")
print(p)
dev.off()
