library(optparse)

#######arguments
option_list <- list( 
    make_option(c("-i", "--input"),metavar="path", dest="otu",help="Specify the path of collapsed bacteria table",default=NULL),
    make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file",default=NULL),
    make_option(c("-e", "--exclude"),metavar="string",dest="ex", help="Specify the numeric variables excluded from plot and seprated by commas in mapping file",default="none"),
    make_option(c("-n", "--number"),metavar="int", dest="num",help="The number of most related species you want to plot, default is 20",default=30),
    make_option(c("-r", "--min-cor"),metavar="float", dest="minr",help="Min correlation coefficent to label significance, default is 0",default=0),
    make_option(c("-p", "--prefix"),metavar="str", dest="prefix",help="The prefix of output files, default if null",default=""),
    make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files",default="./")
    )

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to relate the bacteria species and enviromental factors(numeric), and a heatmap is used to visualize the rank correlation coefficent"))

library(pheatmap)
library(psych)
library(stringr)


if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}

opt$out<-paste(opt$out,"/",opt$prefix,sep="")
########prepare the cor data