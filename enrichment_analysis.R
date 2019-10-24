#!/usr/bin/env Rscript
library(optparse)

#######arguments
option_list <- list( 
    make_option(c("-i", "--input"),metavar="path", dest="input",help="Specify the path of metaboanalyst input table",default=NULL),
    make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files",default="./")
    )

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to relate the bacteria species and enviromental factors(numeric), and a heatmap is used to visualize the rank correlation coefficent"))

library("MetaboAnalystR")
library("getopt")
base_dir<-normalizePath(dirname(get_Rscript_filename()))
hsa_mfn<-paste(base_dir,"/sources/hsa_mfn.rds", sep = "")

if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}

input<-normalizePath(opt$input)

setwd(opt$out)

write("######################Perform data processing", stdout())
mSet<-InitDataObjects("pktable", "stat", FALSE)
## [1] "MetaboAnalyst R objects initialized ..."
mSet<-Read.TextData(mSet, input, "row", "disc")
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-FilterVariable(mSet, "iqr", "F", 25)
## [1] " Further feature filtering based on Interquantile Range"
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "AutoNorm", ratio=FALSE, ratioNum=20)

write("######################Perform OPLS-DA", stdout())
dir.create("OPLS-DA",recursive = T)
setwd("OPLS-DA")
# View the OPLS-DA plot
mSet<-OPLSR.Anal(mSet, reg=TRUE)
mSet<-PlotOPLS2DScore(mSet, "opls_score2d_", "pdf", 300, width=30, 1,2,0.95,0,0)
vip<-data.frame(mSet$analSet$oplsda$vipVn)
colnames(vip)<-"VIP"
write.csv(vip,"opls_VIP.csv")
knitr::include_graphics("opls_score2d_dpi300.pdf")

write("######################Reperform data processing", stdout())
setwd('../')
# Re-perform normalization, without auto-scaling
mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)
# Perform t-test
mSet<-Ttests.Anal(mSet, F, 0.25, FALSE, TRUE)
# Convert the output of the t-test to the proper input for mummichog analysis
mSet<-Convert2Mummichog(mSet)

find_mummi_input<-function(){
    for(file in list.files()){
        if(grepl("mummichog_input",file)){
            return(file)
        }
    }
}

# Read in the ranked peak list (don't forget to change the file name!)
SetPeakFormat("mpt")
mSet<-UpdateInstrumentParameters(mSet, "5", "negative");
mSet<-Read.PeakListData(mSet, find_mummi_input());
mSet<-SanityCheckMummichogData(mSet)

write("######################Perform original mummichog algorithm", stdout())
dir.create("Mummichog",recursive = T)
setwd("Mummichog")
file.copy(from=hsa_mfn, to="./")
# Perform original mummichog algorithm
mSet<-SetPeakEnrichMethod(mSet, "mum")
mSet<-SetMummichogPval(mSet, 0.25)
mSet<-PerformPSEA(mSet, "hsa_mfn", permNum = 1000)
mSet<-PlotPeaks2Paths(mSet, "peaks_to_paths_mummichog_", "pdf", 300, width=30)
file.remove("./hsa_mfn.rds")

write("######################Perform GSEA", stdout())
dir.create("../GSEA",recursive = T)
setwd("../GSEA")
file.copy(from=hsa_mfn, to="./")
# Perform GSEA
mSet<-SetPeakEnrichMethod(mSet, "gsea")
mSet<-PerformPSEA(mSet, "hsa_mfn", permNum = 1000)
mSet<-PlotPeaks2Paths(mSet, "peaks_to_paths_GSEA_", "pdf", 300, width=30)
knitr::include_graphics("peaks_to_paths_GSEA_dpi300.pdf")
file.remove("./hsa_mfn.rds")

write("######################Perform GSEA-mummichog", stdout())
dir.create("../GSEA-mummichog",recursive = T)
setwd("../GSEA-mummichog")
file.copy(from=hsa_mfn, to="./")
# Perform GSEA
mSet<-SetPeakEnrichMethod(mSet, "integ")
mSet<-SetMummichogPval(mSet, 0.25)
mSet<-PerformPSEA(mSet, "hsa_mfn")
mSet<-PlotIntegPaths(mSet, "integ_peaks_", "pdf", 300, width=30)
knitr::include_graphics("integ_peaks_dpi300.pdf")
file.remove("./hsa_mfn.rds")