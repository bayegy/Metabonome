setwd("/home/cheng/pipelines/testdir/testbono/")
df<-read.csv("oplsda/oplsda_features_importance.csv")
plot_volcano(data = df, out_img = "test.pdf",  axis=c(4,8), threshold = c(100,100), label=1)

OPLSDA("iHMP2_48_metaboanalyst_input.csv","oplsda")
tf<-groupCenter("iHMP2_48_metaboanalyst_input.csv")
enrichment_from_peak("iHMP2_48_metaboanalyst_input.csv", "hsa","new")

rownames(tf)


result<-mSet$analSet
summary(mSet$matches.res)
cpd<-mSet$matches.res
max(sapply(cpd, FUN=length))
mSet$mz2cpd_dict

help("PerformPSEA")

help(grepl)

mummichog.lib <- readRDS("hsa_mfn.rds")
cpd.treep <- mummichog.lib$cpd.tree[["positive"]];
cpd.treen <- mummichog.lib$cpd.tree[["negative"]]

cpd.treen[c(1900,1901,902)]

syn<-readRDS("syn_nms.rds")
syn<-syn$syns.list
setwd("/home/cheng/pipelines/testdir/testbono/")
set_db_location("/home/cheng/pipelines/Metabonome/database/")
pathway_analysis("human_cachexia.csv","pathway2")

require_db("hsa.rda")
data.frame(check.names = )
load("hsa.rda")
load("test_topo_database/hsa.rda")

network_map("human_cachexia.csv","test_network2")

rm(list = ls())

t<-c()
writeChar(t, mSet$dataSet$cmpd.mat)
help("writeChar")

remove.packages("MetaboAnalystR")
rm(list = ls())
load("/home/cheng/pipelines/Metabonome/database/hsa/hsa.rda")
df <-  readRDS("pathinteg.impTopo")
