#!/usr/bin/env Rscript
library(optparse)
library(purrr)
#######arguments
option_list <- list(
    make_option(c("-i", "--input"),metavar="path", dest="cpd",help="Specify the path of metabolites abundance table. Required",default=NULL),
    make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file. Required",default=NULL),
    # make_option(c("-q", "--qc-meta"),metavar="path",dest="qc", help="Specify the path of metadata file for QC.",default=NULL),
    make_option(c("-p", "--pca"),metavar="paths",dest="pca", help="comma seprated paths to QC completed compound tables.",default=NULL),
    make_option(c("-c", "--category"),metavar="string",dest="category", help="Category to compare. Required",default=NULL),
    make_option(c("-C", "--colors"),metavar="string",dest="colors", help="Comma seprated group colors.",default=NULL),
    make_option(c("-d", "--database"),metavar="path", dest="db",help="Path to the database directory.",
      # default="/home/cheng/Databases/Metabonome_database/database",
      default="/home/bayegy/Databases/metabolome"
    ),
    make_option(c("-s", "--species"),metavar="string", dest="species",help="The organism type of metabolites. Available now: hsa, mmu, rno. hsa means human; mmu means small mouse; rno means big rat.",default="hsa"),
    make_option(c("-t", "--type"),metavar="string", dest="type",help="The input data type, (target_name, untarget_name, untarget_keggid_name, untarget_kegg_id). default: untarget_keggid_name",default="untarget_keggid_name"),
    make_option(c("-f", "--flow"),metavar="string", dest="flow",help="Experiment work flow type, (lc_ms, gc_ms or gc_q_ms)",default="gc_q_ms"),
    make_option(c("-n", "--normalize"),metavar="logical", dest="normalize",help="normalize the abundance table, T or F. default T",default="T"),
    make_option(c("-r", "--report"),metavar="logical", dest="report",help="generate report, T or F. default T",default="T"),
    make_option(c("-o", "--out-dir"),metavar="directory",dest="out", help="Specify the directory of output files",default="./")
    )
opt <- parse_args(OptionParser(option_list=option_list,description = "The metabonome analysis pipeline"))

#####packages
library("purrr")
library("RColorBrewer")
library("statTarget")
library("getopt")
library("stringr")
library("VennDiagram")
# library("ggbiplot")

#####arguments
query_type = "kegg_id"
out_dir = opt$out # "/home/cheng/pipelines/testdir/cpd_pipline"
db_dir = opt$db # "/home/cheng/pipelines/Metabonome/database/"
cpd_file = opt$cpd # "/home/cheng/pipelines/testdir/cpd_pipline/compound_abundance.csv"
# metadata_file = opt$qc # "/home/cheng/pipelines/testdir/cpd_pipline/metadata.csv"
# run_qc <- !is.null(metadata_file)
map_file = opt$map # "/home/cheng/pipelines/testdir/cpd_pipline/mapping_file.txt"
categories = str_split(opt$category, ',')[[1]]
organism = opt$species
base_dir<-normalizePath(dirname(get_Rscript_filename()))
source(paste(base_dir, "metabonome_core.R", sep = "/"))
pcas <- opt$pca
input_type <- opt$type
exp_flow <- opt$flow
colors<-opt$colors
report<-as.logical(opt$report)
NORMALIZE_DATA <<- as.logical(opt$normalize)
WRITE_INPUT <<- TRUE
#TEST
# setwd("/media/cheng/disk1/Metabonome_Projects/cheli")
# query_type = "kegg_id"
# out_dir = "results" # "/home/cheng/pipelines/testdir/cpd_pipline"
# db_dir = "/home/cheng/Databases/Metabonome_database/database"
# cpd_file = "msms.csv" # "/home/cheng/pipelines/testdir/cpd_pipline/compound_abundance.csv"
# metadata_file = NULL # "/home/cheng/pipelines/testdir/cpd_pipline/metadata.csv"
# run_qc <- !is.null(metadata_file)
# map_file = "map.csv" # "/home/cheng/pipelines/testdir/cpd_pipline/mapping_file.txt"
# categories = "Category"
# organism = "rno"
# base_dir<-"/home/cheng/pipelines/Metabonome/"
# source(paste(base_dir, "metabonome_core.R", sep = "/"))
# pcas <- "data/metabolome_neg.csv,data/metabolome_pos.csv"
# input_type <- "untarget"

#####arguments check
if(!dir.exists(out_dir)){
  dir.create(out_dir)
}
out_dir <- normalizePath(out_dir)
db_dir <- normalizePath(db_dir)
cpd_file <- normalizePath(cpd_file)
map_file <- normalizePath(map_file)
root_dir <- paste(out_dir, "/ResultsMetabonome", sep = "")
pca_files <- choose(is.null(pcas), pcas, normalizePath(str_split(pcas, ',')[[1]]))


for(file in c(pca_files, map_file, cpd_file)){
  system(sprintf("%s/strip_input_file.py %s ,", base_dir, file))
}

sprintf("Results is saved at %s", root_dir)
######dir preparation
moveto_dir <- function(dir="/"){
  dir <- paste(root_dir, dir, sep = "/")
  if(!dir.exists(dir)){
    dir.create(dir, recursive = TRUE)
  }
  setwd(dir = dir)
}

df_map<-function(df, func){
  out<-c()
  for(i in 1:ncol(df)){
    out[i]<-func(df[, i])
  }
  return(out)
}

set_db_location(db_dir)

write("#####################Read Table", stdout())
# cpd_table <- read.table(cpd_file, sep = "\t", header = T, row.names = 1, stringsAsFactors = FALSE)
mapping_table <- read.csv(map_file, header = T, na.strings="", row.names = 1, stringsAsFactors = FALSE)
sample_id_map <- c()
if('Description'%in%colnames(mapping_table)){
  actual_samid <- mapping_table['Description'][,1]
}else{
  actual_samid <- rownames(mapping_table)
}

sample_id_map[rownames(mapping_table)] <- actual_samid
rownames(mapping_table) <- actual_samid

############ write out input #########
moveto_dir()
cpd_table_out <- read.csv(cpd_file, check.names = FALSE, header = T, stringsAsFactors = FALSE)
write.table(
  data.frame(compounds = cpd_table_out[, 1], cpd_table_out[df_map(cpd_table_out, is.numeric)]),
  "compounds.tsv", sep = "\t", row.names=FALSE
)
write.table(
  data.frame(SampleID=rownames(mapping_table), mapping_table[!colnames(mapping_table)%in%c("order")]),
  "mapping_file.tsv", sep = "\t", row.names=FALSE
)


write("#####################Color Assignment", stdout())
if(is.null(colors)){
  pallet <- unique(c(c("#20B2AA", "#DC143C", "#FF9933"), rev(brewer.pal(12,"Paired")),brewer.pal(8,"Set2")[-c(7,8)],brewer.pal(8,"Dark2"),brewer.pal(12,"Set3"),brewer.pal(8,"Accent"),brewer.pal(11,"Spectral")))
}else{
  pallet <- str_split(colors, ",")[[1]]
}

group_df <- mapping_table[categories]
group_df$id <- rownames(group_df)
group_df <- na.omit(melt(group_df, id.vars=c("id")))
# print(group_df)
all_groups <- unique(group_df[,ncol(group_df)])
all_groups <- sort(all_groups)
group_colors <- pallet[1:length(all_groups)]
names(group_colors) <- all_groups
pre_ordered_samples <- rownames(mapping_table)[order(mapping_table["order"][, 1])]

# print(group_colors)

yield_combine <- function(vec){
  vec <- sort(vec)
  out_list <- list()
  index <- 1
  for(i in 2:length(vec)){
    sub_com <- combn(vec, i)
    for(j in 1:ncol(sub_com)){
      out_list[[index]] <- sub_com[, j]
      index = index +1
    }
  }
  return(out_list)
}

select_num<-function(df, rownamecol=NA, save_col=NA){
  if(!is.na(rownamecol)){
    rownames(df) <- df[rownamecol][, 1]
  }
  col_is_num <- c()
  for(c in 1:ncol(df)){
    col_is_num[c] <- is.numeric(df[, c])
  }
  df_num <- df[, col_is_num]
  if(!is.na(save_col)){
    df_num <- data.frame(df[save_col], df_num, check.rows = TRUE, check.names = FALSE)
  }
  return(df_num)
}

select_cpd_df <- function(cpd_label="name", fix_columns = c(1:9)){
  label <- cpd_df[cpd_label][,1]
  selected <- !(is.na(label)|label=="")
  out_df <- cpd_df[selected, ]
  label <- label[selected]
  out_df <- out_df[,-fix_columns]
  out_df <- apply(out_df, 2, function(x){return(base::tapply(x, INDEX = label, sum))})
  return(t(out_df))
}

data_factory <- function(category="Category", cpd_label = "name", type = "pair", cmpType=1, fix_columns = c(1:9)){
  map_df <- mapping_table[category]
  map_df <- na.omit(map_df)
  #print(map_df)
  #print(select_cpd_df(cpd_label = cpd_label, fix_columns = fix_columns))
  sub_cpd_df <- select_cpd_df(cpd_label = cpd_label, fix_columns = fix_columns)
  map_samples <- rownames(map_df)
  cpd_samples <- rownames(sub_cpd_df)

  intersect_samples <- intersect(cpd_samples, map_samples)

  print(sprintf("FOR %s:", category))
  print("Following samples found in mapping file, but not in compound abundance table:")
  print(map_samples[!map_samples%in%cpd_samples])
  print("Following samples found in abundance table, but not in compound mapping file:")
  print(cpd_samples[!cpd_samples%in%map_samples])

  map_df <- map_df[map_samples%in%intersect_samples, ]
  map_df <- data.frame(map_df, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(map_df) <- category
  rownames(map_df) <- map_samples[map_samples%in%intersect_samples]

  sub_cpd_df <- sub_cpd_df[cpd_samples%in%intersect_samples, ]

  # print(map_df)
  # print(sub_cpd_df)

  base_df <- rownames_join(map_df, sub_cpd_df)
  # print(base_df)
  base_df <- base_df[order(base_df[,1]), ]
  base_df <- data.frame(SampleID=rownames(base_df), base_df, check.names = FALSE, stringsAsFactors = FALSE)
  out_list <- list()
  if(type=="pair"){
    groups <- unique(map_df[,1])
    paired_groups <- yield_combine(groups)
    nums <- length(paired_groups)
    cmpType <- choose(length(cmpType)==1, rep(cmpType, nums), cmpType)
    for(i in 1:nums){
      ele_groups <- paired_groups[[i]]
      ele_name <- reduce(ele_groups, paste, sep=".vs.")
      ele_is_paired <- length(ele_groups)==2
      ele_cmpType <- cmpType[i]
      ele_data <- base_df[base_df[,2]%in%ele_groups, ]
      # ��ȥ����ֵ����NA�Ļ�����
      ele_data <- ele_data[, c(T,T,colSums(!is.na(ele_data[,-c(1,2)]))>0)]
      ele_data <- ele_data[, c(T,T,colSums(!ele_data[,-c(1,2)]==0)>0)]
      ele_colors <- group_colors[ele_groups]
      # print(ele_groups)
      # print(ele_colors)
      out_list[[i]] <- list(data=ele_data, name=ele_name, groups=ele_groups, colors=ele_colors, cmpType=ele_cmpType, paired=ele_is_paired)
    }
    return(out_list)
  }else{
    if(type=="single"){
      # print("Here !!!!!!!!!!!!!!!!!!!!!!!!!!!")
      groups <- unique(base_df[, 2])
      for (i in 1:length(groups)) {
        ele_data = base_df[base_df[, 2]==groups[i], ]
        ele_data <- ele_data[, c(T,T,colSums(!is.na(ele_data[,-c(1,2)]))>0)]
        out_list[[i]] <- list(name=groups[i], data=ele_data, colors=group_colors[groups[i]])
      }
    }else{
     out_list <- base_df
    }
  }
  return(out_list)
}


map_table_header <- function(cpd_table){
  sample_index <- colnames(cpd_table)%in%names(sample_id_map)
  colnames(cpd_table)[sample_index] <- sample_id_map[colnames(cpd_table)[sample_index]]
  return(cpd_table)
}



for(category in categories){
  # category <- "Category"
  write("#####################Data Prepare", stdout())
  moveto_dir(paste(category, "/00-DataProcess/", sep = "/"))

  # map
  cpd_table <- read.csv(cpd_file, check.names = FALSE, header = T, stringsAsFactors = FALSE)
  cpd_table <- map_table_header(cpd_table=cpd_table)
  if(input_type=="untarget_keggid_name"){
    # if(run_qc){
    #   dir.create("./01-QA_QC")
    #   metadata <- read.csv(metadata_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    #   map_category <- mapping_table[category]
    #   metadata <- rownames_join(metadata, map_category)
    #   colnames(metadata)[3] <- "class"
    #   write.csv(data.frame(sample=rownames(metadata), metadata, check.names = FALSE), file = "metadata.csv", row.names = FALSE)
    #   shiftCor("metadata.csv", cpd_file, Frule = 0.8, MLmethod = "QCRFSC", QCspan = 0,imputeM = "KNN", plot = TRUE)
    #   cpd_table <- read.csv("statTarget/shiftCor/After_shiftCor/shift_sample_cor.csv", check.names = FALSE, header = T, row.names = 1, stringsAsFactors = FALSE)
    #   system("mv statTarget/shiftCor/* ./01-QA_QC")
    #   # system("mv statTarget/shiftCor/RSDresult/* 01-QA_QC")
    #   system("rm -r statTarget")
    #   cpd_table <- t(cpd_table[, -1])
    #   cpd_names <-  rownames(cpd_table)
    #   ref_df <- cross_reference(cpd_names, type="name")
    #   # rownames(ref_df)<-ref_df[,1]
    #   cpd_df <- data.frame(name=cpd_names, ref_df, cpd_table, check.names=FALSE)
    # }else{
    cpds <- cpd_table[, 1] # First column is keggid
    ref_df <- cross_reference(cpds, type="kegg_id")
    # Second column is metabolite name
    cpd_df <- data.frame(cpd_table[, 2], ref_df, cpd_table[, -c(1, 2)], check.names=FALSE)
    colnames(cpd_df)[1] <- "name"
    # cpd_df["kegg_id"] <- cpds
    # }
    cpd_df<-cpd_df[, -c(2,3)]
    write.csv(cpd_df, "compound_cross_reference_after_QC.csv", row.names = F)
    fix_columns <- c(1:9)
  }else if(input_type=="target_name"){
    # cpd_df <- read.csv(cpd_file, check.names = FALSE, header = T, stringsAsFactors = FALSE)
    cpd_df <- cpd_table
    # colnames(cpd_df)[1:2]<-c("name", "class")
    # fix_columns <- c(1:2)
    colnames(cpd_df)[1]<-c("name")
    fix_columns <- c(1)
  }else{
    # untarget_name
    id_type <- str_split(input_type, "_", n=2)[[1]][2]
    cpd_names <- cpd_table[, 1]
    ref_df <- cross_reference(cpd_names, type=id_type)
    cpd_df <- data.frame(choose(id_type=="name",cpd_names, ref_df[1]), ref_df, cpd_table[, -c(1)], check.names=FALSE)
    cpd_df <- cpd_df[!is.na(cpd_df[, 1]), ]
    colnames(cpd_df)[1] <- "name"
    cpd_df<-cpd_df[, -c(2,3)]
    write.csv(cpd_df, "compound_cross_reference_after_QC.csv", row.names = F)
    fix_columns <- c(1:9)
  }

  rownames(cpd_df) <- cpd_df[, 1]
  write("#####################QC plot", stdout())
  moveto_dir(paste(category, "/00-DataProcess/01-QA_QC", sep = "/"))
  if(length(pca_files)>0){
    for(pca_file in  pca_files){
      print(sprintf("Plot PCA for %s",  pca_file))
      pca_df <- read.csv(pca_file, header=TRUE,  row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
      pca_df <- map_table_header(pca_df)
      pca_df <- t(pca_df)
      sam_gp <- rownames(pca_df)
      is_qc <- grepl("QC", sam_gp)
      sam_gp[is_qc] <- "QC"
      sam_gp[!is_qc] <- "Sample"
      pca.data <- data.frame(Sample=rownames(pca_df), Category=sam_gp, pca_df)
      rownames(pca.data)<-rownames(pca_df)
      # pdf("")
      file=paste(str_replace(pca_file, "\\.[^\\.]+$", ""), "_PCA_QC.pdf", sep="")
      file_withid=paste(str_replace(pca_file, "\\.[^\\.]+$", ""), "_PCA_QC_with_ids.pdf", sep="")
      file=str_replace(file, "^.+/", "")
      file_withid=str_replace(file_withid, "^.+/", "")
      # print(pca.data)
      elipse_pca(table = pca.data, file_biplot = file, file_pc1 = NA, label = F)
      elipse_pca(table = pca.data, file_biplot = file_withid, file_pc1 = NA, label = T)
    }
  }
  write("#####################QA plot", stdout())
  data <- data_factory(category = category, type = "all", fix_columns=fix_columns)
  # print(data)
  ordered_samples <- pre_ordered_samples[pre_ordered_samples%in%data[,1]]
  elipse_pca(table = data, file_pc1 = "AllSample_PC1.pdf", file_biplot = NA, sam_order = ordered_samples)
  datas <- data_factory(category = category, type = "single", fix_columns=fix_columns)
  for (data in datas) {
      # print(data$data)
      elipse_pca(table = data$data, colors = data$colors, file_biplot = paste(data$name, "_PCA_QA.pdf", sep = ""), file_pc1 = NA, label = T)
      elipse_pca(table = data$data, colors = data$colors, file_biplot = paste(data$name, "_PCA_QA_without_ids.pdf", sep = ""), file_pc1 = NA, label = F)
  }

  write("#####################Data Normlization", stdout())
  moveto_dir(paste(category, "/00-DataProcess/02-Normlization", sep = "/"))
  df <- data_factory(category=category, type = "all", fix_columns=fix_columns)
  mSet <- prepare(df)
  mSet <- PlotNormSummary(mSetObj = mSet, "compound_wise_normlization.pdf", format = "pdf", dpi = 300, width = 30)
  mSet <- PlotSampleNormSummary(mSetObj = mSet, "sample_wise_normlization.pdf", format = "pdf", dpi = 300, width = 30)
  write.csv(mSet$dataSet$row.norm, "sample_wise_normalized.csv")
  write.csv(mSet$dataSet$norm, "compound_wise_normalized.csv")

  if("class" %in% colnames(cpd_df)){
    datas <- data_factory(category = category, cpd_label = "class", fix_columns=fix_columns)
    for(data in datas){
      write("#####################ConcentrationSummary", stdout())
      moveto_dir(paste(category, "01-ConcentrationSummary/1-Barplot", data$name, sep = "/"))
      abundance_barplot(data$data, colors = data$colors, by_group_mean = FALSE, prefix = "Compounds_with_Biological_Roles_")   
      abundance_barplot(data$data, colors = data$colors, by_group_mean = TRUE, prefix = "Compounds_with_Biological_Roles_group_mean_")
    }
  }
  
  datas <- data_factory(category = category, fix_columns=fix_columns)
  report_group <- datas[[1]]$name
  sig_features_set <- list()
  for(data in datas){
    write("#####################ConcentrationSummary", stdout())
    moveto_dir(paste(category, "01-ConcentrationSummary/1-Barplot", data$name, sep = "/"))
    abundance_barplot(data$data, colors = data$colors, by_group_mean = FALSE)   
    abundance_barplot(data$data, colors = data$colors, by_group_mean = TRUE, prefix = "group_mean_")   

    moveto_dir(paste(category, "01-ConcentrationSummary/2-Heatmap", data$name, sep = "/"))
    # print(data$data)
    abundance_heatmap(data$data, colors = data$colors, by_group_mean = FALSE, cluster = FALSE)   
    abundance_heatmap(data$data, colors = data$colors, by_group_mean = FALSE, cluster = TRUE, prefix = "clustered_")   
    abundance_heatmap(data$data, colors = data$colors, by_group_mean = TRUE,  cluster = FALSE, prefix = "group_mean_")  
    
    write("#####################PCA Analysis", stdout())
    moveto_dir(paste(category, "02-PCA", sep = "/"))
    PCA(data$data, out_dir = data$name, colors = data$colors)
    
    write("#####################PLSDA Analysis", stdout())
    moveto_dir(paste(category, "03-PLSDA", sep = "/"))
    PLSDA(data$data, out_dir = data$name, colors = data$colors)
    
    write("#####################OPLSDA Analysis", stdout())
    moveto_dir(paste(category, "04-OPLSDA", sep = "/"))
    OPLSDA(data$data, out_dir = data$name, colors = data$colors)   

    write("#####################UnivariateAnalysis Analysis", stdout())
    moveto_dir(paste(category, "05-UnivariateAnalysis", data$name, sep = "/"))
    if(data$paired){
      print("Paired data, plot valcano")
      # Volcano(data$data, colors = data$colors, cmpType = data$cmpType)
      Volcano(data$data, cmpType = data$cmpType) # colors is not needed for volcano
      sig_features <- sig_boxplot(table = data$data, colors = data$colors)
      # print(sig_features)
      if(length(sig_features)>0){
        sig_features_set[[data$name]] <- sig_features
      }
    }else{
      sig_features <- sig_boxplot(table = data$data, colors = data$colors)
    }
    
    
    write("#####################RandomForest Analysis", stdout())
    moveto_dir(paste(category, "06-RandomForest", sep = "/"))
    RF(data$data, out_dir = data$name)

    write("#####################SupportVectorMachine Analysis ", stdout())
    moveto_dir(paste(category, "07-SupportVectorMachine", sep = "/"))
    if(data$paired){
      SVM(data$data, out_dir = data$name)
    }
  }
  
  
  write("#####################Venn Plot ", stdout())
  # if(length(sig_features_set)>1){
  #   moveto_dir(paste(category, "05-UnivariateAnalysis", sep = "/"))
  #   vs_combs <- yield_combine(names(sig_features_set))
  #   for (comb in vs_combs) {
  #     if(length(comb)<6){
  #       file_base <- reduce(comb, paste, sep="_and_")
  #       venn.diagram(sig_features_set[comb],
  #                    filename =paste(file_base, "_venn.tiff", sep = ""),
  #                    imagetype="tiff",alpha= 0.50,lwd =1.2,cat.cex=1.4,
  #                    fill=rainbow(length(comb)),
  #                    margin=0.15)
  #     }
  #   }
  #   system("rm *.log")
  # }
  
  if(!input_type=="target_name"){
    datas <- data_factory(category = category, cpd_label = "kegg_id", fix_columns=fix_columns)
    for(data in datas){
      write("#####################PathwayTopoEnrichment Analysis ", stdout())
      moveto_dir(paste(category, "09-PathwayTopoEnrichment", data$name, sep = "/"))
      pathway_analysis(data$data, out_dir = "./", threshp = 0.05, organism = organism)
      tryCatch({
        if(file.exists("selected_compounds.csv")){
          df <- read.csv("selected_compounds.csv", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
          ora.df <- read.csv("pathway_enrichment_and_topo_results.csv", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
          rownames(ora.df) <- str_replace(rownames(ora.df), '^\\D*', 'map')
          write.csv(ora.df, "pathway_enrichment_and_topo_results.csv")
          gc <- groupCenter(prepare(data$data))
          df <- rownames_join(df, gc)
          write.csv(df, "selected_compounds.csv")
          sig_pathways <- ora.df["Description"][ora.df["Raw p"]<=0.05]
          sig_pathways_id <- rownames(ora.df)[ora.df["Raw p"]<=0.05]
          ele_colors <- data$colors
          ele_colors <- paste(reduce(names(ele_colors), paste, sep=","), reduce(ele_colors, paste, sep=","), sep = ";")
          for(i in 1:length(sig_pathways_id)){
            mapid <- paste("map", str_extract(sig_pathways_id[i], "\\d+"), sep = "")
            if(!(mapid=="map"|mapid=="mapNA")){
              print(mapid)
              sig_path <- str_replace_all(sig_pathways[i], " |,", "_")
              cmd <- sprintf("CColorMap.py -i selected_compounds.csv -m %s -c '%s' -n 'Which-Max' -p '%s'", mapid, ele_colors, sig_path)
              print(cmd)
              system(cmd)
            }
          }
        }
      }, error=function(e){print(e)})
    }
  }
  show_label <- !exp_flow=="lc_ms"
  datas <- data_factory(category = category, type = "single", fix_columns=fix_columns)
  for (data in datas) {
    # data <- datas[[4]]
    write("#####################Correlation Analysis", stdout())
    moveto_dir(paste(category, "08-CorrelationAnalysis", data$name, sep = "/"))
    group_cor_heatmap(data$data[, -c(1,2)], show_names = show_label)
    # group_cor_heatmap(data$data[, -c(1,2)], cluster = TRUE, prefix = "clustered_", show_names = show_label)
  }
  data <- data_factory(category = category, type = "all", fix_columns=fix_columns)
  data <- data[, -c(1, 2)]
  moveto_dir(paste(category, "08-CorrelationAnalysis", sep = "/"))
  group_cor_heatmap(data, prefix = "all_group_", show_names = show_label)
  # group_cor_heatmap(data, cluster = TRUE, prefix = "all_group_clustered_", show_names = show_label)
  moveto_dir(category)
  if(report){
    system(sprintf("cp -rp %s/FiguresTablesForReport ./", base_dir))
    system("mv FiguresTablesForReport/结题报告.html ./")
    system(sprintf("bash %s/prepare_report_images.sh %s", base_dir, report_group))
    system(sprintf("render_report.py 结题报告.html %s %s", input_type, exp_flow))
  }
}

print("metabonome_pipline.R program is finished.")
print(Sys.time())
