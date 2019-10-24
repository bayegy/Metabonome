#!/usr/bin/env Rscript
library("purrr")
library("RColorBrewer")
library("statTarget")
library("getopt")
library("stringr")
#####arguments
out_dir = "/home/cheng/pipelines/testdir/cpd_pipline"
db_dir = "/home/cheng/pipelines/Metabonome/database/"
cpd_file = "/home/cheng/pipelines/testdir/cpd_pipline/compound_abundance.csv"
metadata_file = "/home/cheng/pipelines/testdir/cpd_pipline/metadata.csv"
map_file = "/home/cheng/pipelines/testdir/cpd_pipline/mapping_file.txt"
categories = c("Category")
organism = "hsa"
base_dir<-normalizePath(dirname(get_Rscript_filename()))
source(paste(base_dir, "metabonome_core.R", sep = "/"))
#####arguments check
out_dir <- normalizePath(out_dir)
db_dir <- normalizePath(db_dir)
cpd_file <- normalizePath(cpd_file)
map_file <- normalizePath(map_file)
root_dir <- paste(out_dir, "/ResultsMetabonome", sep = "")


######dir preparation
moveto_dir <- function(dir="/"){
  dir <- paste(root_dir, dir, sep = "/")
  if(!dir.exists(dir)){
    dir.create(dir, recursive = TRUE)
  }
  setwd(dir = dir)
}

set_db_location(db_dir)

write("#####################Read Table", stdout())
# cpd_table <- read.table(cpd_file, sep = "\t", header = T, row.names = 1, stringsAsFactors = FALSE)
mapping_table <- read.table(map_file, sep = "\t", header = T, row.names = 1, stringsAsFactors = FALSE)


write("#####################Color Assignment", stdout())
pallet <- unique(c(c("#20B2AA", "#DC143C", "#FF9933"), rev(brewer.pal(12,"Paired")),brewer.pal(8,"Set2")[-c(7,8)],brewer.pal(8,"Dark2"),brewer.pal(12,"Set3"),brewer.pal(8,"Accent"),brewer.pal(11,"Spectral")))
group_df <- na.omit(melt(mapping_table[categories]))
all_groups <- unique(group_df[,ncol(group_df)])
all_groups <- sort(all_groups)
group_colors <- pallet[1:length(all_groups)]
names(group_colors) <- all_groups



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

select_cpd_df <- function(cpd_label="name"){
  label <- cpd_df[cpd_label][,1]
  selected <- !(is.na(label)|label=="")
  out_df <- cpd_df[selected, ]
  label <- label[selected]
  out_df <- out_df[,-c(1:9)]
  out_df <- apply(out_df, 2, function(x){return(base::tapply(x, INDEX = label, sum))})
  return(t(out_df))
}

data_factory <- function(category="Category", cpd_label = "name", type = "pair", cmpType=1){
  map_df <- mapping_table[category]
  base_df <- rownames_join(map_df, select_cpd_df(cpd_label = cpd_label))
  base_df <- base_df[order(base_df[,1]), ]
  base_df <- data.frame(SampleID=rownames(base_df), base_df, check.names = FALSE)
  out_list <- list()
  if(type=="pair"){
    groups <- unique(map_df[,1])
    paired_groups <- yield_combine(groups)
    nums <- length(paired_groups)
    cmpType <- choose(length(cmpType)==1, rep(cmpType, nums), cmpType)
    for(i in 1:nums){
      ele_groups <- paired_groups[[i]]
      ele_name <- reduce(ele_groups, paste, sep="_")
      ele_is_paired <- length(ele_groups)==2
      ele_cmpType <- cmpType[i]
      ele_data <- base_df[base_df[,2]%in%ele_groups, ]
      ele_colors <- group_colors[ele_groups]
      out_list[[i]] <- list(data=ele_data, name=ele_name, groups=ele_groups, colors=ele_colors, cmpType=ele_cmpType, paired=ele_is_paired)
    }
    return(out_list)
  }else{
    if(type=="single"){
      groups <- unique(base_df[, 2])
      for (i in 1:length(groups)) {
        out_list[[i]] <- list(name=groups[i], data=base_df[base_df[, 2]==groups[i], ][, -c(1,2)])
      }
    }else{
     out_list <- base_df
    }
  }
  return(out_list)
}

for(category in categories){
  category <- "Category"
  write("#####################Data Prepare", stdout())
  moveto_dir(paste(category, "/00-DataProcess", sep = "/"))
  metadata <- read.csv(metadata_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  map_category <- mapping_table[category]
  metadata <- rownames_join(metadata, map_category)
  colnames(metadata)[3] <- "class"
  write.csv(data.frame(sample=rownames(metadata), metadata, check.names = FALSE), file = "metadata.csv", row.names = FALSE)
  shiftCor("metadata.csv", cpd_file, Frule = 0.8, MLmethod = "QCRFSC", QCspan = 0,imputeM = "KNN", plot = TRUE)
  cpd_table <- read.csv("statTarget/shiftCor/After_shiftCor/shift_sample_cor.csv", check.names = FALSE, header = T, row.names = 1, stringsAsFactors = FALSE)
  cpd_table <- t(cpd_table[, -1])
  cpd_names <-  rownames(cpd_table)
  ref_df <- cross_reference(cpd_names)
  rownames(ref_df)<-ref_df[,1]
  cpd_df <- rownames_join(ref_df, cpd_table)
  colnames(cpd_df)[1]<-"name"
  cpd_df<-cpd_df[,-c(2,3)]
  write.csv(cpd_df, "compound_cross_reference_after_QC.csv", row.names = F)
  dir.create("01-QC")
  system("mv statTarget/shiftCor/* 01-QC")
  # system("mv statTarget/shiftCor/RSDresult/* 01-QC")
  system("rm -r statTarget")
  
  
  write("#####################Data Normlization", stdout())
  moveto_dir(paste(category, "/00-DataProcess/02-Normlization", sep = "/"))
  df <- data_factory(type = "all")
  mSet <- prepare(df)
  mSet <- PlotNormSummary(mSetObj = mSet, "compound_wise_normlization.pdf", format = "pdf", dpi = 300, width = 30)
  mSet <- PlotSampleNormSummary(mSetObj = mSet, "sample_wise_normlization.pdf", format = "pdf", dpi = 300, width = 30)
  
  
  datas <- data_factory(category = category, cpd_label = "class")
  for(data in datas){
    write("#####################ConcentrationSummary", stdout())
    moveto_dir(paste(category, "01-ConcentrationSummary/1-Barplot", data$name, sep = "/"))
    abundance_barplot(data$data, colors = data$colors, by_group_mean = FALSE, prefix = "Compounds_with_Biological_Roles_")   
    abundance_barplot(data$data, colors = data$colors, by_group_mean = TRUE, prefix = "Compounds_with_Biological_Roles_group_mean_")
  }
  datas <- data_factory(category = category)
  report_group <- datas[[1]]$name
  for(data in datas){
    write("#####################ConcentrationSummary", stdout())
    moveto_dir(paste(category, "01-ConcentrationSummary/1-Barplot", data$name, sep = "/"))
    abundance_barplot(data$data, colors = data$colors, by_group_mean = FALSE)   
    abundance_barplot(data$data, colors = data$colors, by_group_mean = TRUE, prefix = "group_mean_")   

    moveto_dir(paste(category, "01-ConcentrationSummary/2-Heatmap", data$name, sep = "/"))
    abundance_heatmap(data$data, colors = data$colors, by_group_mean = FALSE, cluster = FALSE)   
    abundance_heatmap(data$data, colors = data$colors, by_group_mean = FALSE, cluster = TRUE, prefix = "clustered_")   
    abundance_barplot(data$data, colors = data$colors, by_group_mean = TRUE, prefix = "group_mean_")  
    
    write("#####################PCA Analysis", stdout())
    moveto_dir(paste(category, "02-PCA", sep = "/"))
    PCA(data$data, out_dir = data$name, colors = data$colors)
    
    write("#####################PLSDA Analysis", stdout())
    moveto_dir(paste(category, "03-PLSDA", sep = "/"))
    PLSDA(data$data, out_dir = data$name, colors = data$colors)
    
    write("#####################OPLSDA Analysis", stdout())
    moveto_dir(paste(category, "04-OPLSDA", sep = "/"))
    OPLSDA(data$data, out_dir = data$name, colors = data$colors)   

    write("#####################FoldChange Analysis", stdout())
    moveto_dir(paste(category, "05-FoldChange", sep = "/"))
    if(data$paired){
      Volcano(data$data, out_dir = data$name, colors = data$colors, cmpType = data$cmpType)
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
  
  datas <- data_factory(category = category, cpd_label = "kegg_id")
  for(data in datas){
    write("#####################PathwayTopoEnrichment Analysis ", stdout())
    moveto_dir(paste(category, "08-PathwayTopoEnrichment", data$name, sep = "/"))
    pathway_analysis(data$data, out_dir = "./", threshp = 0.2, organism = organism)
    if(file.exists("selected_compounds.csv")){
      df <- read.csv("selected_compounds.csv", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
      ora.df <- read.csv("pathway_enrichment_and_topo_results.csv", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
      gc <- groupCenter(data$data)
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
  }
  datas <- data_factory(category = category, type = "single")
  for (data in datas) {
    write("#####################Correlation Analysis", stdout())
    moveto_dir(paste(category, "09-CorrelationAnalysis", data$name, sep = "/"))
    group_cor_heatmap(data$data, cluster = FALSE)
    group_cor_heatmap(data$data, cluster = TRUE, prefix = "clustered_")
  }
  data <- data_factory(category = category, type = "all")
  data <- data[, -c(1, 2)]
  moveto_dir(paste(category, "09-CorrelationAnalysis", sep = "/"))
  group_cor_heatmap(data, cluster = FALSE, prefix = "all_group_")
  group_cor_heatmap(data, cluster = TRUE, prefix = "all_group_clustered_")
  moveto_dir(category)
  system(sprintf("cp -rp %s/FiguresTablesForReport ./", base_dir))
  system("mv FiguresTablesForReport/结题报告.html ./")
  # system(sprintf("cp 08-PathwayTopoEnrichment/%s/%s_dpi300.pdf FiguresTablesForReport/图9-3.pdf", report_group, str_replace_all(REPORT_PATHWAY, " +", "_")))
  system(sprintf("bash %s/prepare_report_images.sh %s", base_dir, report_group))
}


