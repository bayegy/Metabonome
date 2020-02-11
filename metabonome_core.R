library("MetaboAnalystR")
library("scatterplot3d")
library("purrr")
library("ggplot2")
library("ggrepel")
library("stringr")
library("pheatmap")
library("psych")
library("ggpubr")
# library("ggbiplot")
#'util
choose<-function(condition,choice1,choice2){
  if(condition){
    return(choice1)
  }else{
    return(choice2)
  }
}

#'util
or<-function(choice1,choice2){
  if(length(choice1)==1){
      if(is.na(choice1)|choice1==FALSE){
        return(choice2)
      }  
  }else{
      if(length(choice1)==0){
        return(choice2)
      }
  }
  return(choice1)
}

up_to_down<-function(m){
  d2<-ncol(m)
  for(i in 1:d2){
    for (j in 1:d2) {
      if(i>j){
        m[i,j]<-m[j,i]
      }
    }
  }
  return(m)
}

my_comb <- function(set, n){
  list_append<-function(ls, ele){
    len <- length(ls)
    ls[[len+1]] <- ele
    return(ls)
  }
  env1 <- new.env()
  env1$results<-list()
  env1$current_vector <- c()
  len <- length(set)
  ceil <- len - n + 1
  combinations<-function(l, j){
    if(l==n){
      env1$results<-list_append(env1$results, env1$current_vector)
      return(NULL)
    }
    for(i in j:(ceil+l)){
      env1$current_vector <- c(env1$current_vector, set[i])
      combinations(l+1, i+1)
      env1$current_vector <- env1$current_vector[-length(env1$current_vector)]
    }
  }
  combinations(0, 1)
  return(env1$results)
}


set_db_location <- function(location){
  DB_LOCATION <<- normalizePath(location)
}


require_db <- function(dbs){
  file.copy(paste(DB_LOCATION, "/", dbs, sep = ""), "./")
}

prepare<-function(table, anal.type="stat"){
    mSet<-InitDataObjects("conc", anal.type = anal.type, FALSE)
    ## [1] "MetaboAnalyst R objects initialized ..."
    mSet<-Read.TextData(mSet, table, "row", "disc")
    mSet<-SanityCheckData(mSet)
    mSet<-ReplaceMin(mSet);
    # mSet<-FilterVariable(mSet, "iqr", "F", 25)
    ## [1] " Further feature filtering based on Interquantile Range"
    mSet<-PreparePrenormData(mSet)
    mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "AutoNorm", ratio=FALSE, ratioNum=20)
    return(mSet)
}

remove_files<-function(files_pattern=c("t_test", "loadings", "coef", "vip"), mode = "exclude"){
  if(mode == "exclude"){
    for(pattern in files_pattern){
      for(f in list.files()){
        if(grepl(pattern, f)){
          file.remove(f)
        }
      }        
    }
  }else{
    
    grepin <- function(string, patterns){
      for(pattern in patterns){
        if(grepl(pattern, string)){
          return(TRUE)
        }
      }
      return(FALSE)
    }
    
    for(f in list.files()){
      if(!grepin(f, files_pattern)){
        file.remove(f)
      }
    }
  }
}


prepare_out_dir<-function(out_dir, file_prefix){
    PWD <<- getwd()
    if(!dir.exists(out_dir)){dir.create(out_dir,recursive = T)}
    setwd(out_dir)
}


asign_category_map<-function(category, map=NA){
    map_out<-as.character(category)
    unig<-sort(unique(category))
    groups_map<-or(map, rainbow(length(unig)))
    for (i in 1:length(unig)) {
        map_out[map_out==unig[i]]<-groups_map[i]
    }
    return(map_out)
}


rownames_join<-function(x,y){
    # print(x)
    # print(y)
    return(data.frame(x, y[match(rownames(x), rownames(y)), ], check.rows=TRUE, check.names=FALSE))
}



cross_reference<-function(query_compounds, type="name", mask=c("\\'", "\\s")){
  norm_for_match<-function(x){
    x<-tolower(x)
    for(s in mask){
      x<-gsub(s,"",x)
    }
    return(x)
  }
  cpd_map<-readRDS(paste(DB_LOCATION, "/compound_map_file.rds", sep = ""))
  cpd_class <- read.table(paste(DB_LOCATION, "/br08001_first_level.txt", sep = ""), sep = "\t", header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
  cpd_class <- cpd_class[!duplicated(cpd_class[, 1]), ]
  class_cpd <- cpd_class[, 2]
  names(class_cpd) <- cpd_class[, 1]
  # query_table<-read.table(query_table, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings = "")
  if(type=="name"){
    query_name<-sapply(query_compounds, norm_for_match)
    cpd_name<-sapply(cpd_map[, 2], norm_for_match)
    cpd_map<-rbind(cpd_map, NA)
    nr<-nrow(cpd_map)
    find_index<-function(q){
      for(i in 1:length(cpd_name)){
        for(ob in str_split(cpd_name[i],";")[[1]]){
          if(ob==q){
            return(i)
          }
        }
      }
      return(nr)
    }
    matched_index<-sapply(query_name, find_index)
  }else{
    matched_index<-match(query_compounds, cpd_map[type][,1])
  }

  find_results<- data.frame(cpd_map[matched_index, ], check.names=FALSE)
  find_results$class <- class_cpd[find_results["kegg_id"][,1]]
  # find_results$class[is.na(find_results$class)] <- "Others"
  return(find_results)
}


scatter3D<-function(xyz, color_by, out_img="Scatter3D_plot.pdf", colors=NA, shapes=NA, axis_labels=c("Axis.1","Axis.2","Axis.3"), plot_line=FALSE, box=TRUE){
    ######pcoa 3d plot
    pdf(out_img, width=8.4, height=6.4)
    opar<-par(no.readonly=TRUE)
    par(fig=c(0,0.8,0,1))
    scatterplot3d(xyz,mar=c(2.2,2.2,0,0)+1, xlab=axis_labels[1], ylab=axis_labels[2],zlab=axis_labels[3],color=asign_category_map(color_by, colors), grid=TRUE, box=box, type=choose(plot_line,"h","p"), lty.hplot=2, pch=choose(any(is.na(shapes)), 19, asign_category_map(color_by, shapes)))
    par(fig=c(0.8,1,0,1),xpd=TRUE)
    legend("center", legend = sort(unique(color_by)),bty = 'n',xpd = TRUE,horiz = FALSE,col = or(colors, rainbow(length(unique(color_by)))), pch = or(shapes, 19), inset = -0.1)
    par(opar)
    dev.off()
    # write.table(as.matrix(tdata), PCoA_ordtxtname, quote=FALSE, col.names=NA, sep="\t")
}


plot_volcano<-function(data, out_img="volcano.pdf", axis=c(1,2), top=NA, threshold=c(2, 1.3), threshold_both=TRUE, label="rownames", x_lab=NA, y_lab=NA, bg_colors=c("#E1FFFF05", "#99FF9905", "#FFCC9905"), size="integ"){
    xy<-choose(is.numeric(axis), colnames(data)[axis], axis)
    lbl<-choose(is.numeric(label), choose(label, colnames(data)[label], "rownames"), label)
    sz <- choose(is.numeric(size), colnames(data)[size], size)
    x_index<-which(colnames(data)%in%xy[1])
    y_index<-which(colnames(data)%in%xy[2])
    colnames(data)[x_index]<-"plot_volcano_x"
    colnames(data)[y_index]<-"plot_volcano_y"
    span_x<-0.5*sd(data[,x_index])
    span_y<-0.5*sd(data[,y_index])
    d_max_x<-max(abs(data[,x_index]))+span_x
    d_max_y<-max(abs(data[,y_index]))+span_y
    d_min_x<-min(data[, x_index])
    d_min_y<-min(data[, y_index])
    d_min_x<-choose(d_min_x<0, d_min_x-span_x, d_min_x)
    d_min_y<-choose(d_min_y<0, d_min_y-span_y, d_min_y)
    ifelse(lbl=="rownames", data$plot_volcano_lbl<-rownames(data), colnames(data)[colnames(data)%in%lbl]<-"plot_volcano_lbl")
    if(size=="integ"){
      data$plot_volcano_size <- abs(data$plot_volcano_x) + abs(data$plot_volcano_y)
    }else{
      if(!is.na(size)){
        colnames(data)[colnames(data)%in%sz]<-"plot_volcano_size"
      }
    }
    data <- data[order(abs(data[, x_index]), abs(data[, y_index]), decreasing = T),]
    label_data<-choose(is.na(top),
                       choose(threshold_both,
                              data[(abs(data[x_index])>or(threshold[1], 0))&(abs(data[y_index])>or(threshold[2], 0)), ],
                              data[(abs(data[x_index])>or(threshold[1], 0))|(abs(data[y_index])>or(threshold[2], 0)), ]),
                       head(data, top))
    p<-ggplot(data=data,aes(x=plot_volcano_x, y=plot_volcano_y)) +
      geom_hline(yintercept=threshold[2],linetype="dashed") + 
      choose(d_min_y>=0, theme(), geom_hline(yintercept=-threshold[2],linetype="dashed")) +
      geom_vline(xintercept=threshold[1],linetype="dashed") +
      choose(d_min_x>=0, theme(), geom_vline(xintercept=-threshold[1],linetype="dashed")) +
      choose(d_min_x>=0, theme(), geom_vline(xintercept = 0)) +
      choose(d_min_y>=0, theme(), geom_hline(yintercept = 0)) +
      choose(is.na(threshold[1]), theme(), geom_rect(xmin=threshold[1], xmax=d_max_x, ymin=d_min_y, ymax=d_max_y, fill=bg_colors[1])) +
      choose(d_min_x>=0, theme(), geom_rect(xmin=-threshold[1], xmax=-d_max_x, ymin=d_min_y, ymax=d_max_y, fill=bg_colors[1])) +
      choose(is.na(threshold[2]), theme(), geom_rect(ymin=threshold[2], ymax=d_max_y, xmin=d_min_x, xmax=d_max_x, fill=bg_colors[2])) +
      choose(d_min_y>=0, theme(), geom_rect(ymin=-threshold[2], ymax=-d_max_y, xmin=d_min_x, xmax=d_max_x, fill=bg_colors[2])) +
      choose(any(is.na(threshold)), theme(), geom_rect(ymin=threshold[2], ymax=d_max_y, xmin=threshold[1], xmax=d_max_x, fill=bg_colors[3])) +
      choose(any(is.na(threshold))|d_min_x>=0, theme(), geom_rect(ymin=threshold[2], ymax=d_max_y, xmin=-threshold[1], xmax=-d_max_x, fill=bg_colors[3])) +
      choose(any(is.na(threshold))|d_min_y>=0, theme(), geom_rect(ymin=-threshold[2], ymax=-d_max_y, xmin=threshold[1], xmax=d_max_x, fill=bg_colors[3])) +
      choose(any(is.na(threshold))|d_min_y>=0|d_min_x>=0, theme(), geom_rect(ymin=-threshold[2], ymax=-d_max_y, xmin=-threshold[1], xmax=-d_max_x, fill=bg_colors[3])) +
      choose(is.na(size), geom_point(size=3, color="#DC143C", alpha=0.3), geom_point(aes(size=plot_volcano_size), color="#DC143C", alpha=0.3)) +
      geom_text_repel(data=label_data, aes(x=plot_volcano_x,y=plot_volcano_y,label=plot_volcano_lbl), color='black', size=3) +
      theme_bw() + xlab(or(x_lab, xy[1])) + ylab(or(y_lab, xy[2])) + labs(size = size) +
      choose(size=="integ", guides(size=FALSE), theme()) +
      scale_y_continuous(limits=c(d_min_y, d_max_y), expand = c(0, 0)) +
      scale_x_continuous(limits=c(d_min_x, d_max_x), expand = c(0, 0)) +
      theme(panel.grid=element_blank(),
            axis.line = element_line(),
            panel.border =  element_blank(),
            text = element_text(size = 15))
    ggsave(plot = p, out_img, dpi=300, height = 8, width = choose(is.na(size), 8, 10))
}


PCA<-function(table, out_dir, colors=NA){
    mSet<-prepare(table)
    prepare_out_dir(out_dir)
    mSet<-tryCatch(
        {
            UpdateGraphSettings(mSet, colVec=colors, shapeVec=NA)
            mSet<-PCA.Anal(mSet)
            mSet<-PlotPCA2DScore(mSet, "pca_score2d_", "pdf", 300, width=30, 1,2,0.95,show=0, grey.scale =0)
            mSet<-PlotPCA2DScore(mSet, "pca_score2d_with_ids_", "pdf", 300, width=30, 1,2,0.95,show=1, grey.scale =0)
            # mSet<-PlotPCABiplot(mSetObj=mSet, imgName="pca_score_biplot_", format="png", dpi=300, width=30, 1, 2)
            # pca <-  mSet$analSet$pca
            # labels<-paste("PC", c(1, 2, 3), " (", 100*round(pca$variance[c(1, 2, 3)], 3), "%)", sep="")
            # Scatter3D(pca$x[, c(1, 2, 3)], mSet$dataSet$cls, "pca_score3d_plot.pdf", colors=colors, axis_labels=labels)
            mSet<-PlotPCA3DScoreImg(mSet,"pca_score3d_","pdf", 300, 30, 1,2,3, 40)
        },
        error=function(e){print(e)},
        finally={setwd(PWD)}
    )
}


PLSDA<-function(table, out_dir, colors=NA){
    mSet<-prepare(table)
    prepare_out_dir(out_dir)
    mSet<-tryCatch(
        {
            UpdateGraphSettings(mSet, colVec=colors, shapeVec=NA)
            mSet<-Ttests.Anal(mSet = mSet, nonpar = F, threshp = 2, equal.var = F)
            mSet<-PLSR.Anal(mSet, reg=TRUE)
            # 
            tryCatch(
              {
                mSet<-PLSDA.CV(mSet, methodName = "L", compNum = 2, choice = "Q2")
                mSet<-PLSDA.Permut(mSet, type="accu")
                mSet<-PlotPLS.Permutation(mSet = mSet, imgName="plsda_permutation_", format = "pdf", dpi = 300, width = 30)
              },
              error=function(e){print("Permutation failed!");print(e)}
            )            
            # browser()
            remove_files()
            mSet<-PlotPLS2DScore(mSet, "plsda_score2d_", "pdf", 300, width=30, 1,2,0.95,show=0,0)
            mSet<-PlotPLS2DScore(mSet, "plsda_score2d_with_ids_", "pdf", 300, width=30, 1,2,0.95,show=1,0)
            mSet<-PlotPLS3DScoreImg(mSet,"plsda_score3d_","pdf", 300, 30, 1,2,3, 40)
            # 
            imp<-purrr::reduce(list(mSet$analSet$plsr$loadings[,1:3], mSet$analSet$plsda$vip.mat, mSet$analSet$plsda$coef.mat, mSet$analSet$tt$sig.mat), rownames_join)
            imp<-imp[, -5]
            colnames(imp)[1:5]<-c("comp1.loadings", "comp2.loadings", "comp3.loadings", "VIP", "Coefficients")
            imp["-log10(FDR adjusted p)"] <- -log10(imp["FDR"])
            plot_volcano(imp, out_img = "plsda_features_importance_fdr_adjusted.pdf", axis = c("VIP","-log10(FDR adjusted p)"), threshold=c(1, 1.3))
            plot_volcano(imp, out_img = "plsda_features_importance.pdf", axis = c("VIP","-log10(p)"), threshold=c(1, 1.3))
            write.csv(imp, "plsda_features_importance.csv")
            # print(imp)
        },
        error=function(e){print(e)},
        finally={setwd(PWD)}
    )
}


SPLSDA<-function(table, out_dir, colors=NA){
    mSet<-prepare(table)
    prepare_out_dir(out_dir)
    mSet<-tryCatch(
        {
            UpdateGraphSettings(mSet, colVec=colors, shapeVec=NA)
            mSet<-SPLSR.Anal(mSet,c(50,50,50))
            mSet<-PlotSPLS2DScore(mSet, "splsda_score2d_", "pdf", 300, width=30, 1,2,0.95,0,0)
            mSet<-PlotSPLS3DScoreImg(mSet,"splsda_score3d_","pdf", 300, 30, 1,2,3, 40)
            write.csv(mSet$analSet$splsr$error.rate, "model_error_rate.csv")
        },
        error=function(e){print(e)},
        finally={setwd(PWD)}
    )
}



OPLSDA<-function(table, out_dir, colors=NA){
    mSet<-prepare(table);
    prepare_out_dir(out_dir);
    mSet<-tryCatch(
        {
            UpdateGraphSettings(mSet, colVec=colors, shapeVec=NA)
            mSet<-Ttests.Anal(mSet = mSet, nonpar = F, threshp = 2, equal.var = F)
            mSet<-OPLSR.Anal(mSet, reg=TRUE)
            remove_files()
            mSet<-PlotOPLS2DScore(mSet, "oplsda_score2d_", "pdf", 300, width=30, 1,2,0.95,show=0,0)
            mSet<-PlotOPLS2DScore(mSet, "oplsda_score2d_with_ids_", "pdf", 300, width=30, 1,2,0.95,show=1,0)
            # tryCatch(
            #  {
            mSet<-OPLSDA.Permut(mSet, num=100)
            PlotOPLS.Permutation(mSet,"oplsda_permutation_",format="pdf", dpi=300, width=30)
            #  },
            #  error=function(e){print("Permutation failed!");print(e)}
            #)
            # browser()
            imp<-data.frame(mSet$analSet$oplsda$vipVn, mSet$analSet$oplsda$coefficients, check.rows=T)
            colnames(imp)<-c("VIP", "Coefficients");
            load.mat <- cbind(mSet$analSet$oplsda$loadingMN[,1], mSet$analSet$oplsda$orthoLoadingMN[,1])
            colnames(load.mat) <- c("Loading (t1)","OrthoLoading (to1)")
            imp <- reduce(list(load.mat, imp, mSet$analSet$tt$sig.mat), rownames_join)
            imp["-log10(FDR adjusted p)"] <- -log10(imp["FDR"])
            plot_volcano(imp, out_img = "oplsda_features_importance_fdr_adjusted.pdf", axis = c("VIP","-log10(FDR adjusted p)"), threshold=c(1, 1.3))
            plot_volcano(imp, out_img = "oplsda_features_importance.pdf", axis = c("VIP","-log10(p)"), threshold=c(1, 1.3))
            write.csv(imp, "oplsda_features_importance.csv")
            write.csv(mSet$analSet$oplsda$modelDF, file="oplsda_model_fitness_summary.csv")
            
        },
        error=function(e){print(e)},
        finally={setwd(PWD)}
    )
}


# RF<-function(table, out_dir, colors=NA){
#   mSet<-prepare(table);
#   prepare_out_dir(out_dir);
#   mSet<-tryCatch(
#     {
#       UpdateGraphSettings(mSet, colVec=colors, shapeVec=NA)
#       mSet<-Ttests.Anal(mSet = mSet, nonpar = F, threshp = 2, equal.var = F)
#       mSet<-RF.Anal(mSetObj = mSet)
#       remove_files(c("csv"))
#       imp <- reduce(list(mSet$analSet$rf$importance, mSet$analSet$tt$sig.mat), rownames_join)
#       # imp["-log10(FDR adjusted p)"] <- -log10(imp["FDR"])
#       # plot_volcano(imp, out_img = "RF_features_importance.pdf", axis = c("MeanDecreaseAccuracy","-log10(FDR adjusted p)"))
#       write.csv(imp, "RF_features_importance.csv")
#       mSet<-PlotRF.Classify(mSetObj=mSet, imgName="RF_performance", format="pdf", dpi=300, width=30)
#       mSet<-PlotRF.Outlier(mSetObj=mSet, imgName="RF_outliers", format="pdf", dpi=300, width=30)
#       mSet<-PlotRF.VIP(mSetObj=mSet, imgName="RF_features_importance", format="pdf", dpi=300, width=30)
#     },
#     error=function(e){print(e)},
#     finally={setwd(PWD)}
#   )
# }


groupCenter<-function(table, method="mean", check.names=TRUE){
  df <- choose(is.character(table), read.csv(table, quote = "", check.names = F, stringsAsFactors = F), table)
  rownames(df) <- df[, 1]
  gp <- df[, 2]
  df <- df[,-c(1,2)]
  df <- t(base::apply(df, 2, function(x){return(base::tapply(x, INDEX = gp, FUN = choose(method=="mean", base::mean, base::median), na.rm=T))}))
  colnames_cp <- colnames(df)
  colnames(df)<-paste(colnames_cp, choose(method=="mean","-Mean","-Median"), sep = "")
  df<-data.frame(df, check.names = FALSE)
  df["Which-Max"]<-colnames_cp[apply(df, 1, which.max)]
  if(check.names){
    for(m in c("\\(", "\\)", "\\+", "\\[", "\\]")){
      rownames(df) <- str_replace_all(rownames(df), m, "")
    }
  }
  return(df)
}


Volcano <- function(table, out_dir="./", cmpType=0, colors=NA){
  mSet<-prepare(table);
  gp_mean<-groupCenter(table = table);
  prepare_out_dir(out_dir);
  mSet<-tryCatch(
    {
      UpdateGraphSettings(mSet, colVec=colors, shapeVec=NA)
      # mSet<-Ttests.Anal(mSet = mSet, nonpar = F, threshp = 2)
      # mSet<-FC.Anal.unpaired(mSetObj = mSet)
      mSet<-Volcano.Anal(mSet, threshp = 2, fcthresh = 1, cmpType = cmpType, equal.var = F)
      print("################################################check the following compounds names if error occured")
      dif_name <- rownames(gp_mean)[!rownames(gp_mean)%in%rownames(mSet$analSet$volcano$sig.mat)]
      print(dif_name)
      imp <- reduce(list(mSet$analSet$volcano$sig.mat, gp_mean), rownames_join)
      write.csv(imp, "Volcano_features_significance.csv")
      mSet<-Volcano.Anal(mSet, threshp = 0.05, fcthresh = 2, cmpType = cmpType)
      remove_files(c("volcano.csv"))
      # imp["-log10(FDR adjusted p)"] <- -log10(imp["FDR"])
      plot_volcano(imp, out_img = "Volcano_features_importance.pdf", axis = c("log2(FC)","-log10(p)"), threshold = c(1, 1.3))
      PlotVolcano(mSetObj=mSet, imgName="Volcano_features_significance_", plotLbl=1, format="pdf", dpi=300, width=30)
    },
    error=function(e){print(e)},
    finally={setwd(PWD)}
  )
}


RF<-function(table, out_dir){
  mSet<-prepare(table);
  prepare_out_dir(out_dir);
  mSet<-tryCatch(
    {
      # UpdateGraphSettings(mSet, colVec=colors, shapeVec=NA)
      mSet<-Ttests.Anal(mSet = mSet, nonpar = F, threshp = 2, equal.var = F)
      mSet<-RF.Anal(mSetObj = mSet)
      remove_files(c("csv"))
      imp <- reduce(list(mSet$analSet$rf$importance, mSet$analSet$tt$sig.mat), rownames_join)
      # imp["-log10(FDR adjusted p)"] <- -log10(imp["FDR"])
      # plot_volcano(imp, out_img = "RF_features_importance.pdf", axis = c("MeanDecreaseAccuracy","-log10(FDR adjusted p)"))
      write.csv(imp, "RF_features_importance.csv")
      # mSet<-PlotRF.Classify(mSetObj=mSet, imgName="RF_performance", format="pdf", dpi=300, width=30)
      mSet<-PlotRF.Outlier(mSetObj=mSet, imgName="RF_outliers", format="pdf", dpi=300, width=30)
      mSet<-PlotRF.VIP(mSetObj=mSet, imgName="RF_features_importance", format="pdf", dpi=300, width=30)
    },
    error=function(e){print(e)},
    finally={setwd(PWD)}
  )
}



SVM<-function(table, out_dir){
  mSet<-prepare(table);
  prepare_out_dir(out_dir);
  mSet<-tryCatch(
    {
      # UpdateGraphSettings(mSet, colVec=colors, shapeVec=NA)
      mSet<-Ttests.Anal(mSet = mSet, nonpar = F, threshp = 2, equal.var = F)
      # Set the biomarker analysis mode to perform Multivariate exploratory ROC curve analysis ("explore")
      mSet<-SetAnalysisMode(mSet, "explore")
      # Prepare data for biomarker analysis
      mSet<-PrepareROCData(mSet)
      # Perform multivariate ROC curve analysis, using SVM classification and ranking
      mSet<-PerformCV.explore(mSet, cls.method = "svm", rank.method = "svm", lvNum = 2)
      ### OPTION 1 Comparison plot of ROC curves of all models ###
      mSet<-PlotROC(mSet, imgName = "svm_roc_all_models_", format="png", dpi=300, mdl.inx= 0, avg.method = "threshold", show.conf = 0, show.holdout = 0, focus="fpr", cutoff=0.5)
      # Plot predicted class probabilities for each sample for a selected model, not showing labels of wrongly classified samples
      # mSet<-PlotProbView(mSet, imgName = "multi_roc_prob_", format="pdf", dpi=300, width=30, mdl.inx = -1, show = 1, showPred = 0)
      # Plot the predictive accuracy of models with increasing number of features
      # mSet<-PlotAccuracy(mSet, imgName = "svm_accuracy_", format="png", dpi=300)
      # Plot the most important features of a selected model ranked from most to least important
      mSet<-PlotImpVars(mSet, imgName = "svm_features_importance_", format="png", dpi=300, mdl.inx = -1, measure="mean", feat.num=15)
      remove_files(c("csv"))
      imp <- reduce(list(mSet$analSet$multiROC$imp.mat, mSet$analSet$tt$sig.mat), rownames_join)
      # imp["-log10(FDR adjusted p)"] <- -log10(imp["FDR"])
      # plot_volcano(imp, out_img = "RF_features_importance.pdf", axis = c("MeanDecreaseAccuracy","-log10(FDR adjusted p)"))
      write.csv(imp, "SVM_features_importance.csv")
    },
    error=function(e){print(e)},
    finally={setwd(PWD)}
  )
}


enrichment_from_peak<-function(table, out_dir, organism="hsa"){
  # database_path<-normalizePath(database)
  database_path<- paste(DB_LOCATION, "/", organism, "/", organism, "_mfn.rds", sep = "")
  # database_file<-str_extract(database, "[^/]+$")
  # database_name<-str_extract(database_file, "^[^\\.]+")
  database_name <- paste(organism, "_mfn", sep = "") 
  mSet<-prepare(table);
  prepare_out_dir(out_dir);
  mSet<-tryCatch(
    {
      # Re-perform normalization, without auto-scaling
      # mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)
      # Perform t-test
      mSet<-Ttests.Anal(mSet = mSet, nonpar = F, threshp = 0.25, equal.var = F)
      # Convert the output of the t-test to the proper input for mummichog analysis
      mSet<-Convert2Mummichog(mSet)
      find_mummi_input<-function(){
        for(file in list.files()){
          if(grepl("mummichog_input",file)){
            return(file)
          }
        }
      }
      
      # Read in the ranked peak list
      SetPeakFormat("mpt")
      mSet<-UpdateInstrumentParameters(mSet, "5", "negative");
      mSet<-Read.PeakListData(mSet, find_mummi_input());
      mSet<-SanityCheckMummichogData(mSet)
      
      write("######################Perform original mummichog algorithm", stdout())
      dir.create("Mummichog",recursive = T)
      setwd("Mummichog")
      file.copy(from=database_path, to="./")
      # Perform original mummichog algorithm
      mSet<-SetPeakEnrichMethod(mSet, "mum")
      mSet<-SetMummichogPval(mSet, 0.25)
      mSet<-PerformPSEA(mSet, database_name, permNum = 1000)
      # mSet<-PlotPeaks2Paths(mSet, "peaks_to_paths_mummichog_", "pdf", 300, width=30)
      mummi.mat <- mSet$mummi.resmat
      plot_data<-data.frame(x=mummi.mat[,3]/mummi.mat[,4], y=-log10(mummi.mat[,5]))
      colnames(plot_data)<-c("Enrichment Factor", "-log10(p)")
      plot_data<-data.frame(mummi.mat, plot_data, check.names = F)
      plot_volcano(plot_data, axis = c("Enrichment Factor", "-log10(p)"), out_img = "mummichog_enriched_pathways.pdf", threshold = c(NA, 1.3))
      remove_files(c(database_name, ".json", "mummichog_pathway_enrichment.csv"))
      write.csv(plot_data, "mummichog_pathway_enrichment.csv")
      
      write("######################Perform GSEA", stdout())
      dir.create("../GSEA",recursive = T)
      setwd("../GSEA")
      file.copy(from=database_path, to="./")
      # Perform GSEA
      mSet<-SetPeakEnrichMethod(mSet, "gsea")
      mSet<-PerformPSEA(mSet, database_name, permNum = 1000)
      # mSet<-PlotPeaks2Paths(mSet, "peaks_to_paths_GSEA_", "pdf", 300, width=30)
      gsea.mat <- mSet$mummi.gsea.resmat
      plot_data<-data.frame(x = gsea.mat[,2]/gsea.mat[,1], y = -log10(gsea.mat[,3]))
      colnames(plot_data)<-c("Enrichment Factor", "-log10(p)")
      plot_data<-data.frame(gsea.mat, plot_data, check.names = F)
      plot_volcano(plot_data, axis = c("Enrichment Factor", "-log10(p)"), out_img = "gsea_enriched_pathways.pdf", threshold = c(NA, 1.3))
      remove_files(c(database_name,".json", "mummichog_fgsea_pathway_enrichment.csv"))
      write.csv(plot_data, "gsea_pathway_enrichment.csv")
      
      write("######################Perform GSEA-mummichog", stdout())
      dir.create("../GSEA-Mummichog",recursive = T)
      setwd("../GSEA-Mummichog")
      file.copy(from=database_path, to="./")
      # Perform GSEA
      mSet<-SetPeakEnrichMethod(mSet, "integ")
      mSet<-SetMummichogPval(mSet, 0.25)
      mSet<-PerformPSEA(mSet, database_name)
      combo.resmat <- mSet$integ.resmat
      plot_data<-data.frame(x = -log10(combo.resmat[,5]), y = -log10(combo.resmat[,4]), combo.p <- -log10(combo.resmat[,6]))
      colnames(plot_data)<-c("GSEA -log10(p)", "Mummichog -log10(p)", "Combined -log10(p)")
      plot_data<-data.frame(combo.resmat, plot_data, check.names = F)
      plot_volcano(plot_data, axis = c("GSEA -log10(p)", "Mummichog -log10(p)"), out_img = "gsea_mummichog_enriched_pathways.pdf", threshold = c(1.3, 1.3), threshold_both=F, size = "Combined -log10(p)")
      # mSet<-PlotIntegPaths(mSet, "integ_paths_", "pdf", 300, width=30)
      remove_files(c(database_name,".json","csv"))
      write.csv(plot_data, "gsea_mummichog_pathway_enrichment.csv")
      # return(mSet)
    },
    error=function(e){print(e)},
    finally={setwd(PWD)}
  )
}



enrichment_from_compounds <- function(table, out_dir, threshp = 0.05){
  table <- normalizePath(table)
  mSet<-prepare(table);
  prepare_out_dir(out_dir);
  compound_db_path <- paste(DB_LOCATION, "/", "compound_db.rds", sep = "")
  kegg_pathway_path <- paste(DB_LOCATION, "/", "kegg_pathway.rda", sep = "")
  mSet<-tryCatch(
    {
      # UpdateGraphSettings(mSet, colVec=colors, shapeVec=NA)
      mSet<-Ttests.Anal(mSet = mSet, nonpar = F, threshp = threshp, equal.var = F)
      df <- read.csv(table, stringsAsFactors = F, check.names = F)
      df <- df[, c(1,2,which(colnames(df)%in%rownames(mSet$analSet$tt$sig.mat)))]

      write("######################Perform QEA", stdout())
      dir.create("QEA",recursive = T)
      setwd("QEA")
      file.copy(c(compound_db_path, kegg_pathway_path), "./")
      # browser()
      mSet <- prepare(df, "msetqea")
      mSet<-CrossReferencing(mSet, "kegg")
      mSet<-CreateMappingResultTable(mSet)
      mSet<-SanityCheckData(mSet);
      mSet<-SetMetabolomeFilter(mSet, F);
      mSet<-SetCurrentMsetLib(mSet, "kegg_pathway", 2);
      mSet<-CalculateGlobalTestScore(mSet)
      mSet<-PlotQEA.Overview(mSet, "QEA_global_test_score_", "bar", "pdf", 300, width=30)
      # mSet<-PlotQEA.Overview(mSet, "QEA_network_", "net", "pdf", 300, width=30)
      file.remove(c("compound_db.rds", "kegg_pathway.rda"))
      
      write("######################Perform ORA", stdout())
      dir.create("../ORA", recursive = T)
      setwd("../ORA")
      file.copy(c(compound_db_path, kegg_pathway_path), "./")
      tmp.vec <- colnames(df)[-c(1,2)]
      mSet<-InitDataObjects("conc", "msetora", FALSE)
      mSet<-Setup.MapData(mSet, tmp.vec);
      mSet<-CrossReferencing(mSet, "kegg");
      # Create the mapping results table
      mSet<-CreateMappingResultTable(mSet)
      # Set the metabolite filter
      mSet<-SetMetabolomeFilter(mSet, F);
      # Select metabolite set library
      mSet<-SetCurrentMsetLib(mSet, "kegg_pathway", 2);
      # Calculate hypergeometric score, results table generated in your working directory
      mSet<-CalculateHyperScore(mSet)
      # Plot the ORA, bar-graph
      mSet<-PlotORA(mSet, "ORA_global_test_score_", "bar", "pdf", 300, width=30)
      # mSet<-PlotORA(mSet, "ORA_network_", "net", "pdf", 300, width=30)
      file.remove(c("compound_db.rds", "kegg_pathway.rda"))
    },
    error=function(e){print(e)},
    finally={setwd(PWD)}
  ) 
}


topo_analysis <- function(table, out_dir, threshp = 0.5, organism="hsa", cmpType=0){
  table <- normalizePath(table)
  mSet<-prepare(table);
  prepare_out_dir(out_dir);
  mSet<-tryCatch(
    {
      # UpdateGraphSettings(mSet, colVec=colors, shapeVec=NA)
      mSet<-Volcano.Anal(mSet, threshp = threshp, fcthresh = 1.2, cmpType = cmpType)
      imp <- mSet$analSet$volcano$sig.mat
      data_mat <- data.frame(rownames(imp), imp[, 2])
      colnames(data_mat) <- c("#KEGG",	"logFC")
      write.table(data_mat, file = "integ_cmpds.txt", row.names = F, sep = "\t", quote = FALSE)
    
      write("######################Perform Topo", stdout())
      require_db(c("hsa/hsa.rda", "compound_db.rds"))
      mSet<-InitDataObjects("conc", "pathinteg", FALSE)
      #  # Set organism library
      mSet<-SetOrganism(mSet, org = organism)
      cmpdListFile<-"integ_cmpds.txt"
      cmpdList<-readChar(cmpdListFile, file.info(cmpdListFile)$size)
      mSet<-PerformIntegCmpdMapping(mSet, cmpdList, organism, "kegg");
      mSet<-CreateMappingResultTable(mSet)
      mSet<-PrepareIntegData(mSet);
      mSet<-PerformIntegPathwayAnalysis(mSet, "dc", "hyper", "metab")
      # mSet<-PlotORA(mSet, "ORA_global_test_score_", "bar", "pdf", 300, width=30)
      mSet<-PlotKEGGPath(mSet, "Alanine, aspartate and glutamate metabolism",528, 480, "png", NULL)
      # remove_files(c("topo_result_pathway.csv", "volcano.csv"), mode = "keep")
      plot_volcano(mSet$dataSet$path.df, out_img = "topo_pathway_impact.pdf", axis = c("-log10(p)", "Impact"), threshold = c(1.3, NA))
    },
    error=function(e){print(e)},
    finally={setwd(PWD)}
  )
}


network_map <- function(table, out_dir, threshp = 0.5, organism="hsa", paired=TRUE, cmpType=0){
  table <- normalizePath(table)
  mSet<-prepare(table);
  prepare_out_dir(out_dir);
  mSet<-tryCatch(
    {
      # UpdateGraphSettings(mSet, colVec=colors, shapeVec=NA)
      if(paired){
        mSet<-Volcano.Anal(mSet, threshp = threshp, fcthresh = 1.2, cmpType = cmpType)
        imp <- mSet$analSet$volcano$sig.mat
        data_mat <- data.frame(rownames(imp), imp[, 2])
        colnames(data_mat) <- c("#KEGG",	"logFC")
        write.table(data_mat, file = "integ_cmpds.txt", row.names = F, sep = "\t", quote = FALSE)
      }else{
        mSet<-Ttests.Anal(mSet = mSet, nonpar = F, threshp = threshp, equal.var = F)
        imp <- mSet$analSet$tt$sig.mat
        data_mat <- data.frame(rownames(imp), imp[, 1])
        colnames(data_mat) <- c("#KEGG",	"logFC")
        write.table(data_mat, file = "integ_cmpds.txt", row.names = F, sep = "\t", quote = FALSE)
      }
      
      write("######################Perform Topo", stdout())
      require_db(c("hsa/hsa.rda", "MetPriCNet.sqlite", "compound_db.rds"))
      
      mSet<-InitDataObjects("conc", "network", FALSE)
      #  # Set organism library
      mSet<-SetOrganism(mSet, org = organism)
      cmpdListFile<-"integ_cmpds.txt"
      cmpdList<-readChar(cmpdListFile, file.info(cmpdListFile)$size)
      mSet<-PerformIntegCmpdMapping(mSet, cmpdList, organism, "kegg");
      
      mSet<-CreateMappingResultTable(mSet)
      mSet<-PrepareNetworkData(mSet);
      mSet<-GetNetworkGeneMappingResultTable(mSet)
      mSet<-SearchNetDB(mSet, "pheno", "metabo_metabolites", FALSE, 0.5)
      mSet<-CreateGraph(mSet)
      
      # remove_files(c("topo_result_pathway.csv", "volcano.csv"), mode = "keep")
      # plot_volcano(mSet$dataSet$path.df, out_img = "topo_pathway_impact.pdf", axis = c("-log10(p)", "Impact"), threshold = c(1.3, NA))
    },
    error=function(e){print(e)},
    finally={setwd(PWD)}
  )  
}


pathway_analysis <- function(table="/home/cheng/pipelines/testdir/testbono/human_cachexia.csv", out_dir="./", threshp = 0.05, organism="rno"){
  # table <- normalizePath(table)
  mSet<-prepare(table);
  prepare_out_dir(out_dir);
  mSet<-tryCatch(
    {
      # UpdateGraphSettings(mSet, colVec=colors, shapeVec=NA)
      # mSet<-Volcano.Anal(mSet, threshp = threshp, fcthresh = 1.2, cmpType = cmpType)
      # imp <- mSet$analSet$volcano$sig.mat
      mSet <- Ttests.Anal(mSet = mSet, nonpar = F, threshp = 2, equal.var = F)
      imp <- mSet$analSet$tt$sig.mat
      # browser()
      imp <- imp[imp[, 2]<=threshp, ]
      tmp.vec <- rownames(imp)
      write.csv(imp, file = "t_test.csv")
      file.rename("t_test.csv", "selected_compounds.csv")
      write("######################Perform topo and enrichment analysis", stdout())
      require_db(c(sprintf("%s/%s.rda", organism, organism), "compound_db.rds", "syn_nms.rds"))
      mSet<-InitDataObjects("conc", "pathora", FALSE)
      mSet<-Setup.MapData(mSet, tmp.vec);
      mSet<-CrossReferencing(mSet, "kegg");
      mSet<-CreateMappingResultTable(mSet);
      mSet<-SetKEGG.PathLib(mSet, organism)
      mSet<-SetMetabolomeFilter(mSet, F);
      mSet<-CalculateOraScore(mSet, "dgr", "hyperg")
      ora.df <- mSet$analSet$ora.df
      plot_volcano(ora.df, out_img = "pathway_enrichment_and_topo_impact.pdf", axis = c("-log10(p)", "Impact"), threshold = c(1.3, NA), label = "Description")
      mSet<-PlotORA(mSet, "ora_hyper_score_", "bar", "pdf", 300, width=30)
      sig_pathways <- ora.df["Description"][ora.df["Raw p"]<=0.05]
      sig_pathways <- sort(sig_pathways)
      for(pathway in sig_pathways){
        print(pathway)
        PlotKEGGPath(mSet, pathway,30, 30, "pdf", 301)
      }
      # REPORT_PATHWAY <<- sig_pathways[1];
      remove_files(c("\\.rda", "\\.rds"));
    },
    error=function(e){
      print(e);
      remove_files(c("\\.rda", "\\.rds"));
      },
    finally={setwd(PWD)}
  )
}


generate_span <- function(number_list, start = 0){
  step_num <- c()
  current_num <- start
  for(i in 1:length(number_list)){
    current_num <- current_num + number_list[i];
    step_num[i] <- current_num
  }
  flat <- step_num - number_list
  out_list <- list()
  for(i in 1:length(flat)){
    out_list[[i]] <- c(flat[i], step_num[i])
  }
  return(out_list)
}




abundance_barplot <- function(table, out_dir="./", by_group_mean=FALSE, num=20, prefix="", colors = NA){
  prepare_out_dir(out_dir);
  mSet<-tryCatch(
    {
      rownames(table) <- table[, 1]
      table <- table[!is.na(table[, 2]),]
      table <- table[order(table[, 2]), ]
      
      group <- table[, 2]
      label_order<-rownames(table)
      
      otu <- table[, -c(1,2)]
      if(min(otu, na.rm = T)<0){
        otu <- otu - min(otu)
      }

      sum_abundance<-colSums(otu)
      otu<-otu[, order(sum_abundance,decreasing = T)]
      
      otu <- apply(otu, 1, function(x){x/sum(x)})*100
      otu <- t(otu)
      
      
      if(by_group_mean){
        otu<-data.frame(apply(otu,2,function(x){tapply(x,INDEX = group,mean)}),check.names=F)
        otu<-data.frame(t(apply(otu,1,function(x){x/sum(x)}))*100,check.names=F)
        label_order<-ordered(rownames(otu))
        #print(otu)
      }
      
      if(num<ncol(otu)-1){
        otu<-data.frame(otu[,1:num],Other=apply(otu[,(num+1):ncol(otu)],1,sum),check.names = F,check.rows = T)
      }else{
        otu<-data.frame(otu,check.names = F,check.rows = T)
      }
      
      otu_out<-t(otu)/100
      #otu<-data.frame(otu,group)
      otu$id<-rownames(otu)
      p1<-(max(nchar(colnames(otu)))*0.05+0.3)*ceiling(ncol(otu)/17)+2.5
      
      otu<-otu[rev(colnames(otu))]
      otu<-melt(otu,id.vars = "id")
      ban_width <- 0.7
      pallet<-c(rev(brewer.pal(12,"Paired")),brewer.pal(8,"Set2")[-c(7,8)],brewer.pal(8,"Dark2"),brewer.pal(12,"Set3"),brewer.pal(8,"Accent"),brewer.pal(11,"Spectral"))
      p<-ggplot(otu,aes(x=id,y=value))+geom_bar(mapping = aes(fill=variable),stat = "identity",width = ban_width)+
        guides(fill=guide_legend(title = NULL))+
        scale_fill_manual(values = pallet)+
        scale_x_discrete(limits=label_order)+
        xlab("")+ylab("Relative Abundance(%)")+theme_bw()+
        theme(text = element_text(size = 10),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.line = element_line(),panel.border =  element_blank(),
              axis.text.x = element_text(angle = 45,size = 10,hjust = 1))+
        scale_y_continuous(limits=c(0,101), expand = c(0, 0))
      
      if(!by_group_mean){
        uni_groups <- table(group)
        colors <- or(colors, rainbow(length(uni_groups)))
        span_data <- generate_span(uni_groups, start = 1 - (ban_width/2))
        span_len <- length(span_data)
        annotate_data <- data.frame(matrix(nrow = span_len, ncol = 4))
        for(i in 1:span_len){
          ele <- span_data[[i]]
          ele[2] <- ele[2] - (1- ban_width)
          annotate_data[i, ]<-c(ele, c(mean(ele), 110))
        }
        colnames(annotate_data) <- c("x", "xend", "textx", "texty")
        annotate_data$group <- names(uni_groups)
        p <- p + scale_y_continuous(limits=c(0, 115), expand = c(0, 0), breaks = c(0, 20, 40, 60, 80, 100)) + 
          geom_segment(aes(x=x, xend=xend, y=105, yend=105, colour = group), data = annotate_data,  size = 5, lineend = "butt") +
          scale_colour_manual(values = colors) + guides(colour = FALSE) +
          geom_text(aes(x = textx, y = texty, label = group), data = annotate_data, size = 5)
      }
      
      wd<-length(label_order)*0.2+p1
      wd<-ifelse(wd<50,wd,49.9)
      write.table(otu_out, paste(prefix,"table.txt", sep = ""), sep = "\t",quote=FALSE, row.names = TRUE, col.names = NA)
      ggsave(plot = p, paste(prefix,"barplot.pdf",sep = ""), width = wd, height = 7, dpi = 300)
    },
    error=function(e){print(e)},
    finally={setwd(PWD)}
  )  
}


abundance_heatmap <- function(table, out_dir="./", by_group_mean=FALSE, num=30, prefix="", colors = NA, tile_colors=c("#00FF00","#222222","#FF0000"), trans = FALSE, cluster = TRUE, scale_method="std"){
  prepare_out_dir(out_dir);
  mSet<-tryCatch(
    {
      rownames(table) <- table[, 1]
      table <- table[!is.na(table[, 2]),]
      table <- table[order(table[, 2]), ]
      group <- table[2]
      otu <- table[, -c(1,2)]
      if(min(otu, na.rm = T)<0){
        otu <- otu - min(otu)
      }
      p1<-max(nchar(colnames(otu)))
      sum_otu<-colSums(otu)
      sel<-head(order(sum_otu,decreasing = T), num)
      otu<-otu[, sel]
      if(scale_method=="log"){
        otu<-log(otu+1, base=10)
      }else{
        otu<-scale(otu, center = TRUE, scale = TRUE)
      }
      p2<-3+(0.3*dim(otu)[1])+(0.06*p1)
      p3<-dim(otu)[2]
      if(trans){
        otu<-t(otu)
      }
      ht<-ifelse(trans,3+0.4*p3 ,ifelse(p2<50,p2,49.9))
      wd<-ifelse(trans, ifelse(p2<50,p2,49.9),3+0.4*p3)
      cc<-ifelse(trans,F,T)
      cr<-ifelse(trans,T,F)
      if(cluster){
        cc=TRUE
        cr=TRUE
      }
      annotation_row=NA
      annotation_col=NA
      if(!by_group_mean){
        if(trans){
          annotation_row = NA
          annotation_col = group
        }else{
          annotation_row = group
          annotation_col = NA
        }
      }
      if(by_group_mean){
        if(trans){
          otu<-t(apply(otu,1,function(x){tapply(x,INDEX = group[,1],mean)}))
          p2<-4+(0.3*dim(otu)[2])+(0.06*p1)
          p3<-dim(otu)[1]
          wd<-ifelse(p2<50,p2,49.9)
          ht<-2+0.4*p3
        }else{
          otu<-apply(otu,2,function(x){tapply(x,INDEX = group[,1],mean)})
          p2<-4+(0.3*dim(otu)[1])+(0.06*p1)
          p3<-dim(otu)[2]
          ht<-ifelse(p2<50,p2,49.9)
          wd<-2+0.4*p3
        }
      }
      groups_color <- colors
      names(groups_color) <- sort(unique(group[,1]))
      rexp <- sprintf("annotation_color=list(%s=groups_color)", colnames(group)[1])
      eval(parse(text = rexp))
      write.table(t(otu),paste(prefix,"table.txt",sep = ""),sep = "\t",quote=FALSE,col.names=NA)
      pdf(paste(prefix,"heatmap.pdf",sep = ""), height=ht,width=wd)
      pheatmap(otu,annotation_row=annotation_row,
               annotation_col=annotation_col,fontsize=10,border_color = "black",
               color = colorRampPalette(colors = tile_colors)(100),
               cluster_cols=cc,clustering_distance_cols="euclidean",
               cluster_rows=cr,clustering_distance_rows="euclidean",annotation_colors=annotation_color)
      dev.off()
    },
    error=function(e){print(e)},
    finally={setwd(PWD)}
  )  
}



group_cor_heatmap <- function(table, out_dir="./", prefix="", show_names=FALSE){
  prepare_out_dir(out_dir);
  tryCatch(
    {
      cor.mat <- corr.test(table, method ="pearson", adjust="fdr")
      # browser()
      r_mat <- cor.mat$r
      r_mat <- r_mat[, !colSums(is.na(r_mat))==nrow(r_mat)]
      r_mat <- r_mat[!colSums(is.na(t(r_mat)))==ncol(r_mat), ]
      plot_heatmap <- function(cluster,fix=""){
        pdf(paste(prefix, fix, "pearson_cor_heatmap.pdf", sep = ""), height=10, width=11)
        tryCatch({
          pheatmap(r_mat, show_rownames = show_names,  show_colnames = show_names,
                   color = colorRampPalette(colors = c("green", "#555555", "red"))(100),
                   border_color = NA,
                   cluster_rows=cluster, clustering_distance_rows="euclidean",
                   cluster_cols=cluster, clustering_distance_cols="euclidean")
        }, 
        error=function(e){print(e)},
        finally={dev.off()})
      }
      # browser()
      plot_heatmap(cluster = T, fix = "clustered_")
      plot_heatmap(cluster = F)
      write.table(cor.mat$r,"pearson_cor_matrix.txt", sep="\t", col.names=NA)
      write.table(cor.mat$p,"fdr_adjusted_p_value_matrix.txt", sep="\t", col.names=NA)
    },
    error=function(e){print(e)},
    finally={setwd(PWD)}
  )  
}


get_database<- function(){
  url_for_topo_rda <- "https://www.metaboanalyst.ca/resources/libs/kegg/metabolic/"
  url_for_enrich_rda <- "https://www.metaboanalyst.ca/resources/libs/kegg/hsa.rda"
}


sig_boxplot <- function(table="/home/cheng/pipelines/testdir/testbono/human_cachexia.csv", out_dir="./", threshp = 0.05, colors=NA, top=25, pair_wise=FALSE, width=NA, height=NA){
  # table <- normalizePath(table)
  mSet<-prepare(table);
  prepare_out_dir(out_dir);
  tryCatch(
    {
      # UpdateGraphSettings(mSet, colVec=colors, shapeVec=NA)
      # mSet<-Volcano.Anal(mSet, threshp = threshp, fcthresh = 1.2, cmpType = cmpType)
      # imp <- mSet$analSet$volcano$sig.mat
      mSet <- Ttests.Anal(mSet = mSet, nonpar = F, threshp = 2, equal.var = F)
      imp <- mSet$analSet$tt$sig.mat
      # browser()
      imp <- data.frame(imp, check.names = F)
      imp <- imp[imp[, 2]<=threshp, ]
      sig_features <- NULL
      if(nrow(imp)>0){
        imp <- imp[order(imp[, 2]), ]
        sig_features <- rownames(imp)
        tmp.vec <- head(sig_features, top)
        actual_num <- length(tmp.vec)
        datmat <- data.frame(Category=mSet$dataSet$cls, mSet$dataSet$norm[, tmp.vec])
        df <- melt(datmat, id.vars = "Category")
        groups <- unique(datmat[, 1])
        len_g <- length(groups)
        if(pair_wise|len_g==2){
          method <- "t.test"
        }else{
          method <- "anova"
        }
        remove_files(c("t_test","anova"))
        num_col <- ceiling(sqrt(actual_num))
        num_row <- ceiling(actual_num/num_col)
        p <- ggplot(data = df, mapping = aes(x=Category, y=value, fill=Category)) + geom_boxplot() + theme_bw() +
          stat_compare_means(mapping = aes(group=Category), comparisons = choose(pair_wise, my_comb(groups, 2), NULL),
                             label = "p.signif", label.x = (1+len_g)/2, method = method)+
          facet_wrap(facets = "variable", ncol = num_col) + 
          scale_fill_manual(values = or(colors, rainbow(len_g))) +
          theme(legend.position = "none") +
          labs(x="", y="")
        wd <- num_col * 3
        ht <- num_row * 3.5
        ggsave(filename = "Significant_features_boxplot.pdf", plot = p, width = or(width, wd), height = or(height, ht))
      }
    },
    error=function(e){
      print(e);
    },
    finally={setwd(PWD)}
  )
  return(sig_features)
}


elipse_pca <- function(table="/home/cheng/pipelines/testdir/testbono/human_cachexia.csv", file_biplot="PCA.pdf", file_pc1="PC1.pdf", colors=NA, sam_order=NA, out_dir="./", label=T){
    # table <- normalizePath(table)
    mSet<-prepare(table);
    prepare_out_dir(out_dir);
    mSet<-tryCatch(
      {
        df.pca <- prcomp(mSet$dataSet$norm, center = TRUE, scale. = TRUE)
        # pdf("")
        label <- choose(label, rownames(df.pca$x), NA)
        len_g <- length(unique(mSet$dataSet$cls))
        ss <- apply(df.pca$x, 2, var)
        pct <- ss/sum(ss)
        xylab <- paste(c("PC1 (", "PC2 ("), round(pct[1:2]*100, 2), "%)", sep = "")
        if(!is.na(file_biplot)){
          p<-ggplot(data = data.frame(scale(df.pca$x)), aes(x=PC1, y=PC2))+
            theme_classic()+
            labs(x=xylab[1], y=xylab[2])+
            # geom_abline(data = data.frame(slope=c(0,91), intercept=c(0,0)), mapping = aes(slope=slope, intercept=intercept), size=0.1)+
            geom_hline(yintercept = 0, size=0.1) + geom_vline(xintercept = 0, size=0.1) +
            # scale_x_continuous(limits = c(-3, 3)) + scale_y_continuous(limits = c(-3, 3))+
            geom_point(size=4, mapping = aes(color= mSet$dataSet$cls))+
            geom_text_repel(label=label, size=2) +
            scale_color_manual(values = or(colors, rainbow(len_g))) +
            stat_ellipse(level = 0.95, type = "norm", size=0.1, segments = 300)+
            theme(legend.title = element_text(size = 0), text = element_text(size = 6))
          # browser()
          ggsave(plot=p, file=file_biplot, width=5, height=4.3)
        }
        if(!is.na(file_pc1)){
          pcdf <- df.pca$x
          pcdf <- data.frame(Sample=choose(is.na(sam_order),rownames(pcdf), factor(rownames(pcdf), levels = sam_order)), PC1=pcdf[, 1])
          pc1 <- pcdf[,2]
          pc_mean <- mean(pc1)
          pc_sd <- sd(pc1)
          sds <- c(pc_sd*2, pc_sd*3)
          intercept_y <- c(pc_mean-rev(sds), pc_mean, pc_mean+sds)
          col <- c("red", "green", "black", "green", "red")
          ltype <- c("dashed", "dashed", "solid", "dashed", "dashed")
          text_data <- data.frame(x=4, y=intercept_y+0.5, lab=c("Mean-3SD","Mean-2SD","Mean", "Mean+2SD","Mean+3SD"), size=6)
          # browser()
          p <- ggplot(data = pcdf, mapping = aes(x=Sample, y=PC1)) + 
            geom_point(stat = "identity") +
            geom_hline(linetype=ltype, color=col, yintercept=intercept_y) +
            geom_text(data = text_data, mapping = aes(x=x,y=y,label=lab), size=2) +
            theme_classic2() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                  text = element_text(size = 9))
          w1 <- length(pc1) * 0.1 + 4
          w1 <- ifelse(w1>=50, 50, w1)
          ggsave(filename = file_pc1, plot = p, width = w1, height = 6)
        }
      },
      error=function(e){
        print(e);
      },
      finally={setwd(PWD)}
    )
}
