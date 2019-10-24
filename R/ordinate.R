library("MetaboAnalystR")
library("scatterplot3d")
library("purrr")
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


prepare<-function(table){
    mSet<-InitDataObjects("pktable", "stat", FALSE)
    ## [1] "MetaboAnalyst R objects initialized ..."
    mSet<-Read.TextData(mSet, table, "colu", "disc")
    mSet<-SanityCheckData(mSet)
    mSet<-ReplaceMin(mSet);
    mSet<-FilterVariable(mSet, "iqr", "F", 25)
    ## [1] " Further feature filtering based on Interquantile Range"
    mSet<-PreparePrenormData(mSet)
    mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "AutoNorm", ratio=FALSE, ratioNum=20)
    return(mSet)

}

remove_files<-function(files_pattern=c("t_test", "loadings", "coef", "vip")){
    for(pattern in files_pattern){
        for(f in list.files()){
            if(length(grep(pattern, f))>0){
                file.remove(f)
            }
        }        
    }
}


prepare_out_dir<-function(out_dir, file_prefix){
    PWD<<-getwd()
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
    return(data.frame(x, y[match(rownames(x), rownames(y)), ], check.rows=TRUE, check.names=FALSE))
}


Scatter3D<-function(xyz, color_by, to_pdf="Scatter3D_plot.pdf", colors=NA, shapes=NA, axis_labels=c("Axis.1","Axis.2","Axis.3"), plot_line=FALSE, box=TRUE){
    ######pcoa 3d plot
    pdf(to_pdf, width=8.4, height=6.4)
    opar<-par(no.readonly=TRUE)
    par(fig=c(0,0.8,0,1))
    scatterplot3d(xyz,mar=c(2.2,2.2,0,0)+1, xlab=axis_labels[1], ylab=axis_labels[2],zlab=axis_labels[3],color=asign_category_map(color_by, colors), grid=TRUE, box=box, type=choose(plot_line,"h","p"), lty.hplot=2, pch=choose(any(is.na(shapes)), 19, asign_category_map(color_by, shapes)))
    par(fig=c(0.8,1,0,1),xpd=TRUE)
    legend("center", legend = sort(unique(color_by)),bty = 'n',xpd = TRUE,horiz = FALSE,col = or(colors, rainbow(length(unique(color_by)))), pch = or(shapes, 19), inset = -0.1)
    par(opar)
    dev.off()
    # write.table(as.matrix(tdata), PCoA_ordtxtname, quote=FALSE, col.names=NA, sep="\t")
}


PCA<-function(table, out_dir, colors=NA){
    mSet<-prepare(table)
    prepare_out_dir(out_dir)
    mSet<-tryCatch(
        {
            UpdateGraphSettings(mSet, colVec=colors, shapeVec=NA)
            mSet<-PCA.Anal(mSet)
            mSet<-PlotPCA2DScore(mSet, "pca_score2d_", "pdf", 300, width=30, 1,2,0.95,0,0)
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
            mSet<-Ttests.Anal(mSet = mSet, nonpar = F, threshp = 2)
            mSet<-PLSR.Anal(mSet, reg=TRUE)
            mSet<-PLSDA.CV(mSet, methodName = "T", compNum = 2, choice = "Q2")
            mSet<-PLSDA.Permut(mSet)
            remove_files()
            mSet<-PlotPLS2DScore(mSet, "plsda_score2d_", "pdf", 300, width=30, 1,2,0.95,0,0)
            mSet<-PlotPLS3DScoreImg(mSet,"plsda_score3d_","pdf", 300, 30, 1,2,3, 40)
            mSet<-PlotPLS.Permutation(mSet = mSet, imgName="plsda_permutation_", format = "pdf", dpi = 300, width = 30)
            # browser()
            imp<-reduce(list(mSet$analSet$plsr$loadings[,1:3], mSet$analSet$plsda$vip.mat, mSet$analSet$plsda$coef.mat, mSet$analSet$tt$sig.mat), rownames_join)
            imp<-imp[, -5]
            colnames(imp)[1:5]<-c("comp1.loadings", "comp2.loadings", "comp3.loadings", "Vip", "Coefficients")
            write.csv(imp, "plsda_features_importance.csv")

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
            mSet<-PlotSPLS3DScoreImg(mSet,"splsda_score3d_","pdf", 300, 30, 1,2,3, 40, box=TRUE)
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
            mSet<-Ttests.Anal(mSet = mSet, nonpar = F, threshp = 2)
            mSet<-OPLSR.Anal(mSet, reg=TRUE)
            remove_files()            
            mSet<-OPLSDA.Permut(mSet, num=100)
            mSet<-PlotOPLS2DScore(mSet, "oplsda_score2d_", "pdf", 300, width=30, 1,2,0.95,0,0)
            PlotOPLS.Permutation(mSet,"oplsda_permutation_",format="pdf", dpi=300, width=30)
            imp<-data.frame(mSet$analSet$oplsda$vipVn, mSet$analSet$oplsda$coefficients, check.rows=T)
            colnames(imp)<-c("Vip", "Coefficients");
            load.mat <- cbind(mSet$analSet$oplsda$loadingMN[,1], mSet$analSet$oplsda$orthoLoadingMN[,1])
            colnames(load.mat) <- c("Loading (t1)","OrthoLoading (to1)")
            imp <- reduce(list(load.mat, imp, mSet$analSet$tt$sig.mat), rownames_join)
            write.csv(imp, "oplsda_features_importance.csv")
            write.csv(mSet$analSet$oplsda$modelDF, file="oplsda_model_fitness_summary.csv")
            
        },
        error=function(e){print(e)},
        finally={setwd(PWD)}
    )
}