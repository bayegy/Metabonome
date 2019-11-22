#/usr/bin/env bash

report_group=$1

cpfirst 00-DataProcess/01-QA_QC/QA.pdf FiguresTablesForReport/图3-1.pdf
cpfirst 00-DataProcess/01-QA_QC/PC1.pdf FiguresTablesForReport/图3-2.pdf
cpfirst 00-DataProcess/01-QA_QC/QC.pdf FiguresTablesForReport/图3-3.pdf
# cp '00-DataProcess/01-QA_QC/RSDresult/RSD variation.pdf' FiguresTablesForReport/图3-2.pdf
cp 00-DataProcess/02-Normlization/compound_wise_normlization.pdfdpi300.pdf FiguresTablesForReport/图3-4.pdf
cp 01-ConcentrationSummary/1-Barplot/${report_group}/barplot.pdf  FiguresTablesForReport/图4-1.pdf
cp 01-ConcentrationSummary/1-Barplot/${report_group}/Compounds_with_Biological_Roles_barplot.pdf  FiguresTablesForReport/图4-2.pdf
cp 01-ConcentrationSummary/2-Heatmap/${report_group}/clustered_heatmap.pdf  FiguresTablesForReport/图4-3.pdf
cp 02-PCA/${report_group}/pca_score2d_dpi300.pdf FiguresTablesForReport/图5-1.pdf
cp 03-PLSDA/${report_group}/plsda_score2d_dpi300.pdf FiguresTablesForReport/图6-1.pdf
cp 03-PLSDA/${report_group}/plsda_permutation_dpi300.pdf FiguresTablesForReport/图6-2.pdf
cp 03-PLSDA/${report_group}/plsda_features_importance.pdf FiguresTablesForReport/图6-3.pdf
cp 04-OPLSDA/${report_group}/oplsda_permutation_dpi300.pdf FiguresTablesForReport/图6-4.pdf
cp 05-UnivariateAnalysis/${report_group}/Volcano_features_importance.pdf FiguresTablesForReport/图7-1.pdf
cp 05-UnivariateAnalysis/${report_group}/Significant_features_boxplot.pdf FiguresTablesForReport/图7-2.pdf
cp 06-RandomForest/${report_group}/RF_performancedpi300.pdf FiguresTablesForReport/图8-1.pdf
cp 06-RandomForest/${report_group}/RF_features_importancedpi300.pdf FiguresTablesForReport/图8-2.pdf
cp 07-SupportVectorMachine/${report_group}/svm_accuracy_dpi300.png FiguresTablesForReport/图8-3.png
cp 07-SupportVectorMachine/${report_group}/svm_roc_all_models_dpi300.png FiguresTablesForReport/图8-4.png
cp 07-SupportVectorMachine/${report_group}/svm_features_importance_dpi300.png FiguresTablesForReport/图8-5.png
cp 09-PathwayTopoEnrichment/${report_group}/ora_hyper_score_dpi300.pdf FiguresTablesForReport/图10-1.pdf
cp 09-PathwayTopoEnrichment/${report_group}/pathway_enrichment_and_topo_impact.pdf FiguresTablesForReport/图10-2.pdf
cpfirst 09-PathwayTopoEnrichment/${report_group}/dpi301.pdf FiguresTablesForReport/图10-3.pdf
cpfirst 09-PathwayTopoEnrichment/${report_group}/png FiguresTablesForReport/图10-4.png
cp 08-CorrelationAnalysis/all_group_clustered_pearson_cor_heatmap.pdf FiguresTablesForReport/图9-1.pdf


cd FiguresTablesForReport/
if [ -f 图3-1.pdf ];then echo "Converting pdf to png"; for pdfs in *.pdf; do echo $pdfs; base=$(basename $pdfs .pdf); convert  -density 300 -quality 80 $pdfs ${base}.png; rm $pdfs;done;fi;
