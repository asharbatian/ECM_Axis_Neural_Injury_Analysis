###############################################################################
# Script B_BI_with6: WGCNA for BI only, with top-6 module limit
###############################################################################
options(java.parameters = "-Xmx64g")
library(WGCNA)
library(limma)

enableWGCNAThreads(nThreads=8)
set.seed(123)

scriptA_bi_with6 <- "resultsBI_with6_A"
bi_data_path <- file.path(scriptA_bi_with6,"bi_noLPSsaline_data.rds")
if(!file.exists(bi_data_path)) stop("Cannot find: ", bi_data_path)

exprBI <- readRDS(bi_data_path)

# sample info => 21 columns => (1d_BI=5,3d_BI=3,7d_BI=3,1d_sham=4,3d_sham=3,7d_sham=3)
biSamples <- colnames(exprBI)
sampleInfoBI <- data.frame(
  GSM = biSamples,
  Time = rep(c("1d_BI","3d_BI","7d_BI","1d_sham","3d_sham","7d_sham"),
             times=c(5,3,3,4,3,3)),
  stringsAsFactors=FALSE
)

# WGCNA
datExpr_bi <- t(exprBI)
gsg <- goodSamplesGenes(datExpr_bi, verbose=3)
if(!gsg$allOK){
  datExpr_bi <- datExpr_bi[gsg$goodSamples, gsg$goodGenes]
}

softPower <- 6
net_bi <- blockwiseModules(
  datExpr_bi,
  power=softPower,
  networkType="signed",
  TOMType="signed",
  corType="bicor",
  deepSplit=2,
  minModuleSize=30,
  reassignThreshold=0.8,
  mergeCutHeight=0.2,
  numericLabels=TRUE,
  pamRespectsDendro=FALSE,
  maxBlockSize=ncol(datExpr_bi),
  verbose=3
)

moduleLabels_bi <- net_bi$colors
moduleColors_bi <- labels2colors(moduleLabels_bi)

uniqueMods_bi <- sort(unique(moduleColors_bi))
nonGrey_bi <- setdiff(uniqueMods_bi,"grey")

finalColors_bi <- moduleColors_bi
# Force top-6
if(length(nonGrey_bi)>6){
  modSizes <- table(moduleColors_bi)
  top6 <- names(sort(modSizes[nonGrey_bi],decreasing=TRUE))[1:6]
  finalColors_bi <- ifelse(moduleColors_bi %in% top6, moduleColors_bi, "grey")
  for(i in seq_along(top6)){
    finalColors_bi[ finalColors_bi==top6[i] ] <- paste0("M",i)
  }
} else {
  # rename them M1..M(n)
  idx <- 1
  for(m in nonGrey_bi){
    finalColors_bi[ finalColors_bi==m ] <- paste0("M", idx)
    idx <- idx+1
  }
}
finalColors_bi[is.na(finalColors_bi)] <- "grey"
names(finalColors_bi) <- colnames(datExpr_bi)

MEList_bi <- moduleEigengenes(datExpr_bi, colors=finalColors_bi)
MEs_bi    <- MEList_bi$eigengenes
# rename columns from "MExx" => simpler
colnames(MEs_bi) <- sub("^ME","", colnames(MEs_bi))

scriptB_bi_with6 <- "resultsBI_with6_B"
dir.create(scriptB_bi_with6, showWarnings=FALSE, recursive=TRUE)

saveRDS(net_bi,              file.path(scriptB_bi_with6,"bi_noLPSsaline_wgcna_net.rds"))
saveRDS(finalColors_bi,      file.path(scriptB_bi_with6,"bi_noLPSsaline_moduleColors.rds"))
saveRDS(MEs_bi,              file.path(scriptB_bi_with6,"bi_noLPSsaline_eigengenes.rds"))
saveRDS(sampleInfoBI,        file.path(scriptB_bi_with6,"bi_noLPSsaline_sample_info.rds"))

cat("Script B_BI_with6 => WGCNA done with top-6 modules => ", scriptB_bi_with6,"\n")
