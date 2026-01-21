###############################################################################
# Script B_SCI: WGCNA for SCI only (forcing top-6 modules)
###############################################################################
options(java.parameters = "-Xmx64g")
library(WGCNA)
library(limma)

enableWGCNAThreads(nThreads=8)
set.seed(123)

scriptA_sci <- "results_scriptA_SCI"
sci_data_path <- file.path(scriptA_sci,"sci_data.rds")
if(!file.exists(sci_data_path)) stop("Cannot find: ", sci_data_path)

exprSCI <- readRDS(sci_data_path)

# sample info
sciSamples <- colnames(exprSCI)
sampleInfoSCI <- data.frame(
  GSM = sciSamples,
  Time = rep(c("Control","30min","4h","1d","3d","7d","28d"),
             times=c(6,9,9,9,9,9,9)),
  stringsAsFactors=FALSE
)

datExpr <- t(exprSCI)  # samples x genes
gsg <- goodSamplesGenes(datExpr, verbose=3)
if(!gsg$allOK){
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

softPower <- 6
net_sci <- blockwiseModules(
  datExpr,
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
  maxBlockSize=ncol(datExpr),
  verbose=3
)

moduleLabels <- net_sci$colors
moduleColors <- labels2colors(moduleLabels)
uniqueMods <- sort(unique(moduleColors))
nonGreyMods<- setdiff(uniqueMods,"grey")

finalColors <- moduleColors
# Force top 6 modules
if(length(nonGreyMods)>6){
  modSizes <- table(moduleColors)
  selMods  <- names(sort(modSizes[nonGreyMods],decreasing=TRUE))[1:6]
  finalColors <- ifelse(moduleColors %in% selMods, moduleColors, "grey")
  for(i in seq_along(selMods)){
    finalColors[ finalColors==selMods[i] ] <- paste0("M",i)
  }
} else {
  # if <=6, rename each color => M1..M(n)
  idx <- 1
  for(m in nonGreyMods){
    finalColors[ finalColors==m ] <- paste0("M", idx)
    idx <- idx+1
  }
}
finalColors[is.na(finalColors)] <- "grey"
names(finalColors) <- colnames(datExpr)

MEList <- moduleEigengenes(datExpr, colors=finalColors)
MEs    <- MEList$eigengenes
# rename columns
colnames(MEs) <- sub("^ME","",colnames(MEs))

scriptB_sci <- "results_scriptB_SCI"
dir.create(scriptB_sci, showWarnings=FALSE, recursive=TRUE)

saveRDS(net_sci,       file.path(scriptB_sci,"sci_wgcna_net.rds"))
saveRDS(finalColors,   file.path(scriptB_sci,"sci_wgcna_moduleColors.rds"))
saveRDS(MEs,           file.path(scriptB_sci,"sci_wgcna_eigengenes.rds"))
saveRDS(sampleInfoSCI, file.path(scriptB_sci,"sci_sample_info.rds"))
cat("Script B_SCI done: WGCNA for SCI only.\n")
