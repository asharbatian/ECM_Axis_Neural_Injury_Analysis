###############################################################################
# Script C_SCI: Differential Expression + GSEA for SCI
###############################################################################
options(java.parameters = "-Xmx64g")
library(limma)
library(fgsea)

scriptA_sci <- "results_scriptA_SCI"
exprSCI <- readRDS(file.path(scriptA_sci,"sci_data.rds"))

scriptB_sci <- "results_scriptB_SCI"
finalColors <- readRDS(file.path(scriptB_sci,"sci_wgcna_moduleColors.rds"))
sampleInfoSCI <- readRDS(file.path(scriptB_sci,"sci_sample_info.rds"))

# 1. DE for SCI
sciLevels <- c("Control","X30min","X4h","X1d","X3d","X7d","X28d")
sampleInfoSCI$Time <- factor(make.names(sampleInfoSCI$Time), levels=sciLevels)
designSCI <- model.matrix(~0 + Time, data=sampleInfoSCI)
colnames(designSCI) <- sciLevels

fitSCI <- lmFit(exprSCI, designSCI)
contrastsSCI <- makeContrasts(
  "X30min-Control"=X30min-Control,
  "X4h-Control"   =X4h-Control,
  "X1d-Control"   =X1d-Control,
  "X3d-Control"   =X3d-Control,
  "X7d-Control"   =X7d-Control,
  "X28d-Control"  =X28d-Control,
  levels=designSCI
)
fitSCI2 <- eBayes(contrasts.fit(fitSCI, contrastsSCI))
sciContrasts <- colnames(contrastsSCI)

# 2. GSEA
modules <- setdiff(unique(finalColors),"grey")

sciNES <- list()
for(con in sciContrasts){
  res <- topTable(fitSCI2, coef=con, number=Inf, sort.by="none")
  ranks <- sort(setNames(res$t, rownames(res)), decreasing=TRUE)
  for(mod in modules){
    geneSet <- names(finalColors)[ finalColors==mod ]
    geneSet <- intersect(geneSet, names(ranks))
    if(length(geneSet)<5) next
    fgRes <- fgsea(pathways=list(mySet=geneSet), stats=ranks, nperm=10000)
    sciNES[[mod]][[con]] <- fgRes$NES[1]
  }
}

scriptC_sci <- "results_scriptC_SCI"
dir.create(scriptC_sci, showWarnings=FALSE, recursive=TRUE)
saveRDS(sciNES, file.path(scriptC_sci,"sciNES.rds"))

message("Script C_SCI done. DE + GSEA for SCI saved in ", scriptC_sci)
