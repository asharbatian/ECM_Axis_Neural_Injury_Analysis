###############################################################################
# Script C_BI_with6: DE + GSEA for BI only, omitting LPS/Saline
###############################################################################
options(java.parameters = "-Xmx64g")
library(limma)
library(fgsea)

scriptA_bi_with6 <- "resultsBI_with6_A"
exprBI <- readRDS(file.path(scriptA_bi_with6,"bi_noLPSsaline_data.rds"))

scriptB_bi_with6 <- "resultsBI_with6_B"
finalColors_bi <- readRDS(file.path(scriptB_bi_with6,"bi_noLPSsaline_moduleColors.rds"))
sampleInfoBI   <- readRDS(file.path(scriptB_bi_with6,"bi_noLPSsaline_sample_info.rds"))

# factor levels => c("X1d_BI","X3d_BI","X7d_BI","X1d_sham","X3d_sham","X7d_sham")
biLevels <- c("X1d_BI","X3d_BI","X7d_BI","X1d_sham","X3d_sham","X7d_sham")
sampleInfoBI$Time <- factor(make.names(sampleInfoBI$Time), levels=biLevels)
designBI <- model.matrix(~0 + Time, data=sampleInfoBI)
colnames(designBI) <- biLevels

fitBI <- lmFit(exprBI, designBI)
# define e.g. (X1d_BI - X1d_sham), (X3d_BI - X1d_sham), ...
contrastsBI <- makeContrasts(
  "X1d_BI - X1d_sham" = X1d_BI - X1d_sham,
  "X3d_BI - X1d_sham" = X3d_BI - X1d_sham,
  "X7d_BI - X1d_sham" = X7d_BI - X1d_sham,
  levels = designBI
)
fitBI2 <- eBayes(contrasts.fit(fitBI, contrastsBI))
biContrasts <- colnames(contrastsBI)

modules_bi <- setdiff(unique(finalColors_bi),"grey")

biNES <- list()
for(con in biContrasts){
  res <- topTable(fitBI2, coef=con, number=Inf, sort.by="none")
  ranks <- sort(setNames(res$t, rownames(res)), decreasing=TRUE)
  for(mod in modules_bi){
    geneSet <- names(finalColors_bi)[ finalColors_bi==mod ]
    geneSet <- intersect(geneSet, names(ranks))
    if(length(geneSet)<5) next
    fgRes <- fgsea(pathways=list(mySet=geneSet), stats=ranks, nperm=10000)
    biNES[[mod]][[con]] <- fgRes$NES[1]
  }
}

scriptC_bi_with6 <- "resultsBI_with6_C"
dir.create(scriptC_bi_with6, showWarnings=FALSE, recursive=TRUE)
saveRDS(biNES, file.path(scriptC_bi_with6,"biNES_noLPSsaline_with6.rds"))

message("Script C_BI_with6 done => DE+GSEA => ", scriptC_bi_with6)
