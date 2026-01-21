Script A — Fetch & Preprocess (01_bi_noLPSsaline_fetch_preprocess.R)
###############################################################################
# Script A_BI_with6: Fetch and Preprocess (Mouse 430 2.0) for BI 
# Omitting "1d_LPS" and "1d_Saline" entirely
###############################################################################
options(java.parameters = "-Xmx64g")
library(GEOquery)
library(limma)
library(AnnotationDbi)
library(mouse4302.db)

scriptA_bi_with6 <- "resultsBI_with6_A"
dir.create(scriptA_bi_with6, showWarnings=FALSE, recursive=TRUE)

if(!dir.exists("data")){
  dir.create("data", showWarnings=FALSE, recursive=TRUE)
}

fetch_geo_data <- function(gse_id, samples, destdir="data"){
  message("Downloading ", gse_id, " ...")
  gset <- tryCatch({
    getGEO(gse_id, GSEMatrix=TRUE, AnnotGPL=TRUE, destdir=destdir)
  }, error=function(e){
    stop("Error fetching ", gse_id, ": ", e$message)
  })
  if(inherits(gset,"ExpressionSet")){
    eset <- gset
  } else if(is.list(gset)){
    esets <- Filter(function(x) inherits(x,"ExpressionSet"), gset)
    if(length(esets)<1)
      stop("No ExpressionSet found for ", gse_id)
    eset <- esets[[1]]
  } else {
    stop("getGEO() returned neither an ExpressionSet nor a list.")
  }
  
  exprs_data <- exprs(eset)
  missing <- setdiff(samples, colnames(exprs_data))
  if(length(missing)>0){
    stop("Some samples not found in ",gse_id,": ",paste(missing,collapse=", "))
  }
  exprs_data <- exprs_data[,samples]
  exprs_norm <- normalizeBetweenArrays(exprs_data, method="quantile")
  return(exprs_norm)
}

# Updated BI sample set => Omit LPS and Saline
bi_samples_noLPSsaline <- c(
  # 1 day Brain Injury
  "GSM866313","GSM866314","GSM866316","GSM866318","GSM866320",
  # 3 day Brain Injury
  "GSM866331","GSM866332","GSM866333",
  # 7 day Brain Injury
  "GSM866334","GSM866335","GSM866336",
  # 1 day sham
  "GSM866315","GSM866317","GSM866319","GSM866321",
  # 3 day sham
  "GSM866337","GSM866338","GSM866339",
  # 7 day sham
  "GSM866340","GSM866341","GSM866342"
)

exprBI_affy <- fetch_geo_data("GSE35338", bi_samples_noLPSsaline)

annotate_and_collapse <- function(expr_mat){
  affyProbes <- rownames(expr_mat)
  probe2sym <- mapIds(mouse4302.db, keys=affyProbes, column="SYMBOL",
                      keytype="PROBEID", multiVals="first")
  variances <- apply(expr_mat,1,var)
  annotDF <- data.frame(PROBE=affyProbes,
                        SYMBOL=ifelse(is.na(probe2sym),NA,probe2sym),
                        VARIANCE=variances,
                        stringsAsFactors=FALSE)
  annotDF <- annotDF[!is.na(annotDF$SYMBOL), ]
  
  bestProbes <- tapply(annotDF$PROBE, annotDF$SYMBOL, function(pg){
    subDF <- annotDF[ annotDF$PROBE %in% pg, ]
    maxVarRow <- subDF[ which.max(subDF$VARIANCE), ]
    maxVarRow$PROBE
  })
  
  expr_sub <- expr_mat[ bestProbes, ]
  rownames(expr_sub) <- names(bestProbes)
  return(expr_sub)
}

exprBI <- annotate_and_collapse(exprBI_affy)
cat("After annotation, BI dimension:", dim(exprBI), "\n")

saveRDS(exprBI, file.path(scriptA_bi_with6,"bi_noLPSsaline_data.rds"))
message("Script A_BI_with6 done => preprocessed & annotated => ", scriptA_bi_with6)
