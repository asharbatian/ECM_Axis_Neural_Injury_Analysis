###############################################################################
# Script A_SCI: Fetch and Preprocess (Mouse 430 2.0) for SCI Only
###############################################################################
options(java.parameters = "-Xmx64g")
library(GEOquery)
library(limma)
library(AnnotationDbi)
library(mouse4302.db)  # GPL1261 annotation

scriptA_sci <- "results_scriptA_SCI"
dir.create(scriptA_sci, showWarnings=FALSE, recursive=TRUE)

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
    stop("getGEO() didn't return ExpressionSet or list.")
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

# SCI sample IDs
sci_samples <- c(
  # Control
  "GSM119833","GSM119834","GSM119786","GSM119849","GSM119776","GSM119842",
  # 30min
  "GSM118664","GSM118672","GSM119764","GSM119783","GSM119846","GSM119855",
  "GSM119804","GSM119805","GSM119840",
  # 4h
  "GSM119766","GSM119780","GSM119800","GSM119784","GSM119785","GSM119822",
  "GSM119799","GSM119826","GSM119853",
  # 1d
  "GSM119779","GSM119811","GSM120018","GSM119768","GSM119821","GSM119835",
  "GSM119772","GSM119798","GSM119806",
  # 3d
  "GSM119795","GSM119818","GSM119819","GSM119815","GSM119837","GSM119856",
  "GSM119774","GSM119790","GSM119817",
  # 7d
  "GSM119831","GSM119832","GSM119843","GSM119769","GSM119823","GSM119848",
  "GSM119775","GSM119791","GSM119841",
  # 28d
  "GSM119765","GSM119794","GSM119830","GSM119814","GSM119836","GSM119847",
  "GSM119773","GSM119788","GSM119789"
)

exprSCI_affy <- fetch_geo_data("GSE5296", sci_samples)

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
  
  expr_sub <- expr_mat[bestProbes, ]
  rownames(expr_sub) <- names(bestProbes)
  return(expr_sub)
}

exprSCI <- annotate_and_collapse(exprSCI_affy)
cat("After annotation, SCI dimension:", dim(exprSCI), "\n")

saveRDS(exprSCI, file.path(scriptA_sci, "sci_data.rds"))
message("Script A_SCI done. Preprocessed + Annotated data saved in ", scriptA_sci)
