###############################################################################
# Script D_BI_with6: Visualization + Enrichment for BI (No LPS/Saline), top-6
###############################################################################
options(java.parameters = "-Xmx64g")
library(ggplot2)
library(tidyr)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(pheatmap)

scriptB_bi_with6 <- "resultsBI_with6_B"
finalColors_bi <- readRDS(file.path(scriptB_bi_with6,"bi_noLPSsaline_moduleColors.rds"))
MEs_bi <- readRDS(file.path(scriptB_bi_with6,"bi_noLPSsaline_eigengenes.rds"))
sampleInfoBI <- readRDS(file.path(scriptB_bi_with6,"bi_noLPSsaline_sample_info.rds"))

scriptC_bi_with6 <- "resultsBI_with6_C"
biNES <- readRDS(file.path(scriptC_bi_with6,"biNES_noLPSsaline_with6.rds"))

modules_bi <- setdiff(unique(finalColors_bi),"grey")

# Bubble heatmap of NES
biContrasts <- c("X1d_BI - X1d_sham","X3d_BI - X1d_sham","X7d_BI - X1d_sham")
bi_data <- do.call(rbind, lapply(modules_bi, function(mod){
  x <- biNES[[mod]]
  if(is.null(x)) return(NULL)
  svals <- sapply(biContrasts, function(cc) if(!is.null(x[[cc]])) x[[cc]] else NA)
  data.frame(Module=mod, t(svals))
}))
if(!is.null(bi_data) && nrow(bi_data)>0){
  colnames(bi_data)[-1] <- biContrasts
  bi_data_long <- tidyr::pivot_longer(bi_data, cols=-Module, names_to="TimePoint", values_to="NES")
  bi_data_long <- bi_data_long[!is.na(bi_data_long$NES), ]
  p_bubble <- ggplot(bi_data_long, aes(x=TimePoint, y=Module)) +
    geom_point(aes(size=abs(NES), color=NES)) +
    scale_size_area(max_size=15) +
    scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0) +
    theme_minimal() +
    labs(title="BI Module Enrichment", size="|NES|", color="NES")
  ggsave("BI_module_enrichment_bubble_with6.png", p_bubble, width=8, height=6)
}

# Single time-course NES (all modules)
nes_list <- list()
for(mod in modules_bi){
  x <- biNES[[mod]]; if(is.null(x)) next
  df_mod <- data.frame(
    Module=mod, Contrast=biContrasts,
    NES=sapply(biContrasts, function(cc) if(!is.null(x[[cc]])) x[[cc]] else NA)
  )
  df_mod <- df_mod[!is.na(df_mod$NES), ]
  nes_list[[mod]] <- df_mod
}
all_nes_df <- do.call(rbind, nes_list)
if(!is.null(all_nes_df) && nrow(all_nes_df)>0){
  all_nes_df$Contrast <- factor(all_nes_df$Contrast, levels=biContrasts)
  p_bi_allline <- ggplot(all_nes_df, aes(x=Contrast, y=NES, group=Module, color=Module)) +
    geom_line() + geom_point() + theme_minimal() +
    labs(title="BI NES Time Course (all modules)", x="Contrast", y="NES")
  ggsave("BI_NES_timecourse_allModules_with6.png", p_bi_allline, width=7, height=5)
}

# ME–ME correlation heatmap
if(!is.null(MEs_bi) && ncol(MEs_bi)>1){
  cor_mat <- cor(MEs_bi, use="pairwise.complete.obs")
  pheat <- pheatmap(cor_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                    main="Module Eigengene Correlation (BI)", fontsize=10)
  ggsave("BI_module_correlation_pattern_with6.png", plot=pheat, width=7, height=6)
}

# Gene lists + GO/KEGG
for(mod in modules_bi){
  modGenes <- names(finalColors_bi)[ finalColors_bi==mod ]
  write.table(modGenes, file=paste0("bi_with6_module_",mod,"_genes.txt"),
              row.names=FALSE, col.names=FALSE, quote=FALSE)
  gene.df <- bitr(modGenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  egIDs <- unique(gene.df$ENTREZID)
  if(length(egIDs) >= 5){
    ego <- enrichGO(gene=egIDs, OrgDb=org.Mm.eg.db, keyType="ENTREZID", ont="BP",
                    pAdjustMethod="BH", qvalueCutoff=0.05, readable=TRUE)
    if(nrow(as.data.frame(ego))>0){
      write.csv(as.data.frame(ego), paste0("BI_",mod,"_GO_BP_with6.csv"), row.names=FALSE)
    }
    ekegg <- enrichKEGG(gene=egIDs, organism='mmu', pAdjustMethod='BH', qvalueCutoff=0.05)
    if(nrow(as.data.frame(ekegg))>0){
      write.csv(as.data.frame(ekegg), paste0("BI_",mod,"_KEGG_with6.csv"), row.names=FALSE)
    }
  }
}

message("Script D_BI_with6 done.")
