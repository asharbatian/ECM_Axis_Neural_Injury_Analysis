###############################################################################
# Script D_SCI: Visualization + Enrichment for SCI
###############################################################################
options(java.parameters = "-Xmx64g")
library(ggplot2)
library(tidyr)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(pheatmap)

scriptB_sci <- "results_scriptB_SCI"
finalColors <- readRDS(file.path(scriptB_sci,"sci_wgcna_moduleColors.rds"))
MEs_sci     <- readRDS(file.path(scriptB_sci,"sci_wgcna_eigengenes.rds"))

scriptC_sci <- "results_scriptC_SCI"
sciNES <- readRDS(file.path(scriptC_sci,"sciNES.rds"))

modules <- setdiff(unique(finalColors),"grey")

# Bubble heatmap of NES
sciContrasts <- c("X30min-Control","X4h-Control","X1d-Control","X3d-Control","X7d-Control","X28d-Control")
sci_data <- do.call(rbind, lapply(modules, function(mod){
  x <- sciNES[[mod]]
  if(is.null(x)) return(NULL)
  vals <- sapply(sciContrasts, function(cc) if(!is.null(x[[cc]])) x[[cc]] else NA)
  data.frame(Module=mod, t(vals))
}))
if(!is.null(sci_data) && nrow(sci_data)>0){
  colnames(sci_data)[-1] <- sciContrasts
  sci_long <- tidyr::pivot_longer(sci_data, cols=-Module, names_to="TimePoint", values_to="NES")
  sci_long <- sci_long[!is.na(sci_long$NES), ]
  p_sci_bubble <- ggplot(sci_long, aes(x=TimePoint, y=Module)) +
    geom_point(aes(size=abs(NES), color=NES)) +
    scale_size_area(max_size=15) +
    scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0) +
    theme_minimal() +
    labs(title="SCI Module Enrichment", size="|NES|", color="NES")
  ggsave("SCI_module_enrichment_bubble.png", p_sci_bubble, width=8, height=6)
}

# Single time-course NES (all modules)
nes_list <- list()
for(mod in modules){
  x <- sciNES[[mod]]; if(is.null(x)) next
  df_mod <- data.frame(
    Module=mod, Contrast=sciContrasts,
    NES=sapply(sciContrasts,function(cc) if(!is.null(x[[cc]])) x[[cc]] else NA)
  )
  df_mod <- df_mod[!is.na(df_mod$NES), ]
  nes_list[[mod]] <- df_mod
}
all_sci_nes <- do.call(rbind, nes_list)
if(!is.null(all_sci_nes) && nrow(all_sci_nes)>0){
  all_sci_nes$Contrast <- factor(all_sci_nes$Contrast, levels=sciContrasts)
  p_sci_allline <- ggplot(all_sci_nes, aes(x=Contrast, y=NES, group=Module, color=Module))+
    geom_line() + geom_point()+ theme_minimal()+
    labs(title="SCI NES Time Course (all modules)", x="Contrast", y="NES")
  ggsave("SCI_NES_timecourse_allModules.png", p_sci_allline, width=7, height=5)
}

# ME–ME correlation heatmap
if(!is.null(MEs_sci) && ncol(MEs_sci)>1){
  cor_mat <- cor(MEs_sci, use="pairwise.complete.obs")
  p_out <- pheatmap(cor_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                    main="Module Eigengene Correlation (SCI)", fontsize=10)
  ggsave("SCI_module_correlation_pattern.png", plot=p_out, width=7, height=6)
}

# Gene lists + GO/KEGG
for(mod in modules){
  modGenes <- names(finalColors)[ finalColors==mod ]
  write.table(modGenes, file=paste0("sci_module_",mod,"_genes.txt"),
              row.names=FALSE, col.names=FALSE, quote=FALSE)
  gene.df <- bitr(modGenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  egIDs <- unique(gene.df$ENTREZID)
  if(length(egIDs) < 5) next
  ego <- enrichGO(gene=egIDs, OrgDb=org.Mm.eg.db, keyType="ENTREZID", ont="BP",
                  pAdjustMethod="BH", qvalueCutoff=0.05, readable=TRUE)
  if(nrow(as.data.frame(ego))>0){
    write.csv(as.data.frame(ego), paste0("SCI_",mod,"_GO_BP.csv"), row.names=FALSE)
  }
  ekegg <- enrichKEGG(gene=egIDs, organism='mmu', pAdjustMethod='BH', qvalueCutoff=0.05)
  if(nrow(as.data.frame(ekegg))>0){
    write.csv(as.data.frame(ekegg), paste0("SCI_",mod,"_KEGG.csv"), row.names=FALSE)
  }
}

message("Script D_SCI done.")
