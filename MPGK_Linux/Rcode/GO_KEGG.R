library(ieugwasr)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(ggplot2)
library(enrichplot)
library(RColorBrewer)
library(openxlsx)
library(Hmisc)

args <- commandArgs(TRUE)
wd <- args[1]
dis_id <- args[2]
p_cutoff <- args[3]

setwd(wd)

df_dis <- read.table(dis_id, header = TRUE, sep = " ")
df_select <- subset(df_dis, pval < p_cutoff)
SNPs <- c(df_select$SNP)
RSIDinfo <- variants_rsid(rsid = SNPs)
gene_list <- RSIDinfo$geneinfo
gene_lists <- gene_list[-which(gene_list == ".")]
genes <- strsplit(gene_lists, "\\|")
n <- length(genes)
gene_ll <- c()
for (i in 1:n){
  gene_ll <- append(gene_ll, genes[[i]][1])
}
g <- strsplit(gene_ll, ":")
gene_df <- data.frame(matrix(unlist(g), nrow = length(g), byrow = TRUE))
colnames(gene_df) <- c("gene", "id")
genes <- c(gene_df$gene)
genes_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

#GO analysis
go_all <- enrichGO(gene = genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = 'ALL', pAdjustMethod = 'fdr', pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
go_all_df <- go_all@result
go_cc <- enrichGO(gene = genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = 'CC', pAdjustMethod = 'fdr', pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
go_cc_df <- go_cc@result
go_bp <- enrichGO(gene = genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = 'BP', pAdjustMethod = 'fdr', pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
go_bp_df <- go_bp@result
go_mf <- enrichGO(gene = genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = 'MF', pAdjustMethod = 'fdr', pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
go_mf_df <- go_mf@result
sheets <- list("GO" = go_all_df, "CC" = go_cc_df, "BP" = go_bp_df, "MF" = go_mf_df)
write.xlsx(sheets, "./GO/go_enrich.xlsx")

#GO visualization
#CC
cc_dag <- goplot(go_cc) # enrichplot package
ggsave("CC_DAG.png", plot = cc_dag, path = wd, dpi = 300)
ego_cc <- pairwise_termsim(go_cc)
cc_ego <- emapplot(ego_cc, color = "p.adjust", layout = "kk", showCategory = 20)
ggsave("CC_map.png", plot = cc_ego, path = wd, dpi = 300)

#BP
bp_dag <- goplot(go_bp) # enrichplot package
ggsave("BP_DAG.png", plot = bp_dag, path = wd, dpi = 300)
ego_bp <- pairwise_termsim(go_bp)
bp_ego <- emapplot(ego_bp, color = "p.adjust", layout = "kk", showCategory = 20)
ggsave("BP_map.png", plot = bp_ego, path = wd, dpi = 300)

#MF
mf_dag <- goplot(go_mf) # enrichplot package
ggsave("MF_DAG.png", plot = mf_dag, path = wd, dpi = 300)
ego_mf <- pairwise_termsim(go_mf)
mf_ego <- emapplot(ego_mf, color = "p.adjust", layout = "kk", showCategory = 20)
ggsave("MF_map.png", plot = mf_ego, path = wd, dpi = 300)

#barplot
cc <- go_cc_df[1:20,]
bp <- go_bp_df[1:20,]
mf <- go_mf_df[1:20,]

cc$Description <- capitalize(cc$Description)
bp$Description <- capitalize(bp$Description)
mf$Description <- capitalize(mf$Description)

CC <- cc
CC$Description <- str_trunc(CC$Description, width = 60, side = "right")
BP <- bp
BP$Description <- str_trunc(BP$Description, width = 60, side = "right")
MF <- mf
MF$Description <- str_trunc(MF$Description, width = 60, side = "right")

go_theme <- theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), legend.title = element_text(size = 12), legend.text = element_text(size = 10))
CC$term <- factor(CC$Description, levels = rev(CC$Description))
BP$term <- factor(BP$Description, levels = rev(BP$Description))
MF$term <- factor(MF$Description, levels = rev(MF$Description))

GO_Bar <- function(x){
  y <- get(x)
  ggplot(data = y, aes(x = Count, y = term, fill = -log10(pvalue))) + scale_y_discrete(labels = function(y) str_wrap(y, width = 50)) + geom_bar(stat = "identity", width = 0.8) + labs(x = "Gene Number", y = "Description", title = paste0(x, " of GO enrichment barplot")) + theme_bw() + go_theme
}
p_cc <- GO_Bar("CC") + scale_fill_distiller(palette = "YlOrRd", direction = 1)
p_bp <- GO_Bar("BP") + scale_fill_distiller(palette = "YlOrBr", direction = 1)
p_mf <- GO_Bar("MF") + scale_fill_distiller(palette = "YlGnBu", direction = 1)
ggsave("CC_bar.png", plot = p_cc, path = wd)
ggsave("BP_bar.png", plot = p_bp, path = wd)
ggsave("MF_bar.png", plot = p_mf, path = wd)

#dotplot
rf_cc <- apply(CC, 1, function(x){
  GeneRatio <- eval(parse(text = x["GeneRatio"]))
  BgRatio <- eval(parse(text = x["BgRatio"]))
  RF_CC <- round(GeneRatio/BgRatio, 2)
  RF_CC
})
CC$Rich_Factor <- rf_cc
rf_bp <- apply(BP, 1, function(x){
  GeneRatio <- eval(parse(text = x["GeneRatio"]))
  BgRatio <- eval(parse(text = x["BgRatio"]))
  RF_BP <- round(GeneRatio/BgRatio, 2)
  RF_BP
})
BP$Rich_Factor <- rf_bp
rf_mf <- apply(MF, 1, function(x){
  GeneRatio <- eval(parse(text = x["GeneRatio"]))
  BgRatio <- eval(parse(text = x["BgRatio"]))
  RF_MF <- round(GeneRatio/BgRatio, 2)
  RF_MF
})
MF$Rich_Factor <- rf_mf

GO_Dot <- function(x){
  y = get(x)
  ggplot(data = y, aes(x = Rich_Factor, y = term)) + geom_point(aes(size = Count, color = -log10(pvalue))) + scale_y_discrete(labels = function(y) str_wrap(y, width = 50)) + labs(x = "Rich Factor", y = "Description", title = paste0(x, " of GO enrichment Dotplot"), size = "Gene Number") + theme_bw() + go_theme
}
dp_cc <- GO_Dot("CC") + scale_color_distiller(palette = "YlOrRd", direction = 1)
dp_bp <- GO_Dot("BP") + scale_color_distiller(palette = "YlOrBr", direction = 1)
dp_mf <- GO_Dot("MF") + scale_color_distiller(palette = "YlGnBu", direction = 1)
ggsave("CC_dot.png", plot = dp_cc, path = wd)
ggsave("BP_dot.png", plot = dp_bp, path = wd)
ggsave("MF_dot.png", plot = dp_mf, path = wd)

#KEGG analysis 
kegg <- enrichKEGG(gene = genes_entrez$ENTREZID, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
kegg_df <- kegg@result
sheet <- list("KEGG" = kegg_df)
write.xlsx(sheet, "./KEGG/kegg_enrich.xlsx")

#KEGG visualization
KEGG <- kegg_df[1:20,]

KEGG$pathway <- factor(KEGG$Description, levels = rev(KEGG$Description))
kegg_theme <- theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), legend.title = element_text(size = 12), legend.text = element_text(size = 10))
p_kegg <- ggplot(data = KEGG, aes(x = Count, y = pathway, fill = -log10(pvalue))) + scale_fill_distiller(palette = "YlGnBu", direction = 1) + geom_bar(stat = "identity", width = 0.8) + theme_bw() + labs(x = "Number of Gene", y = "Pathway", title = "KEGG enrichment barplot") + kegg_theme
dp_kegg <- ggplot(data = KEGG, aes(x = Count, y = reorder(pathway, Count))) + geom_point(aes(size = Count, color = -log10(pvalue))) + scale_color_distiller(palette = "YlGnBu", direction = 1) +  theme_bw() + labs(x = "Number of Gene", y = "Pathway", title = "KEGG enrichment dotplot", size = "Count") + kegg_theme
ggsave("KEGG_bar.png", plot = p_kegg, path = wd)
ggsave("KEGG_dot.png", plot = dp_kegg, path = wd)

rm(list = ls())
