title: 'Skyrocket'

# Libraries =================
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(progeny)
library(readr)
library(stringr)
library(readxl)
library(harmony)
library(EnhancedVolcano)
library(genekitr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(pheatmap)
library(tradeSeq)
library(monocle)
library(ggpointdensity)
library(slingshot)
library(patchwork)
library(scales)
library(PseudotimeDE)
library(NMF)
library(ggalluvial)

# Loding Data =================

# Total 
load('/DATA07/home/shlee/SKYRocket/SKYRocket_TotalMainAnnotation.robj')

# CD8 T cell (Total)
load('/DATA07/home/shlee/SKYRocket/SKYRocket_CD8TcellAnnotation_recent1.robj')

# CD8 T cell Post
load('/DATA07/home/shlee/SKYRocket/SKYRocket_CD8PostFinal.robj')

# CD8 T cell Trajectory 
load('/DATA07/home/shlee/SKYRocket/CD8TcellTotal_FinalTrajectoryCDS.robj')

# CD4 T cell
load('/DATA07/home/shlee/SKYRocket/SKYRocket_CD4TcellAnnotation_recent1.robj')

# CD4 T cell Trajecroty 
load('/DATA07/home/shlee/SKYRocket/CD4_TrajectoryFinal_CB_CDS.robj')
load('/DATA07/home/shlee/SKYRocket/CD4_TrajectoryFinal_NCB_CDS.robj')
load('/DATA07/home/shlee/SKYRocket/CD4_TrajectoryFinal_Post_CDS.robj')

# Myeloid
load('/DATA07/home/shlee/SKYRocket/SKYRocket_MyeloidAnnotation_recent1.robj')

# CAF
load('/DATA07/home/shlee/SKYRocket/SKYRocket_CAFAnnotation_recent.robj')
load('/DATA07/home/shlee/SKYRocket/Skyrocket_CAFMyeloid_Ref.robj')

#CellChat
load('/DATA07/home/shlee/SKYRocket/SkyrocketCellchatPost_Ref1.robj')
load('/DATA07/home/shlee/SKYRocket/SkyrocketCellchatPostMergeNCB_Ref1.robj')
load('/DATA07/home/shlee/SKYRocket/SkyrocketCellchatPostMergeCB_Ref1.robj')

####################################### Main Figure ############################################

#Figure1A

#Figure1B
use_colors <- c(
  `CD8 Tn` = "red",
  `CD8 Tcm` ="royalblue",
  `CD8 Tem` = "#39B600",
  `CD8 Trm` = "#984EA3",
  `CD8 Temra` = "#F88008",
  `CD8 Tpex` = "#FB3B9B",
  `CD8 Tex` = "#A7FFC1"
)

DimPlot(Skyroket_obj_CD8Tcell1, group.by = 'CD8T_Celltype1', cols = use_colors, label = T, repel = T)

#Figure1C
dat <- table(Skyroket_obj_CD8Tcell1$sample_id, Skyroket_obj_CD8Tcell1$CD8T_Celltype1)
dat <- prop.table(table(Skyroket_obj_CD8Tcell1$sample_id, Skyroket_obj_CD8Tcell1$CD8T_Celltype1),1)
dat <- as.data.frame(dat)
colnames(dat) <- c("sample_id", "Celltype", "Freq")

dat$PrePost <- 'none'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT12002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT39001'] <- 'Pre' 
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40002'] <- 'Post'

dat$PrePost <- factor(dat$PrePost, levels = c('Pre', 'Post'))

ggplot(dat, aes(x = Celltype, y = Freq, fill = PrePost)) + 
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA,
               position = position_dodge(width = 0.75)) + 
  scale_fill_manual(values = c("Pre" = "#3E9BCD", "Post" = "#E37A73")) +
  ylab("Relative frequency") + 
  xlab("") + 
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 9, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "top"
  ) + 
  labs(fill = "Response")

#Figure1D
Idents(Skyroket_obj_CD8Tcell1) <- 'PrePost'
markers.cluster1 <- FindMarkers(Skyroket_obj_CD8Tcell1, ident.1 = 'Post', ident.2 = 'Pre', min.pct = 0.25, logfc.threshold = 0.1, verbose = FALSE)
deg <- markers.cluster1
deg_filtered <- deg[deg$p_val_adj < 0.05, ]

up_genes <- deg_filtered[deg_filtered$avg_log2FC > 0, ]
down_genes <- deg_filtered[deg_filtered$avg_log2FC < 0, ]

up_genes <- rownames(up_genes)
down_genes <- rownames(down_genes)

up_entrez <- bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_entrez <- bitr(down_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ego_up <- enrichGO(gene = up_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   readable = TRUE)

ego_down <- enrichGO(gene = down_entrez$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pvalueCutoff = 0.05,
                     readable = TRUE)

selected_pathways_up <- c("T cell differentiation",
                          "antigen receptor-mediated signaling pathway",
                          "cellular response to radiation"
)

selected_pathways_down <- c("cytoplasmic translation",
                            "natural killer cell mediated cytotoxicity",
                            "chemotaxis",
                            "taxis",
                            
)

up_df <- ego_up@result %>%
  filter(Description %in% selected_pathways_up) %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = -log10(p.adjust),
         Group = "On_Treatment")

down_df <- ego_down@result %>%
  filter(Description %in% selected_pathways_down) %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = log10(p.adjust),  
         Group = "Baseline")

plot_df <- rbind(up_df, down_df)

ggplot(plot_df, aes(x = reorder(Description, EnrichmentScore), 
                    y = EnrichmentScore, fill = Group)) +
  geom_col() +
  coord_flip() +"Pre" = 
  scale_fill_manual(values = c("Pre" = "#3E9BCD", "Post" = "#E37A73")) +
  labs(x = NULL, y = "Enrichment Score (log10(Pvalue))", 
       title = "Baseline vs On_Treatment") +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold")
  )

#Figure1E
load('/DATA07/home/shlee/SKYRocket/SKYRocket_CD8PostFinal.robj')

radiation_gene1 <-list(c('PARP1','GATA3','YY1','NSMCE3','EP300',
                         'ZBTB1','DHX36','ATM','COPS9','CREBBP',
                         'XRCC6','NLRP1','PIK3R1','XRCC5'))

Skyroket_obj_CD8Tcell1_Post_Final <- AddModuleScore(Skyroket_obj_CD8Tcell1_Post_Final, features = radiation_gene1, name = "Radiation_Score")

Idents(Skyroket_obj_CD8Tcell1_Post_Final) <- 'CD8T_Celltype1'
DotPlot(Skyroket_obj_CD8Tcell1_Post_Final, features = "Radiation_Score1", scale.min = 0, scale.max = 15, col.min = -1.5, col.max =1.5) +
  labs(title = "Radiation Score by Response",
       x = "Module Score", y = "Pre vs Post") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(face = "bold")) +
  scale_color_gradient2(low = "blue", high = "red", midpoint = 0) 

#Figure1F
Skyroket_obj_CD8Tcell1_Pre <- subset(Skyroket_obj_CD8Tcell1, PrePost == 'Pre')
Skyroket_obj_CD8Tcell1_Pre@meta.data <- droplevels(Skyroket_obj_CD8Tcell1_Pre@meta.data)

Idents(Skyroket_obj_CD8Tcell1_Pre) <- 'CD8T_Celltype1'

DotPlot(Skyroket_obj_CD8Tcell1_Pre, features = c('TIGIT', 'PDCD1'), scale.max = 50, scale.min = 0)  

#Figure1G 
Skyroket_obj_CD8Tcell1_Tpex <- subset(Skyroket_obj_CD8Tcell1, CD8T_Celltype1 == 'CD8 Tpex')
Skyroket_obj_CD8Tcell1_Tpex@meta.data <- droplevels(Skyroket_obj_CD8Tcell1_Tpex@meta.data)

Idents(Skyroket_obj_CD8Tcell1_Tpex) <- 'PrePost'
markers.cluster1 <- FindMarkers(Skyroket_obj_CD8Tcell1_Tpex, ident.1 = "Post", ident.2 = "Pre", min.pct = 0.25, logfc.threshold = 0., verbose = FALSE)
deg <- markers.cluster1

deg$p_val_adj

deg_filtered <- deg[deg$p_val_adj < 0.05, ]

deg_filtered <- deg[deg$p_val < 0.05, ]

up_genes <- deg_filtered[deg_filtered$avg_log2FC > 0, ]
down_genes <- deg_filtered[deg_filtered$avg_log2FC < 0, ]

up_genes <- rownames(up_genes)
down_genes <- rownames(down_genes)

up_entrez <- bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_entrez <- bitr(down_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ego_up <- enrichGO(gene = up_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   readable = TRUE)

ego_down <- enrichGO(gene = down_entrez$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pvalueCutoff = 0.05,
                     readable = TRUE)


up_df <- ego_up@result %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = -log10(p.adjust),
         Group = "Post Treatment")

down_df <- ego_down@result %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = log10(p.adjust),  
         Group = "Treatment Naive")

selected_pathways_up <- c("T cell activation",
                          "T cell receptor signaling pathway",
                          "response to radiation"
)

selected_pathways_down <- c("negative regulation of apoptotic signaling pathway",
                            "negative regulation of immune response",
                            "maintenance of cell number",
)

up_df <- ego_up@result %>%
  filter(Description %in% selected_pathways_up) %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = -log10(p.adjust),
         Group = "On_Treatment")

down_df <- ego_down@result %>%
  filter(Description %in% selected_pathways_down) %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = log10(p.adjust),  
         Group = "Baseline")

plot_df <- rbind(up_df, down_df)

ggplot(plot_df, aes(x = reorder(Description, EnrichmentScore), 
                    y = EnrichmentScore, fill = Group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Baseline" = "#3E9BCD", "On_Treatment" = "#E37A73")) +
  labs(x = NULL, y = "Enrichment Score (log10(Pvalue))", 
       title = "Baseline vs On_Treatment") +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold")
  )

#Figure1H 
Skyroket_obj_CD8Tcell1_Tex <- subset(Skyroket_obj_CD8Tcell1, CD8T_Celltype1 == 'CD8 Tex')
Skyroket_obj_CD8Tcell1_Tex@meta.data <- droplevels(Skyroket_obj_CD8Tcell1_Tex@meta.data)

Idents(Skyroket_obj_CD8Tcell1_Tex) <- 'PrePost'
markers.cluster1 <- FindMarkers(Skyroket_obj_CD8Tcell1_Tex, ident.1 = "Post", ident.2 = "Pre", min.pct = 0.25, logfc.threshold = 0., verbose = FALSE)
deg <- markers.cluster1

deg_filtered <- deg[deg$p_val_adj < 0.05, ]

up_genes <- deg_filtered[deg_filtered$avg_log2FC > 0, ]
down_genes <- deg_filtered[deg_filtered$avg_log2FC < 0, ]

up_genes <- rownames(up_genes)
down_genes <- rownames(down_genes)

up_entrez <- bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_entrez <- bitr(down_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ego_up <- enrichGO(gene = up_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   readable = TRUE)

ego_down <- enrichGO(gene = down_entrez$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pvalueCutoff = 0.05,
                     readable = TRUE)

up_df <- ego_up@result %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = -log10(p.adjust),
         Group = "Post Treatment")

down_df <- ego_down@result %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = log10(p.adjust),  
         Group = "Treatment Naive")

selected_pathways_up <- c('T cell receptor signaling pathway',
                          'antigen receptor-mediated signaling pathway',
                          'immune response-activating signaling pathway',
                        )

selected_pathways_down <- c('canonical NF-kappaB signal transduction',
                            'negative regulation of leukocyte mediated cytotoxicity',
                            'negative regulation of immune effector process',
                           )

up_df <- ego_up@result %>%
  filter(Description %in% selected_pathways_up) %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = -log10(p.adjust),
         Group = "On_Treatment")

down_df <- ego_down@result %>%
  filter(Description %in% selected_pathways_down) %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = log10(p.adjust),  
         Group = "Baseline")

plot_df <- rbind(up_df, down_df)

ggplot(plot_df, aes(x = reorder(Description, EnrichmentScore), 
                    y = EnrichmentScore, fill = Group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Baseline" = "#3E9BCD", "On_Treatment" = "#E37A73")) +
  labs(x = NULL, y = "Enrichment Score (log10(Pvalue))", 
       title = "Baseline vs On_Treatment") +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold")
  )

#Figure 1I
Skyroket_obj_CD8Tcell1_CB <- subset(Skyroket_obj_CD8Tcell1, Response == 'CB')
Skyroket_obj_CD8Tcell1_CB@meta.data <- droplevels(Skyroket_obj_CD8Tcell1_CB@meta.data)

cytotoxic_genes <- c('GZMA', 'GZMB', 'GZMK', 'GNLY', 'IFNG', 'PRF1', 'NKG7')

Skyroket_obj_CD8Tcell1_CB <- AddModuleScore(Skyroket_obj_CD8Tcell1_CB, features = list(cytotoxic_genes), name = "Cytotoxic_Score")

df_CB <- Skyroket_obj_CD8Tcell1_CB@meta.data
df_CB$PrePost <- factor(df_CB$PrePost, levels = c("Pre", "Post"))

ggplot(df_CB, aes(x = PrePost, y = Cytotoxic_Score1, fill = PrePost)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = 19, outlier.size = 0.7, fill = 'white', color = 'black') +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.x = 1.35,  
    label.y = max(df_NCB$Cytotoxic_Score1, na.rm = TRUE) + 0.2  
  ) +
  labs(title = "CD8T – Cytotoxic GeneSet", y = "Score", x = "") +
  scale_fill_manual(values = c("Pre" = "#3E9BCD", "Post" = "#E37A73")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,  size = 12),
    axis.text.y = element_text(size = 12)
  )

Skyroket_obj_CD8Tcell1_NCB <- subset(Skyroket_obj_CD8Tcell1, Response == 'NCB')
Skyroket_obj_CD8Tcell1_NCB@meta.data <- droplevels(Skyroket_obj_CD8Tcell1_NCB@meta.data)

cytotoxic_genes <- c('GZMA', 'GZMB', 'GZMK', 'GNLY', 'IFNG', 'PRF1', 'NKG7')

Skyroket_obj_CD8Tcell1_NCB <- AddModuleScore(Skyroket_obj_CD8Tcell1_NCB, features = list(cytotoxic_genes), name = "Cytotoxic_Score")

df_NCB <- Skyroket_obj_CD8Tcell1_NCB@meta.data
df_NCB$PrePost <- factor(df_NCB$PrePost, levels = c("Pre", "Post"))


ggplot(df_NCB, aes(x = PrePost, y = Cytotoxic_Score1, fill = PrePost)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = 19, outlier.size = 0.7, fill = 'white', color = 'black') +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.x = 1.35,  
    label.y = max(df_NCB$Cytotoxic_Score1, na.rm = TRUE) + 0.2  
  ) +
  labs(title = "CD8T – Cytotoxic GeneSet", y = "Score", x = "") +
  scale_fill_manual(values = c("Pre" = "#3E9BCD", "Post" = "#E37A73")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,  size = 12),
    axis.text.y = element_text(size = 12)
  )

#Figure 1J
reduced_data <- reducedDims(cds)$UMAP  
cell_meta_data <- colData(cds)

new_dats <- data.frame(
  X = reduced_data[, 1],  
  Y = reduced_data[, 2],  
  T_Celltype = cell_meta_data$CD8T_Celltype1,  
  PrePost= cell_meta_data$PrePost,
  Response = cell_meta_data$Response
)

ggplot(new_dats, aes(x = X, y = Y)) + 
  geom_point(color = "grey", size = 0.5) + theme_bw(base_size = 14) + 
  geom_pointdensity(data = new_dats[new_dats$PrePost == "Pre",], adjust = 1.3) + 
  scale_color_viridis() + theme(panel.border = element_rect(color = "black", linewidth = 1),
                                panel.grid.major = element_blank(),  
                                panel.grid.minor = element_blank(),
                                axis.text = element_text(size = 15,face = "bold", color = "black"))


ggplot(new_dats, aes(x = X, y = Y)) + 
  geom_point(color = "grey", size = 0.5) + theme_bw(base_size = 14) + 
  geom_pointdensity(data = new_dats[new_dats$PrePost == "Post",], adjust = 1.2) + 
  scale_color_viridis() + theme(panel.border = element_rect(color = "black", linewidth = 1),
                                panel.grid.major = element_blank(),  
                                panel.grid.minor = element_blank(),
                                axis.text = element_text(size = 15,face = "bold", color = "black"))

#Figure 1K
sce1 <- readRDS('/DATA07/home/shlee/SKYRocket/CD8TcellPrePost_Tradeseq1s.rds')
deg <- diffEndTest(sce1)

deg <- deg[!is.na(deg$pvalue) & deg$pvalue < 0.05, ]

up_genes <- deg[deg$logFC1_2 > 0, ]
down_genes <- deg[deg$logFC1_2 < 0, ]

up_genes <- rownames(up_genes)
down_genes <- rownames(down_genes)

up_entrez <- bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_entrez <- bitr(down_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ego_up <- enrichGO(gene = up_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   readable = TRUE)

ego_down <- enrichGO(gene = down_entrez$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pvalueCutoff = 0.05,
                     readable = TRUE)

selected_pathways_up <- c('positive regulation of T cell migration',
                          'response to transforming growth factor beta',
                          'cellular response to fibroblast growth factor stimulus',
                          'innate immune response-activating signaling pathway'
                         )
selected_pathways_down <- c('regulation of type I interferon production',
                            'regulation of protein catabolic process',
                            'positive regulation of nitric oxide biosynthetic process',
                            'mitochondrial transport'
                           )

up_df <- ego_up@result %>%
  filter(Description %in% selected_pathways_up) %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = -log10(p.adjust),
         Group = "On_Treatment")

down_df <- ego_down@result %>%
  filter(Description %in% selected_pathways_down) %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = log10(p.adjust),  
         Group = "Baseline")

plot_df <- rbind(up_df, down_df)

ggplot(plot_df, aes(x = reorder(Description, EnrichmentScore), 
                    y = EnrichmentScore, fill = Group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Baseline" = "#3E9BCD", "On_Treatment" = "#E37A73")) +
  labs(x = NULL, y = "Enrichment Score (log10(Pvalue))", 
       title = "On_Treatment vs Baseline") +
  theme_minimal(base_size = 12) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

#Figure 2A 
use_colors <- c(
  `CD4 Tn` = "red",
  `CD4 Tcm` ="royalblue",
  `CD4 Tem` = "#39B600",
  `CD4 IL17` = "#984EA3",
  `CD4 CXCL13` = "#FB3B9B",
  `CD4 Treg` = "#A7FFC1"
)

DimPlot(Skyroket_obj_CD4Tcell, group.by = 'CD4T_Celltype1', cols = use_colors, label = T, repel = T)

#Figure 2B
selected_cell_type <- "CD4 CXCL13"
total_cells_per_sample <- as.data.frame(table(Skyroket_obj_CD4Tcell$sample_id))
filtered_data <- subset(Skyroket_obj_CD4Tcell, CD4T_Celltype1 == selected_cell_type)
cell_counts_per_sample <- as.data.frame(table(filtered_data$sample_id))

dat <- merge(cell_counts_per_sample, total_cells_per_sample, by = "Var1", suffixes = c("_selected", "_total"))
dat$Freq <- dat$Freq_selected / dat$Freq_total
identical(cell_counts_per_sample$Var1, total_cells_per_sample$Var1)
colnames(dat) <- c("sample_id", "Selected_Count", "Total_Count", "Freq")

dat$PrePost <- 'none'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT12002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT39001'] <- 'Pre' 
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40002'] <- 'Post'

dat$PrePost <- factor(dat$PrePost, levels = c('Pre', 'Post'))

dat$Response <- 'none'
dat$Response[dat$sample_id == 'SKYSCRHBOT04001'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT04002'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT12002'] <- 'NCB'
dat$Response[dat$sample_id == 'SKYSCRHBOT25001'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT25002'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT39001'] <- 'NCB' 
dat$Response[dat$sample_id == 'SKYSCRHBOT40001'] <- 'NCB'
dat$Response[dat$sample_id == 'SKYSCRHBOT40002'] <- 'NCB'

dat$Response <- factor(dat$Response, levels = c('NCB', 'CB'))

table(dat$Response, dat$PrePost)
unique(dat$Response)
dat
ggplot(dat, aes(x = Response, y = Freq, fill = PrePost)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Pre" = "#3E9BCD", "Post" = "#E37A73")) +
  geom_signif(data = dat, test = "wilcox.test", comparisons = list(c("Pre", "Post")), 
              map_signif_level = F, step_increase = 0.1, y_position = 0.09, color = "black") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 11, face = "bold"),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "top")

selected_cell_type <- "CD4 Treg"
total_cells_per_sample <- as.data.frame(table(Skyroket_obj_CD4Tcell$sample_id))
filtered_data <- subset(Skyroket_obj_CD4Tcell, CD4T_Celltype1 == selected_cell_type)
cell_counts_per_sample <- as.data.frame(table(filtered_data$sample_id))

dat <- merge(cell_counts_per_sample, total_cells_per_sample, by = "Var1", suffixes = c("_selected", "_total"))
dat$Freq <- dat$Freq_selected / dat$Freq_total
identical(cell_counts_per_sample$Var1, total_cells_per_sample$Var1)
colnames(dat) <- c("sample_id", "Selected_Count", "Total_Count", "Freq")

dat$PrePost <- 'none'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT12002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT39001'] <- 'Pre' 
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40002'] <- 'Post'

dat$PrePost <- factor(dat$PrePost, levels = c('Pre', 'Post'))

dat$Response <- 'none'
dat$Response[dat$sample_id == 'SKYSCRHBOT04001'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT04002'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT12002'] <- 'NCB'
dat$Response[dat$sample_id == 'SKYSCRHBOT25001'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT25002'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT39001'] <- 'NCB' 
dat$Response[dat$sample_id == 'SKYSCRHBOT40001'] <- 'NCB'
dat$Response[dat$sample_id == 'SKYSCRHBOT40002'] <- 'NCB'

dat$Response <- factor(dat$Response, levels = c('NCB', 'CB'))

table(dat$Response, dat$PrePost)
unique(dat$Response)
dat
ggplot(dat, aes(x = Response, y = Freq, fill = PrePost)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Pre" = "#3E9BCD", "Post" = "#E37A73")) +
  geom_signif(data = dat, test = "wilcox.test", comparisons = list(c("Pre", "Post")), 
              map_signif_level = F, step_increase = 0.1, y_position = 0.09, color = "black") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 11, face = "bold"),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "top")

# Figure 2C
Skyroket_obj_CD4Tcell_Pre <- subset(Skyroket_obj_CD4Tcell, PrePost == 'Pre')
Skyroket_obj_CD4Tcell_Pre@meta.data <- droplevels(Skyroket_obj_CD4Tcell_Pre@meta.data)

Idents(Skyroket_obj_CD4Tcell_Pre) <- 'CD4T_Celltype1'
DotPlot(Skyroket_obj_CD4Tcell_Pre, features = c('TIGIT', 'PDCD1'), scale.max = 80, scale.min = 0) 

# Figure 2D 
reduced_data <- reducedDims(cds_CD4CB)$UMAP  
cell_meta_data <- colData(cds_CD4CB)

new_dats <- data.frame(
  X = reduced_data[, 1],  
  Y = reduced_data[, 2],  
  T_Celltype = cell_meta_data$CD4T_Celltype1,  
  PrePost= cell_meta_data$PrePost,
  Response = cell_meta_data$Response
)

ggplot(new_dats, aes(x = X, y = Y)) + 
  geom_point(color = "grey", size = 0.5) + theme_bw(base_size = 14) + 
  geom_pointdensity(data = new_dats[new_dats$PrePost == "Pre",], adjust =4.1) + 
  scale_color_viridis() + theme(panel.border = element_rect(color = "black", linewidth = 1),
                                panel.grid.major = element_blank(),  
                                panel.grid.minor = element_blank(),
                                axis.text = element_text(size = 15,face = "bold", color = "black"))

ggplot(new_dats, aes(x = X, y = Y)) + 
  geom_point(color = "grey", size = 0.5) + theme_bw(base_size = 14) + 
  geom_pointdensity(data = new_dats[new_dats$PrePost == "Post",], adjust = 1.2) + 
  scale_color_viridis() + theme(panel.border = element_rect(color = "black", linewidth = 1),
                                panel.grid.major = element_blank(),  
                                panel.grid.minor = element_blank(),
                                axis.text = element_text(size = 15,face = "bold", color = "black"))

reduced_data <- reducedDims(cds_CD4NCB)$UMAP  
cell_meta_data <- colData(cds_CD4NCB)

new_dats <- data.frame(
  X = reduced_data[, 1],  
  Y = reduced_data[, 2],  
  T_Celltype = cell_meta_data$CD4T_Celltype1,  
  PrePost= cell_meta_data$PrePost,
  Response = cell_meta_data$Response
)

ggplot(new_dats, aes(x = X, y = Y)) + 
  geom_point(color = "grey", size = 0.5) + theme_bw(base_size = 14) + 
  geom_pointdensity(data = new_dats[new_dats$PrePost == "Pre",], adjust =2.4) + 
  scale_color_viridis() + theme(panel.border = element_rect(color = "black", linewidth = 1),
                                panel.grid.major = element_blank(),  
                                panel.grid.minor = element_blank(),
                                axis.text = element_text(size = 15,face = "bold", color = "black"))

ggplot(new_dats, aes(x = X, y = Y)) + 
  geom_point(color = "grey", size = 0.5) + theme_bw(base_size = 14) + 
  geom_pointdensity(data = new_dats[new_dats$PrePost == "Post",], adjust = 1.5) + 
  scale_color_viridis() + theme(panel.border = element_rect(color = "black", linewidth = 1),
                                panel.grid.major = element_blank(),  
                                panel.grid.minor = element_blank(),
                                axis.text = element_text(size = 15,face = "bold", color = "black"))

#Figure 2E
sce1 <- readRDS('/DATA07/home/shlee/SKYRocket/CD4TcellCBPrePost_Tradeseq1s.rds')
deg <- diffEndTest(sce1)

deg <- deg[!is.na(deg$pvalue) & deg$pvalue < 0.05, ]

up_genes <- deg[deg$logFC1_2 > 0, ]
down_genes <- deg[deg$logFC1_2 < 0, ]

up_genes <- rownames(up_genes)
down_genes <- rownames(down_genes)

up_entrez <- bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_entrez <- bitr(down_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ego_up <- enrichGO(gene = up_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   readable = TRUE)

ego_down <- enrichGO(gene = down_entrez$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pvalueCutoff = 0.05,
                     readable = TRUE)

selected_pathways_up <- c('positive regulation of NLRP3 inflammasome complex assembly',
                          'inflammasome-mediated signaling pathway',
                          'type I interferon production',
                          'response to UV',
                          'pattern recognition receptor signaling pathway',
                          )

selected_pathways_down <- c('T cell differentiation',
                            'transforming growth factor beta receptor signaling pathway',
                            'positive regulation of interleukin-10 production',
                            'positive regulation of fatty acid metabolic process',
                            'negative regulation of T cell activation'
                            )

up_df <- ego_up@result %>%
  filter(Description %in% selected_pathways_up) %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = -log10(p.adjust),
         Group = "On_Treatment")

down_df <- ego_down@result %>%
  filter(Description %in% selected_pathways_down) %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = log10(p.adjust),  
         Group = "Baseline")

plot_df <- rbind(up_df, down_df)

ggplot(plot_df, aes(x = reorder(Description, EnrichmentScore), 
                    y = EnrichmentScore, fill = Group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Baseline" = "#3E9BCD", "On_Treatment" = "#E37A73")) +
  labs(x = NULL, y = "Enrichment Score (log10(Pvalue))", 
       title = "On_Treatment vs Baseline") +
  theme_minimal(base_size = 12) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

sce1 <- readRDS('/DATA07/home/shlee/SKYRocket/CD4TcellNCBPrePost_Tradeseq1s.rds')
deg <- diffEndTest(sce1)

deg <- deg[!is.na(deg$pvalue) & deg$pvalue < 0.05, ]

up_genes <- deg[deg$logFC1_2 > 0, ]
down_genes <- deg[deg$logFC1_2 < 0, ]

up_genes <- rownames(up_genes)
down_genes <- rownames(down_genes)

up_entrez <- bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_entrez <- bitr(down_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ego_up <- enrichGO(gene = up_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   readable = TRUE)

ego_down <- enrichGO(gene = down_entrez$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pvalueCutoff = 0.05,
                     readable = TRUE)

up_df <- ego_up@result %>%
  filter(Description %in% selected_pathways_up) %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = -log10(p.adjust),
         Group = "On_Treatment")

down_df <- ego_down@result %>%
  filter(Description %in% selected_pathways_down) %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = log10(p.adjust),  
         Group = "Baseline")

selected_pathways_up <- c('regulation of substrate adhesion-dependent cell spreading',
                          'immune response-regulating cell surface receptor signaling pathway',
                          'receptor recycling',
                          'leukocyte migration',
                          'response to tumor necrosis factor'
                         )

selected_pathways_down <- c('positive regulation of leukocyte proliferation',
                            'positive regulation of T cell activation', 
                            'response to type I interferon',
                            'positive regulation of lymphocyte activation',
                            'response to interferon-beta',
                           )  

plot_df <- rbind(up_df, down_df)

ggplot(plot_df, aes(x = reorder(Description, EnrichmentScore), 
                    y = EnrichmentScore, fill = Group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Baseline" = "#3E9BCD", "On_Treatment" = "#E37A73")) +
  labs(x = NULL, y = "Enrichment Score (log10(Pvalue))", 
       title = "On_Treatment vs Baseline") +
  theme_minimal(base_size = 12) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Figure 2G 
Skyroket_obj_CD4Tcell_Treg <- subset(Skyroket_obj_CD4Tcell, CD4T_Celltype1 == 'CD4 Treg')
Skyroket_obj_CD4Tcell_Treg@meta.data <- droplevels(Skyroket_obj_CD4Tcell_Treg@meta.data)
Skyroket_obj_CD4Tcell_Treg_Post <- subset(Skyroket_obj_CD4Tcell_Treg, PrePost == 'Post')
Skyroket_obj_CD4Tcell_Treg_Post@meta.data <- droplevels(Skyroket_obj_CD4Tcell_Treg_Post@meta.data)

dents(Skyroket_obj_CD4Tcell_Treg_Post) <- 'Response'
markers.cluster1 <- FindMarkers(Skyroket_obj_CD4Tcell_Treg_Post, ident.1 = "CB", ident.2 = "NCB", min.pct = 0.1, logfc.threshold = 0.1, verbose = FALSE)
deg <- markers.cluster1

deg_filtered <- deg[deg$p_val_adj < 0.05, ]
up_genes <- deg_filtered[deg_filtered$avg_log2FC > 0, ]
down_genes <- deg_filtered[deg_filtered$avg_log2FC < 0, ]

up_genes <- rownames(up_genes)
down_genes <- rownames(down_genes)

up_entrez <- bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_entrez <- bitr(down_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ego_up <- enrichGO(gene = up_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   readable = TRUE)

ego_down <- enrichGO(gene = down_entrez$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pvalueCutoff = 0.05,
                     readable = TRUE)

up_df <- ego_up@result %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = -log10(p.adjust),
         Group = "CB")

down_df <- ego_down@result %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = log10(p.adjust),  
         Group = "NCB")

selected_pathways_up <- c('response to interferon-beta',
                          'dendritic cell chemotaxis',
                          'cellular response to interferon-beta',
                          'T cell receptor signaling pathway',
                          'chemokine-mediated signaling pathway'
)

selected_pathways_down <- c('cellular response to tumor necrosis factor',
                            'interleukin-12 production',
                            'inflammatory response to antigenic stimulus',
                            'regulation of mononuclear cell proliferation',
                            'lymphocyte proliferation',
                            'regulation of inflammatory response',
)

up_df <- ego_up@result %>%
  filter(Description %in% selected_pathways_up) %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = -log10(p.adjust),
         Group = "CB")

down_df <- ego_down@result %>%
  filter(Description %in% selected_pathways_down) %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(EnrichmentScore = log10(p.adjust),  
         Group = "NCB")

plot_df <- rbind(up_df, down_df)

ggplot(plot_df, aes(x = reorder(Description, EnrichmentScore), 
                    y = EnrichmentScore, fill = Group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("NCB" = "#3E9BCD", "CB" = "#E37A74")) +
  labs(x = NULL, y = "log10(P-value)", 
       title = "Post Treg NCB vs Post Treg CB") +
  theme_minimal(base_size = 12) +
  scale_y_continuous(limits = c(min(plot_df$EnrichmentScore) - 1, 
                                max(plot_df$EnrichmentScore) + 1)) +
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        panel.grid = element_blank()
  )

# Figure 2H
suppression_genes <- c('FOXP3', 'CTLA4', 'CCL8',  'TNFRSF1B', 'TNFRSF8', 'IL10', 'TNFRSF18', 'TNFRSF9', 'TNFRSF4')
Skyroket_obj_CD4Tcell_Treg <- AddModuleScore(Skyroket_obj_CD4Tcell_Treg, features = list(suppression_genes), name = "suppression_Score")

df_Treg <- Skyroket_obj_CD4Tcell_Treg@meta.data
df_Treg$PrePost <- factor(df_Treg$PrePost, levels = c("Pre", "Post"))

ggplot(df_Treg, aes(x = PrePost, y = suppression_Score1, fill = PrePost)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = 19, outlier.size = 0.7, fill = 'white', color = 'black') +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.x = 1.35,  
    label.y = max(df_Treg$suppression_Score1, na.rm = TRUE) + 0.4 
  ) +
  labs(title = "CD4 Treg – Suppression GeneSet", y = "Score", x = "") +
  scale_fill_manual(values = c("Pre" = "#3E9BCD", "Post" = "#E37A73")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12)
  )

suppression_genes <- c('FOXP3', 'CTLA4', 'CCL8',  'TNFRSF1B', 'TNFRSF8', 'IL10', 'TNFRSF18', 'TNFRSF9', 'TNFRSF4')
Skyroket_obj_CD4Tcell_Treg_Post <- AddModuleScore(Skyroket_obj_CD4Tcell_Treg_Post, features = list(suppression_genes), name = "suppression_Score")

df_Treg <- Skyroket_obj_CD4Tcell_Treg_Post@meta.data
df_Treg$Response <- factor(df_Treg$Response, levels = c("NCB", "CB"))

ggplot(df_Treg, aes(x = Response, y = suppression_Score1, fill = Response)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = 19, outlier.size = 0.7, fill = 'white', color = 'black') +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.x = 1.35,  
    label.y = max(df_Treg$suppression_Score1, na.rm = TRUE) + 0.4 
  ) +
  labs(title = "CD4 Treg – Suppression GeneSet", y = "Score", x = "") +
  scale_fill_manual(values = c("NCB" = "#3E9BCD", "CB" = "#E37A73")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12)
  )

#Figure 2I
Idents(Skyroket_obj_CD4Tcell_Treg) <- 'PrePost'

VlnPlot(Skyroket_obj_CD4Tcell_Treg, features = 'FOXP3',  cols = c("#3E9BCD","#E37A73"))
VlnPlot(Skyroket_obj_CD4Tcell_Treg, features = 'CTLA4',  cols = c("#3E9BCD","#E37A73"))
VlnPlot(Skyroket_obj_CD4Tcell_Treg, features = 'TNFRSF4',  cols = c("#3E9BCD","#E37A73"))
VlnPlot(Skyroket_obj_CD4Tcell_Treg, features = 'TNFRSF18',  cols = c("#3E9BCD","#E37A73"))

#Figure 3A
cds_CD4Post$pseudotime <- pseudotime(cds_CD4Post)
data.pseudo <- as.data.frame(colData(cds_CD4Post))
ggplot(data.pseudo, aes(pseudotime, 
                        reorder(CD4T_Celltype1, pseudotime, median), 
                        fill = CD4T_Celltype1)) +
                        geom_boxplot() +  
                        theme_classic()

Skyroket_obj_CD4Tcell_Post <- subset(Skyroket_obj_CD4Tcell, PrePost == 'Post')
Skyroket_obj_CD4Tcell_Post@meta.data <- droplevels(Skyroket_obj_CD4Tcell_Post@meta.data)

Transition_genes <- list(c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2", "TNFRSF18", "ENTPD1"))

Skyroket_obj_CD4Tcell_Post <- AddModuleScore(
  object = Skyroket_obj_CD4Tcell_Post,
  features = Transition_genes,
  name = c("TransitionScore")
)

pseudotime_df <- data.frame(
  pseudotime = pseudotime(cds_CD4Post),  
  TransitionScore = Skyroket_obj_CD4Tcell_Post$TransitionScore1,
  group = Skyroket_obj_CD4Tcell_Post1$Response  
)

ggplot(pseudotime_df, aes(x = pseudotime, y = TransitionScore, color = group)) +
  geom_smooth(method = "loess", se = FALSE, size = 1) +
  labs(y = "Module score", x = "Pseudotime",
       title = "Trajectory of Tem → Treg transition") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(color = "black", face = "bold", size = 13),
    axis.text.y = element_text(color = "black", face = "bold", size = 13), 
    axis.title.x = element_text(color = "black", face = "bold"),
    axis.title.y = element_text(color = "black", face = "bold") 
  )

#Figure 3B 
selectK(cellchat_NCB, pattern = "incoming")
nPatterns = 5
cellchat_NCB <- identifyCommunicationPatterns(cellchat_NCB, pattern = "incoming", k = nPatterns, width = 7, height = 13, font.size = 4)

selectK(cellchat_CB, pattern = "incoming")
nPatterns = 5
cellchat_CB <- identifyCommunicationPatterns(cellchat_CB, pattern = "incoming", k = nPatterns, width = 7, height = 13, font.size = 4)

#Figure 3C

#Figure 4D 
par(mfrow = c(1,2))
pathways.show <- c("ICAM")
netVisual_aggregate(cellchat_NCB, signaling = pathways.show, layout = "circle", vertex.weight = groupSize_NCB,
                    color.use = color_use_NCB, vertex.label.cex = 1, weight.scale = TRUE, sources.use = NULL,  targets.use = NULL,
                    vertex.receiver = vertex.receiver, arrow.size = 0.15,  thresh = 0.05, edge.width.max = 5)
netVisual_aggregate(cellchat_CB, signaling = pathways.show, layout = "circle", vertex.weight = groupSize_CB,
                    color.use = color_use_CB, vertex.label.cex = 1, weight.scale = TRUE, sources.use = NULL,  targets.use = NULL,
                    vertex.receiver = vertex.receiver, arrow.size = 0.15, thresh = 0.05, edge.width.max = 5)

pathways.show <- c("GALECTIN")
netVisual_aggregate(cellchat_NCB, signaling = pathways.show, layout = "circle", vertex.weight = groupSize_NCB,
                    color.use = color_use_NCB, vertex.label.cex = 1, weight.scale = TRUE, sources.use = NULL,  targets.use = NULL,
                    vertex.receiver = vertex.receiver, arrow.size = 0.15,  thresh = 0.05, edge.width.max = 5)
netVisual_aggregate(cellchat_CB, signaling = pathways.show, layout = "circle", vertex.weight = groupSize_CB,
                    color.use = color_use_CB, vertex.label.cex = 1, weight.scale = TRUE, sources.use = NULL,  targets.use = NULL,
                    vertex.receiver = vertex.receiver, arrow.size = 0.15, thresh = 0.05, edge.width.max = 5)

pathways.show <- c("CD86")
netVisual_aggregate(cellchat_NCB, signaling = pathways.show, layout = "circle", vertex.weight = groupSize_NCB,
                    color.use = color_use_NCB, vertex.label.cex = 1, weight.scale = TRUE, sources.use = NULL,  targets.use = NULL,
                    vertex.receiver = vertex.receiver, arrow.size = 0.15,  thresh = 0.05, edge.width.max = 5)
netVisual_aggregate(cellchat_CB, signaling = pathways.show, layout = "circle", vertex.weight = groupSize_CB,
                    color.use = color_use_CB, vertex.label.cex = 1, weight.scale = TRUE, sources.use = NULL,  targets.use = NULL,
                    vertex.receiver = vertex.receiver, arrow.size = 0.15, thresh = 0.05, edge.width.max = 5)

#Figure 3E 
SkyrocketCellchatPost_Ref1$ind <- "FALSE"
SkyrocketCellchatPost_Ref1$ind[SkyrocketCellchatPost_Ref1$CellType1s %in%c('ECM-remodeling TAMs', 'CXCL10+ Macrophage', 'cDC2', 'CD86+ Antigen-presenting CAF', 'ECM-CAFs')] <- TRUE
SkyrocketCellchatPost_Ref_L <- subset(SkyrocketCellchatPost_Ref1, subset = ind == TRUE)
SkyrocketCellchatPost_Ref_L@meta.data <- droplevels(SkyrocketCellchatPost_Ref_L@meta.data)

Idents(SkyrocketCellchatPost_Ref_L) <- 'CellType1s'
DotPlot(SkyrocketCellchatPost_Ref_L, features = c('ICAM1', 'ICAM2', 'LGALS9', 'CD86'), scale.min = 0, scale.max = 100) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#Figure 3F
SkyrocketCellchatPost_Ref1$ind <- "FALSE"
SkyrocketCellchatPost_Ref1$ind[SkyrocketCellchatPost_Ref1$CellType1s %in%c('CD8 Tpex', 'CD8 Tex', 'CD4 Tem', 'CD4 Treg', 'CD4 CXCL13')] <- TRUE
SkyrocketCellchatPost_Ref_R <- subset(SkyrocketCellchatPost_Ref1, subset = ind == TRUE)
SkyrocketCellchatPost_Ref_R@meta.data <- droplevels(SkyrocketCellchatPost_Ref_R@meta.data)

Idents(SkyrocketCellchatPost_Ref_R) <- 'CellType1s'
DotPlot(SkyrocketCellchatPost_Ref_R, features = c('ITGAL', 'CD44', 'CTLA4', 'CD28'), scale.min = 0, scale.max = 80) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#Figure 3G 

Skyroket_obj_CD4Tcell_Treg$PrePost_Response <- 'none'
Skyroket_obj_CD4Tcell_Treg$PrePost_Response[Skyroket_obj_CD4Tcell_Treg$sample_id == 'SKYSCRHBOT04001'] <- 'CB_Pre'
Skyroket_obj_CD4Tcell_Treg$PrePost_Response[Skyroket_obj_CD4Tcell_Treg$sample_id == 'SKYSCRHBOT04002'] <- 'CB_Post'
Skyroket_obj_CD4Tcell_Treg$PrePost_Response[Skyroket_obj_CD4Tcell_Treg$sample_id == 'SKYSCRHBOT12002'] <- 'NCB_Post'
Skyroket_obj_CD4Tcell_Treg$PrePost_Response[Skyroket_obj_CD4Tcell_Treg$sample_id == 'SKYSCRHBOT25001'] <- 'CB_Pre'
Skyroket_obj_CD4Tcell_Treg$PrePost_Response[Skyroket_obj_CD4Tcell_Treg$sample_id == 'SKYSCRHBOT25002'] <- 'CB_Post'
Skyroket_obj_CD4Tcell_Treg$PrePost_Response[Skyroket_obj_CD4Tcell_Treg$sample_id == 'SKYSCRHBOT39001'] <- 'NCB_Pre' 
Skyroket_obj_CD4Tcell_Treg$PrePost_Response[Skyroket_obj_CD4Tcell_Treg$sample_id == 'SKYSCRHBOT40001'] <- 'NCB_Pre'
Skyroket_obj_CD4Tcell_Treg$PrePost_Response[Skyroket_obj_CD4Tcell_Treg$sample_id == 'SKYSCRHBOT40002'] <- 'NCB_Post'

Skyroket_obj_CD4Tcell_Treg$PrePost_Response <- factor(Skyroket_obj_CD4Tcell_Treg$PrePost_Response, levels = c('CB_Pre', 'CB_Post', 'NCB_Pre', 'NCB_Post'))

genes_all1 <- c("TIGIT","CTLA4","PDCD1","ICOS","TNFRSF1B","TNFRSF18","TNFRSF4","TNFRSF9", "CCR4", 'ENTPD1')

unique(Skyroket_obj_CD4Tcell_Treg$PrePost_Response)
Idents(Skyroket_obj_CD4Tcell_Treg) <- "PrePost_Response"
expression_lv = as.data.frame(AverageExpression(Skyroket_obj_CD4Tcell_Treg), rownames = NULL)
expression_lv = data.frame(gene = rownames(expression_lv), 'CB_Pre' = expression_lv$RNA.CB.Pre,
                           'CB_Post' = expression_lv$RNA.CB.Post, 'NCB_Pre' = expression_lv$RNA.NCB.Pre,
                           'CD4 NCB_Post' = expression_lv$RNA.NCB.Post)

chemo_mac = subset(expression_lv, subset = expression_lv$gene %in% genes_all1)
head(chemo_mac)

row.names(chemo_mac) = chemo_mac$gene
chemo_mac[, 1] = NULL

# Scale the matrix from -1 to 1
scale_to_minus1_1 <- function(x) {
  return((2 * (x - min(x))) / (max(x) - min(x)) - 1)
}

z_chemo_mac = t(apply(chemo_mac, 1, scale_to_minus1_1))
new_mat = z_chemo_mac

colnames(new_mat) <- c('CB_Pre', 'CB_Post', 'NCB_Pre', 'NCB_Post')
library(pheatmap)

new_mat_ordered <- new_mat[genes_all1, ]

pheatmap(t(new_mat_ordered), main = 'TNK DEG', 
         cluster_rows = F, cluster_cols = F, border_color = 'Black', angle_col = 90)

################################### Supplementary Figure #######################################

#Supplemenraty Figure 1A 
use_colors <- c(
  `Epithelial` = "red",
  `NK_T` ="royalblue",
  `Fibroblast` = "#39B600",
  `Myeloid` = "#984EA3",
  `B_Plasma` = "#F88008",
  `SMC` = "yellow",
  `Endothelial` = "#9D5A2F",
  `Mast` = "#F781BF")

genes <- c('EPCAM', 'KRT7', 'KRT8',     # Epithelial
           'PLVAP', 'VWF', 'MYCT1',       # Endothelial
           'DCN', 'LUM', 'COL1A1',    # Fibroblast
           'RERGL', 'CDH6', 'MYH11',    # SMC
           'CD3D', 'CD3E', 'GNLY',      # NK_T
           'IGKV1-39', 'IGHV4-34', 'IGHV1-24',     # B_Plasma
           'LYZ', 'CD14', 'C1QA',    # Myeloid
           'CTSG', 'TPSAB1', 'KIT' # Mast
          )

SKY_PrePost_Meta_sub$New <- factor(SKY_PrePost_Meta_sub$New, levels = c('Epithelial', 'Endothelial', 'Fibroblast', 'SMC', 'NK_T', 'B_Plasma', 'Myeloid', 'Mast'))

Idents(SKY_PrePost_Meta_sub) <- 'New'
DotPlot(SKY_PrePost_Meta_sub, features = genes, cols = c("blue", "red")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_gradient(low = "blue", high = "red") + 
  coord_flip()

#Supplemenraty Figure 1B
DimPlot(SKY_PrePost_Meta_sub, group.by = 'PrePost', cols = c(Pre = "grey", Post = "#9D5A1F"))

#Supplemenraty Figure 1C
DimPlot(SKY_PrePost_Meta_sub, group.by = 'Response', cols = c(NCB = "#4977B0", CB = "#C21F25"))

#Supplemenraty Figure 1D
FeaturePlot(SKY_PrePost_Meta_sub, features = c('CD3D', 'NKG7'), max.cutoff = 1.8) # NK-T
FeaturePlot(SKY_PrePost_Meta_sub, features = c('LYZ', 'CD14'), max.cutoff = 1.8) # Myeloid
FeaturePlot(SKY_PrePost_Meta_sub, features = c('COL10A1', 'SFRP2'), max.cutoff = 1.8) # Fibroblast
FeaturePlot(SKY_PrePost_Meta_sub, features = c('HDC', 'RHEX'), max.cutoff = 1.8) # Mast cells
FeaturePlot(SKY_PrePost_Meta_sub, features = c('CD19', 'FCRLA'), max.cutoff = 1.8) # B cells
FeaturePlot(SKY_PrePost_Meta_sub, features = c('EPCAM', 'KRT6A'), max.cutoff = 1.8) # Epithelial
FeaturePlot(SKY_PrePost_Meta_sub, features = c('ACKR1', 'EMCN'), max.cutoff = 1.8) # Endothelial
FeaturePlot(SKY_PrePost_Meta_sub, features = c('FHL5', 'RGS5'), max.cutoff = 1.8) # SMC

#Supplemenraty Figure 2A 
gg_dat <- SKY_PrePost_Meta_sub@meta.data

stage_order <- c("Pre", "Post")
gg_dat <- gg_dat %>%
  arrange(match(PrePost, stage_order))

gg_dat$sample_id1 <- as.character(gg_dat$sample_id1)
gg_dat$sample_id1[gg_dat$sample_id1 == 'PA01_Pre'] <- 'PA01_Baseline'
gg_dat$sample_id1[gg_dat$sample_id1 == 'PA03_Pre'] <- 'PA03_Baseline'
gg_dat$sample_id1[gg_dat$sample_id1 == 'PA04_Pre'] <- 'PA04_Baseline'
gg_dat$sample_id1[gg_dat$sample_id1 == 'PA05_Pre'] <- 'PA05_Baseline'
gg_dat$sample_id1[gg_dat$sample_id1 == 'PA01_Post'] <- 'PA01_Ontreatment'
gg_dat$sample_id1[gg_dat$sample_id1 == 'PA02_Post'] <- 'PA02_Ontreatment'
gg_dat$sample_id1[gg_dat$sample_id1 == 'PA03_Post'] <- 'PA03_Ontreatment'
gg_dat$sample_id1[gg_dat$sample_id1 == 'PA05_Post'] <- 'PA05_Ontreatment'
  
gg_dat$sample_id1 <- factor(gg_dat$sample_id1, levels = c('PA01_Baseline', 'PA03_Baseline', 'PA04_Baseline', 'PA05_Baseline', 'PA01_Ontreatment', 'PA02_Ontreatment', 'PA03_Ontreatment', 'PA05_Ontreatment'))
gg_dat$New <- factor(gg_dat$New, levels = c('Epithelial', 'NK_T', 'Fibroblast', 'Myeloid', 'B_Plasma', 'SMC', 'Endothelial', 'Mast'))

ggplot(gg_dat, aes(x = sample_id1, fill = New)) + 
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = use_colors) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, face = "bold", 
                               color = "black", size = 12),
    axis.text.y = element_text(face = "bold", 
                               color = "black", size = 12)
  )

ggplot(gg_dat, aes(x = PrePost, fill = New)) + 
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = use_colors) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, face = "bold", 
                               color = "black", size = 12),
    axis.text.y = element_text(face = "bold", 
                               color = "black", size = 12)
  )

#Supplemenraty Figure 2B 
dat <- table(SKY_PrePost_Meta_sub$sample_id, SKY_PrePost_Meta_sub$New)
dat <- prop.table(table(SKY_PrePost_Meta_sub$sample_id, SKY_PrePost_Meta_sub$New),1)
dat <- as.data.frame(dat)
colnames(dat) <- c("sample_id", "Celltype", "Freq")

dat$PrePost <- 'none'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT12002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT39001'] <- 'Pre' 
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40002'] <- 'Post'

dat$PrePost <- factor(dat$PrePost, levels = c('Pre', 'Post'))

ggplot(dat, aes(x = Celltype, y = Freq, fill = PrePost)) + 
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA, position = position_dodge(width = 0.75)) + 
  ylab("Relative frequency") + 
  xlab("") + 
  theme_bw() +
  scale_fill_manual(values = c("Pre" = "#3E9BCD", "Post" = "#E37A73")) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 9, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "top"
  ) + 
  labs(fill = "Response")

#Supplemenraty Figure 2C
TNKmarkers1_CD8S <- c('CCR7', 'SELL', 'TCF7', 'LEF1', 'NELL2', 'IL7R', 'ANXA1', 'FOS', 'JUN', 'ZNF683', 'TOB1',
                      'GZMK', 'GPR183', 'EOMES', 'ITM2C', 'GZMA', 'GZMB', 'GZMH', 'CX3CR1', 'FGFBP2', 'SPON2',
                      'PRF1', 'NKG7', 'IFNG', 'CCL4', 'CCL5' ,'CCL4L2', 'TIGIT', 'CTLA4', 'PDCD1', 'HAVCR2', 'CXCL13')

Idents(Skyroket_obj_CD8Tcell1) <- "CD8T_Celltype1"
expression_lv = as.data.frame(AverageExpression(Skyroket_obj_CD8Tcell1), rownames = NULL)
expression_lv = data.frame(gene = rownames(expression_lv), 'CD8 Tn' = expression_lv$RNA.CD8.Tn,
                           'CD8 Tcm' = expression_lv$RNA.CD8.Tcm, 'CD8 Tem' = expression_lv$RNA.CD8.Tem,
                           'CD8 Trm' = expression_lv$RNA.CD8.Trm, 'CD8 Temra' = expression_lv$RNA.CD8.Temra,
                           'CD8 Tpex' = expression_lv$RNA.CD8.Tpex, 'CD8 Tex' = expression_lv$RNA.CD8.Tex
)

chemo_mac = subset(expression_lv, subset = expression_lv$gene %in% TNKmarkers1_CD8S)

row.names(chemo_mac) = chemo_mac$gene
chemo_mac[, 1] = NULL

scale_to_minus1_1 <- function(x) {
  return((2 * (x - min(x))) / (max(x) - min(x)) - 1)
}

z_chemo_mac = apply(chemo_mac, 1, scale)
new_mat = t(z_chemo_mac)

z_chemo_mac = t(apply(chemo_mac, 1, scale_to_minus1_1))
new_mat = z_chemo_mac

colnames(new_mat) <- c('CD8 Tn', 'CD8 Tcm', 'CD8 Tem', 'CD8 Trm', 'CD8 Temra', 'CD8 Tpex', 'CD8 Tex')

TNKmarkers1_CD8S <- rev(TNKmarkers1_CD8S)
new_mat_ordered <- new_mat[TNKmarkers1_CD8S, ]

pheatmap(t(new_mat_ordered), main = 'TNK DEG', 
         cluster_rows = F, cluster_cols = F, border_color = 'Black', angle_col = 90) + coord_flip()

#Supplemenraty Figure 2D
Skyroket_obj_CD8Tcell1$ind <- "FALSE"
Skyroket_obj_CD8Tcell1$ind[Skyroket_obj_CD8Tcell1$CD8T_Celltype1 %in% c('CD8 Tpex', 'CD8 Tex')] <- TRUE
Skyroket_obj_CD8Tcell_Tex <- subset(Skyroket_obj_CD8Tcell1, subset = ind == TRUE)
Skyroket_obj_CD8Tcell_Tex@meta.data <- droplevels(Skyroket_obj_CD8Tcell_Tex@meta.data)

Idents(Skyroket_obj_CD8Tcell_Tex) <- "CD8T_Celltype1"
expression_lv = as.data.frame(AverageExpression(Skyroket_obj_CD8Tcell_Tex), rownames = NULL)
expression_lv = data.frame(gene = rownames(expression_lv), 'CD8 Tpex' = expression_lv$RNA.CD8.Tpex, 'CD8 Tex' = expression_lv$RNA.CD8.Tex
)

TNKmarkers1_Tpex <- c('TCF7', 'CXCL13')
chemo_mac = subset(expression_lv, subset = expression_lv$gene %in% TNKmarkers1_Tpex)
head(chemo_mac)

row.names(chemo_mac) = chemo_mac$gene
chemo_mac[, 1] = NULL

# Scale the matrix from -1 to 1
scale_to_minus1_1 <- function(x) {
  return((2 * (x - min(x))) / (max(x) - min(x)) - 1)
}

z_chemo_mac = apply(chemo_mac, 1, scale)
new_mat = t(z_chemo_mac)

z_chemo_mac = t(apply(chemo_mac, 1, scale_to_minus1_1))
new_mat = z_chemo_mac

colnames(new_mat) <- c('CD8 Tpex', 'CD8 Tex')

new_mat_ordered <- new_mat[TNKmarkers1_Tpex, ]

pheatmap(new_mat_ordered, main = 'TNK DEG', 
         cluster_rows = F, cluster_cols = F, border_color = 'Black', angle_col = 90)

#Supplemenraty Figure 2E
Skyroket_obj_CD8Tcell1_Post_Final_Tpex <- subset(Skyroket_obj_CD8Tcell1_Post_Final, CD8T_Celltype1 == 'CD8 Tpex')
Skyroket_obj_CD8Tcell1_Post_Final_Tpex@meta.data <- droplevels(Skyroket_obj_CD8Tcell1_Post_Final_Tpex@meta.data)

radiation_gene <-list(c('PARP1','GATA3','YY1','NSMCE3','EP300',
                         'ZBTB1','DHX36','ATM','COPS9','CREBBP',
                         'XRCC6','NLRP1','PIK3R1','XRCC5'))

Skyroket_obj_CD8Tcell1_Post_Final_Tpex <- AddModuleScore(Skyroket_obj_CD8Tcell1_Post_Final_Tpex, features = radiation_gene, name = "Radiation_Score")

Idents(Skyroket_obj_CD8Tcell1_Post_Final_Tpex)<- 'Response'

DotPlot(Skyroket_obj_CD8Tcell1_Post_Final_Tpex, features = "Radiation_Score1", scale.max = 25, scale.min = 0) +
  labs(title = "Radiation Score by Response",
       x = "Module Score", y = "Pre vs Post") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(face = "bold")) +
  scale_color_gradient2(low = "blue", high = "red", midpoint = 0) 

#Supplemenraty Figure 2F
Skyroket_obj_CD8Tcell1_Tpex <- subset(Skyroket_obj_CD8Tcell1, CD8T_Celltype1 == 'CD8 Tpex')
Skyroket_obj_CD8Tcell1_Tpex@meta.data <- droplevels(Skyroket_obj_CD8Tcell1_Tpex@meta.data)

Idents(Skyroket_obj_CD8Tcell1_Tpex) <- 'PrePost'
markers.cluster1 <- FindMarkers(Skyroket_obj_CD8Tcell1_Tpex, ident.1 = "Post", ident.2 = "Pre", min.pct = 0.1, logfc.threshold = 0.1, verbose = FALSE)
data <- markers.cluster1
EnhancedVolcano(data,
                lab = rownames(data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pointSize = 0.3,
                labSize = 3.5,
                drawConnectors = TRUE,
                widthConnectors = 0.1,
                pCutoff = 0.05,
                FCcutoff = 0.25,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = F,                
                xlim = c(-5, 5),
                ylim = c(0, 13))

#Supplemenraty Figure 2G
Skyroket_obj_CD8Tcell1_Tex <- subset(Skyroket_obj_CD8Tcell1, CD8T_Celltype1 == 'CD8 Tex')
Skyroket_obj_CD8Tcell1_Tex@meta.data <- droplevels(Skyroket_obj_CD8Tcell1_Tex@meta.data)

Idents(Skyroket_obj_CD8Tcell1_Tex) <- 'PrePost'
markers.cluster1 <- FindMarkers(Skyroket_obj_CD8Tcell1_Tex, ident.1 = "Post", ident.2 = "Pre", min.pct = 0.1, logfc.threshold = 0.1, verbose = FALSE)
data <- markers.cluster1
EnhancedVolcano(data,
                lab = rownames(data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pointSize = 0.3,
                labSize = 3.5,
                drawConnectors = TRUE,
                widthConnectors = 0.1,
                pCutoff = 0.05,
                FCcutoff = 0.25,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = F,
                xlim = c(-3.5, 3.5),
                ylim = c(0, 45))
Skyroket_obj_CD8Tcell1

#Supplemenraty Figure 2H 
selected_cell_type <- 'CD8 Tpex'
total_cells_per_sample <- as.data.frame(table(Skyroket_obj_CD8Tcell1$sample_id))
filtered_data <- subset(Skyroket_obj_CD8Tcell1, CD8T_Celltype1 == selected_cell_type)
cell_counts_per_sample <- as.data.frame(table(filtered_data$sample_id))
dat <- merge(cell_counts_per_sample, total_cells_per_sample, by = "Var1", suffixes = c("_selected", "_total"))
dat$Freq <- dat$Freq_selected / dat$Freq_total
identical(cell_counts_per_sample$Var1, total_cells_per_sample$Var1)

colnames(dat) <- c("sample_id", "Selected_Count", "Total_Count", "Freq")

dat$PrePost <- 'none'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT12002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT39001'] <- 'Pre' 
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40002'] <- 'Post'

dat$PrePost <- factor(dat$PrePost, levels = c('Pre', 'Post'))

dat$Response <- 'none'
dat$Response[dat$sample_id == 'SKYSCRHBOT04001'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT04002'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT12002'] <- 'NCB'
dat$Response[dat$sample_id == 'SKYSCRHBOT25001'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT25002'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT39001'] <- 'NCB' 
dat$Response[dat$sample_id == 'SKYSCRHBOT40001'] <- 'NCB'
dat$Response[dat$sample_id == 'SKYSCRHBOT40002'] <- 'NCB'

dat$Response <- factor(dat$Response, levels = c('NCB', 'CB'))

ggplot(dat, aes(x = Response, y = Freq, fill = PrePost)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Pre" = "#3E9BCD",
                               "Post" = "#E37A73")) +   
  geom_signif(data = dat, test = "wilcox.test", comparisons = list(c("Pre", "Post")), 
              map_signif_level = F, step_increase = 0.1, y_position = 0.09, color = "black") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 9, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "top") 

selected_cell_type <- 'CD8 Tex'
total_cells_per_sample <- as.data.frame(table(Skyroket_obj_CD8Tcell1$sample_id))
filtered_data <- subset(Skyroket_obj_CD8Tcell1, CD8T_Celltype1 == selected_cell_type)
cell_counts_per_sample <- as.data.frame(table(filtered_data$sample_id))
dat <- merge(cell_counts_per_sample, total_cells_per_sample, by = "Var1", suffixes = c("_selected", "_total"))
dat$Freq <- dat$Freq_selected / dat$Freq_total
identical(cell_counts_per_sample$Var1, total_cells_per_sample$Var1)

colnames(dat) <- c("sample_id", "Selected_Count", "Total_Count", "Freq")

dat$PrePost <- 'none'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT12002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT39001'] <- 'Pre' 
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40002'] <- 'Post'

dat$PrePost <- factor(dat$PrePost, levels = c('Pre', 'Post'))

dat$Response <- 'none'
dat$Response[dat$sample_id == 'SKYSCRHBOT04001'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT04002'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT12002'] <- 'NCB'
dat$Response[dat$sample_id == 'SKYSCRHBOT25001'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT25002'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT39001'] <- 'NCB' 
dat$Response[dat$sample_id == 'SKYSCRHBOT40001'] <- 'NCB'
dat$Response[dat$sample_id == 'SKYSCRHBOT40002'] <- 'NCB'

dat$Response <- factor(dat$Response, levels = c('NCB', 'CB'))

ggplot(dat, aes(x = Response, y = Freq, fill = PrePost)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Pre" = "#3E9BCD",
                               "Post" = "#E37A73")) +   
  geom_signif(data = dat, test = "wilcox.test", comparisons = list(c("Pre", "Post")), 
              map_signif_level = F, step_increase = 0.1, y_position = 0.09, color = "black") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 9, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "top") 

#Supplemenraty Figure 3A 
dat <- table(Skyroket_obj_CD4Tcell$sample_id, Skyroket_obj_CD4Tcell$CD4T_Celltype1)
dat <- prop.table(table(Skyroket_obj_CD4Tcell$sample_id, Skyroket_obj_CD4Tcell$CD4T_Celltype1),1)
dat <- as.data.frame(dat)
colnames(dat) <- c("sample_id", "Celltype", "Freq")

dat$PrePost <- 'none'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT12002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT39001'] <- 'Pre' 
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40002'] <- 'Post'

dat$PrePost <- factor(dat$PrePost, levels = c('Pre', 'Post'))

ggplot(dat, aes(x = Celltype, y = Freq, fill = PrePost)) + 
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA, position = position_dodge(width = 0.75)) + 
  scale_fill_manual(values = c("Pre" = "#3E9BCD", "Post" = "#E37A73")) +
  ylab("Relative frequency") + 
  xlab("") + 
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 9, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "top"
  ) + 
  labs(fill = "Response")

#Supplemenraty Figure 3B
TNKmarkers1_CD4 <- c('SELL', 'TCF7', 'CCR7', 'LEF1', 'FHIT', 'NOSIP', 'IL7R', 'CD40LG', 'ANXA1','FOS', 'FOSB',
                     'JUN', 'GPR183', 'ZNF683', 'GZMK', 'GZMA', 'GZMH', 'CXCR3', 'XCL1', 'XCL2', 'STAT1', 
                     'STAT4', 'IFNG', 'PRF1', 'NKG7', 'CXCL13', 'CD200', 'IL21', 'CXCR5', 'BCL6', 'NMB', 
                     'ICA1', 'GNG4', 'EBI3', 'ITM2A', 'IL2RA', 'FOXP3', 'TNFRSF4', 'TNFRSF9',
                     'TNFRSF18', 'IKZF2', 'IL1R2', 'IL1R1', 'LAIR2', 'HAVCR2', 'PDCD1','CTLA4', 'TIGIT', 'IL17A', 'IL23R', 'IL22')

Idents(Skyroket_obj_CD4Tcell) <- "CD4T_Celltype1"
expression_lv = as.data.frame(AverageExpression(Skyroket_obj_CD4Tcell), rownames = NULL)
expression_lv = data.frame(gene = rownames(expression_lv), 'CD4 Tn' = expression_lv$RNA.CD4.Tn,
                           'CD4 Tcm' = expression_lv$RNA.CD4.Tcm, 'CD4 Tem' = expression_lv$RNA.CD4.Tem,
                           'CD4 IL17' = expression_lv$RNA.CD4.IL17, 'CD4 CXCL13' = expression_lv$RNA.CD4.CXCL13,
                           'CD4 Treg' = expression_lv$RNA.CD4.Treg
)

chemo_mac = subset(expression_lv, subset = expression_lv$gene %in% TNKmarkers1_CD4)
head(chemo_mac)

row.names(chemo_mac) = chemo_mac$gene
chemo_mac[, 1] = NULL

# Scale the matrix from -1 to 1
scale_to_minus1_1 <- function(x) {
  return((2 * (x - min(x))) / (max(x) - min(x)) - 1)
}

z_chemo_mac = apply(chemo_mac, 1, scale)
new_mat = t(z_chemo_mac)

z_chemo_mac = t(apply(chemo_mac, 1, scale_to_minus1_1))
new_mat = z_chemo_mac

colnames(new_mat) <- c('CD4 Tn', 'CD4 Tcm', 'CD4 Tem', 'CD4 IL17', 'CD4 CXCL13', 'CD4 Treg')

TNKmarkers1_CD4 <- rev(TNKmarkers1_CD4)
new_mat_ordered <- new_mat[TNKmarkers1_CD4, ]

pheatmap(t(new_mat_ordered), main = 'TNK DEG', 
         cluster_rows = F, cluster_cols = F, border_color = 'Black', angle_col = 90)

#Supplemenraty Figure 4A
use_colors <- c(
  `Monocyte` = "#E76BF3",
  `CXCL10+ Macrophage` ="royalblue",
  `SPP1+ Macrophage` = "#984EA3",
  `FOLR2+ Macrophage` = "#39B600",
  `TREM2+ Macrophage` = "red",
  `Proliferating Macrophage` = "yellow",
  `ECM-remodeling TAMs` = "#FB3B9B",
  `cDC1` = "#A7FFC1" ,
  `cDC2` =  "#F88008" ,
  `pDC` = "#9D5A2F"
)

DimPlot(Skyroket_obj_Myeloid1s, group.by = 'M_Celltype', cols = use_colors, label = T, repel = T)

#Supplemenraty Figure 4B
Myeloidmarker <- c('FCN1', 'VCAN', 'S100A9', 'S100A12', 'FCGR3A', 'CDKN1C', 'MTSS1', 
                   'CCR2','CLEC10A', 'CX3CR1', 'CSF1R', 'FCGBP', 'PLD4', 'FCER1A',
                   'CXCL10', 'CXCL9', 'CXCL11', 'GBP5', 'STAT1', 'IDO1', 'SPP1',
                   'MMP9', 'MMP19', 'EMP1', 'TREM2','CCL18', 'GPNMB', 'LGALS3', 'CSTB', 'APOC1',
                   'FABP5', 'FOLR2', 'STAB1', 'F13A1', 'SELENOP','MRC1', 'SLC40A1', 'RNASE1', 'MAF', 'LYVE1', 
                   'VSIG4', 'CD163', 'TIMD4', 'CLEC4F', 'SIGLEC1','CD5L', 'TLR4',
                   'FTL', 'MKI67', 'TOP2A', 'COL1A1', 'COL1A2', 'COL3A1', 'XCR1', 'CLEC9A', 'CD1C', 'CLEC10A')

Idents(Skyroket_obj_Myeloid1s) <- "M_Celltype"
expression_lv = as.data.frame(AverageExpression(Skyroket_obj_Myeloid1s), rownames = NULL)
expression_lv = data.frame(gene = rownames(expression_lv), 'Monocyte' = expression_lv$RNA.Monocyte,
                           'CXCL10+ Macrophage' = expression_lv$RNA.CXCL10..Macrophage, 'SPP1+ Macrophage' = expression_lv$RNA.SPP1..Macrophage,
                           'FOLR2+ Macrophage' = expression_lv$RNA.FOLR2..Macrophage, 'TREM2+ Macrophage' = expression_lv$RNA.TREM2..Macrophage,
                           'Proliferating Macrophage' = expression_lv$RNA.Proliferating.Macrophage, 'ECM.remodeling.TAMs' = expression_lv$RNA.ECM.remodeling.TAMs,
                           'cDC1' = expression_lv$RNA.cDC1, 'cDC2' = expression_lv$RNA.cDC2)

chemo_mac = subset(expression_lv, subset = expression_lv$gene %in% Myeloidmarker)
head(chemo_mac)

row.names(chemo_mac) = chemo_mac$gene
chemo_mac[, 1] = NULL

# Scale the matrix from -1 to 1
scale_to_minus1_1 <- function(x) {
  return((2 * (x - min(x))) / (max(x) - min(x)) - 1)
}

z_chemo_mac = t(apply(chemo_mac, 1, scale_to_minus1_1))
new_mat = z_chemo_mac

colnames(new_mat) <- c('Monocyte', 'CXCL10+ Macrophage', 'SPP1+ Macrophage', 'FOLR2+ Macrophage', 'TREM2+ Macrophage',
                       'Proliferating Macrophage', 'ECM-remodeling TAMs', 'cDC1', 'cDC2')

new_mat_ordered <- new_mat[Myeloidmarker, ]

pheatmap(new_mat_ordered, main = 'TNK DEG', 
         cluster_rows = F, cluster_cols = F, border_color = 'Black', angle_col = 90)

#Supplemenraty Figure 4C 
dat <- table(Skyroket_obj_Myeloid1s$sample_id, Skyroket_obj_Myeloid1s$M_Celltype)
dat <- prop.table(table(Skyroket_obj_Myeloid1s$sample_id, Skyroket_obj_Myeloid1s$M_Celltype),1)
dat <- as.data.frame(dat)
colnames(dat) <- c("sample_id", "Celltype", "Freq")

dat$PrePost <- 'none'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT12002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT39001'] <- 'Pre' 
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40002'] <- 'Post'

dat$PrePost <- factor(dat$PrePost, levels = c('Pre', 'Post'))

ggplot(dat, aes(x = Celltype, y = Freq, fill = PrePost)) + 
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA, position = position_dodge(width = 0.75)) + 
  scale_fill_manual(values = c("Pre" = "#3E9BCD", "Post" = "#E37A73")) +
  ylab("Relative frequency") + 
  xlab("") + 
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 11, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "top"
  ) + 
  labs(fill = "Response")

#Supplemenraty Figure 4D
Idents(Skyroket_obj_Myeloid1s) <- 'PrePost'
markers.cluster1 <- FindMarkers(Skyroket_obj_Myeloid1s, ident.1 = "Post", ident.2 = "Pre", min.pct = 0.25, logfc.threshold = 0.15, verbose = FALSE)
data <- markers.cluster1
EnhancedVolcano(data,
                lab = rownames(data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pointSize = 0.3,
                labSize = 2.5,
                drawConnectors = TRUE,
                widthConnectors = 0.1,
                pCutoff = 0.05,
                FCcutoff = 0.25,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = F)

#Supplemenraty Figure 4E
Skyroket_obj_Myeloid1s1 <- subset(Skyroket_obj_Myeloid1s, M_Celltype != 'Proliferating Macrophage')
Skyroket_obj_Myeloid1s1@meta.data <- droplevels(Skyroket_obj_Myeloid1s1@meta.data)

Idents(Skyroket_obj_Myeloid1s1) <- 'PrePost'
markers.cluster1 <- FindMarkers(Skyroket_obj_Myeloid1s1, ident.1 = "Post", ident.2 = "Pre", min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
deg <- markers.cluster1

deg_filtered <- deg[deg$p_val_adj < 0.05, ]

up_genes <- deg_filtered[deg_filtered$avg_log2FC > 0, ]
down_genes <- deg_filtered[deg_filtered$avg_log2FC < 0, ]

up_genes <- rownames(up_genes)
down_genes <- rownames(down_genes)

up_entrez <- bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_entrez <- bitr(down_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ego_up <- enrichGO(gene = up_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   readable = TRUE)

ego_down <- enrichGO(gene = down_entrez$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pvalueCutoff = 0.05,
                     readable = TRUE)

selected_pathways_up <- c('antigen processing and presentation of exogenous peptide antigen',
                          'cellular response to epidermal growth factor stimulus',
                          'MHC class II protein complex assembly',
                          'response to radiation'
                         )

selected_pathways_down <- c('positive regulation of canonical NF-kappaB signal transduction',
                            'regulation of innate immune response',
                            'canonical NF-kappaB signal transduction',
                            'leukocyte cell-cell adhesion'
                           )

up_df <- ego_up@result %>%
  filter(Description %in% selected_pathways_up) %>%
  arrange(p.adjust) %>%
  mutate(EnrichmentScore = -log10(p.adjust),
         Group = "On-Treatment")

down_df <- ego_down@result %>%
  filter(Description %in% selected_pathways_down) %>%
  arrange(p.adjust) %>%
  mutate(EnrichmentScore = log10(p.adjust),  
         Group = "Baseline")

plot_df <- rbind(up_df, down_df)

ggplot(plot_df, aes(x = reorder(Description, EnrichmentScore), 
                    y = EnrichmentScore, fill = Group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Baseline" = "#3E9BCD", "On-Treatment" = "#E37A73")) +
  labs(x = NULL, y = "log10(P-value)", 
       title = "On-Treatment vs Baseline") +
  theme_minimal(base_size = 12) +
  scale_y_continuous(limits = c(min(plot_df$EnrichmentScore) - 7, 
                                max(plot_df$EnrichmentScore) + 7)) +
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        panel.grid = element_blank()
  )

#Supplemenraty Figure 4F
Skyroket_obj_Myeloid1s1_Post <- subset(Skyroket_obj_Myeloid1s1, PrePost == 'Post')
Skyroket_obj_Myeloid1s1_Post@meta.data <- droplevels(Skyroket_obj_Myeloid1s1_Post@meta.data)

radiation_gene <-list(c('COL3A1', 'MT-ND3', 'FOS', 'PARP1', 'IKBIP', 'NABP2', 'BABAM1', 'COL6A3', 'PPP1CC', 'PER1', 'NSMCE3', 'COL6A2', 'INIP', 'RRM1', 'AKT1'))

Skyroket_obj_Myeloid1s1_Post <- AddModuleScore(Skyroket_obj_Myeloid1s1_Post, features = radiation_gene, name = "Radiation_Score")

Idents(Skyroket_obj_Myeloid1s1_Post) <- 'M_Celltype'
DotPlot(Skyroket_obj_Myeloid1s1_Post, features = "Radiation_Score1", scale.max = 75, scale.min = 0) +
  labs(title = "Radiation Score by PrePost",
       x = "Module Score", y = "Pre vs Post") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(face = "bold")) +
  scale_color_gradient2(low = "blue", high = "red", midpoint = 0)

Idents(Skyroket_obj_Myeloid1s1_Post) <- 'Response'
DotPlot(Skyroket_obj_Myeloid1s1_Post, features = "Radiation_Score1", scale.max = 75, scale.min = 0) +
  labs(title = "DNA Regulatory Score by PrePost",
       x = "Module Score", y = "Pre vs Post") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(face = "bold")) +
  scale_color_gradient2(low = "blue", high = "red", midpoint = 0)

#Supplemenraty Figure 4G
selected_cell_type <- 'ECM-remodeling TAMs'
total_cells_per_sample <- as.data.frame(table(Skyroket_obj_Myeloid1s$sample_id))
filtered_data <- subset(Skyroket_obj_Myeloid1s, M_Celltype == selected_cell_type)
cell_counts_per_sample <- as.data.frame(table(filtered_data$sample_id))
dat <- merge(cell_counts_per_sample, total_cells_per_sample, by = "Var1", suffixes = c("_selected", "_total"))
dat$Freq <- dat$Freq_selected / dat$Freq_total
identical(cell_counts_per_sample$Var1, total_cells_per_sample$Var1)

colnames(dat) <- c("sample_id", "Selected_Count", "Total_Count", "Freq")

dat$PrePost <- 'none'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT12002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT39001'] <- 'Pre' 
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40002'] <- 'Post'

dat$PrePost <- factor(dat$PrePost, levels = c('Pre', 'Post'))

dat$Response <- 'none'
dat$Response[dat$sample_id == 'SKYSCRHBOT04001'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT04002'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT12002'] <- 'NCB'
dat$Response[dat$sample_id == 'SKYSCRHBOT25001'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT25002'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT39001'] <- 'NCB' 
dat$Response[dat$sample_id == 'SKYSCRHBOT40001'] <- 'NCB'
dat$Response[dat$sample_id == 'SKYSCRHBOT40002'] <- 'NCB'

dat$Response <- factor(dat$Response, levels = c('NCB', 'CB'))

table(dat$Response, dat$PrePost)
unique(dat$Response)
dat
ggplot(dat, aes(x = Response, y = Freq, fill = PrePost)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Pre" = "#3E9BCD", "Post" = "#E37A73")) +
  geom_signif(data = dat, test = "wilcox.test", comparisons = list(c("Pre", "Post")), 
              map_signif_level = F, step_increase = 0.1, y_position = 0.09, color = "black") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size =13),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "top")

#Supplemenraty Figure 5A
use_colors <- c(
  `CD86+ Antigen-presenting CAF` = "#FB3B9B",
  `CD86- Antigen-presenting CAF` = "red",
  `Myofibroblatic CAF` ="royalblue",
  `Inflammatory CAF` = "#984EA3",
  `Developmental CAF` = "#39B600",
  `ECM-CAFs` = "#E76BF3"
)

DimPlot(Skyroket_obj_CAF1s1, group.by = 'C_Celltype1s', cols = use_colors, label = T, repel = T)

#Supplemenraty Figure 5B
CAFmarker1 <- c('CCND1', 'RCAN1', 'FST', 
                'SCX', 'CRTAC1','MIA',   
                'CST1', 'CCDC102B','HEY1',
                'STMN1', 'PTTG1', 'MKI67',
                'CFD', 'PLA2G2A', 'DACT2',
                'EPYC', 'ACTG2', 'SYT1') 

Skyroket_obj_CAF1s1$C_Celltype1s <- factor(Skyroket_obj_CAF1s1$C_Celltype1s, levels = c('ECM-CAFs','CD86- Antigen-presenting CAF', 'CD86+ Antigen-presenting CAF',
                                                                                      'Developmental CAF', 'Inflammatory CAF', 'Myofibroblatic CAF'))
Idents(Skyroket_obj_CAF1s1) <- "C_Celltype1s"
expression_lv = as.data.frame(AverageExpression(Skyroket_obj_CAF1s1), rownames = NULL)
expression_lv = data.frame(gene = rownames(expression_lv), 'ECM-CAFs' = expression_lv$RNA.ECM.CAFs,
                           'CD86- Antigen-presenting CAF' = expression_lv$RNA.CD86..Antigen.presenting.CAF,
                           'CD86+ Antigen-presenting CAF' = expression_lv$RNA.CD86..Antigen.presenting.CAF.1,
                           'Developmental CAF' = expression_lv$RNA.Developmental.CAF,
                           'Inflammatory CAF' = expression_lv$RNA.Inflammatory.CAF,
                           'Myofibroblatic CAF' = expression_lv$RNA.Myofibroblatic.CAF)

chemo_mac = subset(expression_lv, subset = expression_lv$gene %in% CAFmarker1)

row.names(chemo_mac) = chemo_mac$gene
chemo_mac[, 1] = NULL

scale_to_minus1_1 <- function(x) {
  return((2 * (x - min(x))) / (max(x) - min(x)) - 1)
}

z_chemo_mac = t(apply(chemo_mac, 1, scale_to_minus1_1))
new_mat = z_chemo_mac

colnames(new_mat) <- c('ECM-CAFs', 'CD86- Antigen-presenting CAF', 'CD86+ Antigen-presenting CAF',  'Developmental CAF', 'Inflammatory CAF', 'Myofibroblatic CAF')

CAFmarker1 <- rev(CAFmarker1)
new_mat_ordered <- new_mat[CAFmarker1, ]

pheatmap(new_mat_ordered, main = 'TNK DEG', 
         cluster_rows = F, cluster_cols = F, border_color = 'Black', angle_col = 90)

#Supplemenraty Figure 5C
dat <- table(Skyroket_obj_CAF1s1$sample_id, Skyroket_obj_CAF1s1$C_Celltype1s)
dat <- prop.table(table(Skyroket_obj_CAF1s1$sample_id, Skyroket_obj_CAF1s1$C_Celltype1s),1)
dat <- as.data.frame(dat)
colnames(dat) <- c("sample_id", "Celltype", "Freq")

dat$PrePost <- 'none'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT12002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40002'] <- 'Post'

dat$PrePost <- factor(dat$PrePost, levels = c('Pre', 'Post'))

ggplot(dat, aes(x = Celltype, y = Freq, fill = PrePost)) + 
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA, position = position_dodge(width = 0.75)) + 
  scale_fill_manual(values = c("Pre" = "#3E9BCD", "Post" = "#E37A73")) +
  ylab("Relative frequency") + 
  xlab("") + 
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 11, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = 13, face = "bold"),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "top"
  ) + 
  labs(fill = "Prepost")

#Supplemenraty Figure 5D
Idents(Skyroket_obj_CAF1s1) <- 'PrePost'
markers.cluster1 <- FindMarkers(Skyroket_obj_CAF1s1, ident.1 = "Post", ident.2 = "Pre", min.pct = 0.25, logfc.threshold = 0.1, verbose = FALSE)
data <- markers.cluster1
EnhancedVolcano(data,
                lab = rownames(data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pointSize = 0.5,
                labSize = 2.5,
                drawConnectors = TRUE,
                widthConnectors = 0.1,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = F)

#Supplemenraty Figure 5E
Idents(Skyroket_obj_CAF1s1) <- 'PrePost'
markers.cluster1 <- FindMarkers(Skyroket_obj_CAF1s1, ident.1 = "Post", ident.2 = "Pre", min.pct = 0.25, logfc.threshold = 0.1, verbose = FALSE)
deg <- markers.cluster1

deg_filtered <- deg[deg$p_val_adj < 0.05, ]

up_genes <- deg_filtered[deg_filtered$avg_log2FC > 0, ]
down_genes <- deg_filtered[deg_filtered$avg_log2FC < 0, ]

up_genes <- rownames(up_genes)
down_genes <- rownames(down_genes)

up_entrez <- bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_entrez <- bitr(down_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ego_up <- enrichGO(gene = up_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   readable = TRUE)

ego_down <- enrichGO(gene = down_entrez$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pvalueCutoff = 0.05,
                     readable = TRUE)

up_df <- ego_up@result %>%
  filter(Description %in% selected_pathways_up) %>%
  arrange(p.adjust) %>%
  mutate(EnrichmentScore = -log10(p.adjust),
         Group = "On-Treatment")

down_df <- ego_down@result %>%
  filter(Description %in% selected_pathways_down) %>%
  arrange(p.adjust) %>%
  mutate(EnrichmentScore = log10(p.adjust),  
         Group = "Baseline")

plot_df <- rbind(up_df, down_df)

ggplot(plot_df, aes(x = reorder(Description, EnrichmentScore), 
                    y = EnrichmentScore, fill = Group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Baseline" = "#3E9BCD", "On-Treatment" = "#E37A74")) +
  labs(x = NULL, y = "log10(P-value)", 
       title = "On-Treatment vs Baseline") +
  theme_minimal(base_size = 12) +
  scale_y_continuous(limits = c(min(plot_df$EnrichmentScore) - 7, 
                                max(plot_df$EnrichmentScore) + 7)) +
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        panel.grid = element_blank()
  )

selected_pathways_up <- c('mitochondrion organization',
                          'collagen metabolic process',
                          'extracellular matrix organization',
                          'response to UV',
                          )

selected_pathways_down <- c('response to transforming growth factor beta',
                            'response to interferon-alpha',
                            'response to hypoxia',
                            'canonical NF-kappaB signal transduction'
                           )
                           
#Supplemenraty Figure 5F
Skyroket_obj_CAF1s1_Post <- subset(Skyroket_obj_CAF1s1, PrePost == 'Post')
Skyroket_obj_CAF1s1_Post@meta.data <- droplevels(Skyroket_obj_CAF1s1_Post@meta.data)

radiation_gene1 <-list(c('MMP2','CDKN1A','DDB2','BAX','BBC3','CRYAB','RUVBL2','MDM2','PIK3R1','PCNA','N4BP1',
                         'CREBBP','HRAS','TP53INP1','TSPYL5','ATR','GRB2','RAD1','AKT2','NMT2','SWI5','HUS1',
                         'AQP1','CCND2','TAF1','NLRP1','ERCC1','RHBDD1','SMPD1','ERCC4','EIF2AK4'))

Skyroket_obj_CAF1s1_Post <- AddModuleScore(Skyroket_obj_CAF1s1_Post, features = radiation_gene1, name = "Radiation_Score")

Skyroket_obj_CAF1s1_Post$C_Celltype1s <- factor(Skyroket_obj_CAF1s1_Post$C_Celltype1s, levels = c('CD86- Antigen-presenting CAF', 'CD86+ Antigen-presenting CAF', 'ECM-CAFs',
                                                                                                  'Myofibroblatic CAF', 'Developmental CAF', 'Inflammatory CAF'))
Idents(Skyroket_obj_CAF1s1_Post)<- 'C_Celltype1s'
DotPlot(Skyroket_obj_CAF1s1_Post, features = "Radiation_Score1", scale.max = 80, scale.min = 0) +
  labs(title = "Radiation Score by PrePost",
       x = "Module Score", y = "Pre vs Post") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(face = "bold")) +
  scale_color_gradient2(low = "blue", high = "red", midpoint = 0)

Idents(Skyroket_obj_CAF1s1_Post)<- 'Response'
DotPlot(Skyroket_obj_CAF1s1_Post, features = "Radiation_Score1", scale.max = 80, scale.min = 0) +
  labs(title = "DNA Regulatory Score by PrePost",
       x = "Module Score", y = "Pre vs Post") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(face = "bold")) +
  scale_color_gradient2(low = "blue", high = "red", midpoint = 0)

#Supplemenraty Figure 5G
selected_cell_type <- 'ECM-CAFs'
total_cells_per_sample <- as.data.frame(table(Skyroket_obj_CAF1s1$sample_id))
filtered_data <- subset(Skyroket_obj_CAF1s1, C_Celltype1s == selected_cell_type)
cell_counts_per_sample <- as.data.frame(table(filtered_data$sample_id))
dat <- merge(cell_counts_per_sample, total_cells_per_sample, by = "Var1", suffixes = c("_selected", "_total"))
dat$Freq <- dat$Freq_selected / dat$Freq_total
identical(cell_counts_per_sample$Var1, total_cells_per_sample$Var1)

colnames(dat) <- c("sample_id", "Selected_Count", "Total_Count", "Freq")

dat$PrePost <- 'none'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT12002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40002'] <- 'Post'

dat$PrePost <- factor(dat$PrePost, levels = c('Pre', 'Post'))

dat$Response <- 'none'
dat$Response[dat$sample_id == 'SKYSCRHBOT04001'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT04002'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT12002'] <- 'NCB'
dat$Response[dat$sample_id == 'SKYSCRHBOT25001'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT25002'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT40001'] <- 'NCB'
dat$Response[dat$sample_id == 'SKYSCRHBOT40002'] <- 'NCB'

dat$Response <- factor(dat$Response, levels = c('NCB', 'CB'))

ggplot(dat, aes(x = Response, y = Freq, fill = PrePost)) +
  geom_boxplot() +
  geom_signif(data = dat, test = "wilcox.test", comparisons = list(c("Pre", "Post")), 
              map_signif_level = F, step_increase = 0.1, y_position = 0.09, color = "black") +
  scale_fill_manual(values = c("Pre" = "#3E9BCD", "Post" = "#E37A73")) +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 13),
    axis.text.y = element_text(color = "black", size = 13, face = "bold"),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "top")

#Supplemenraty Figure 6A
Skyrocket_CAFMyeloid_RefPre <- subset(Skyrocket_CAFMyeloid_Ref, PrePost == 'Pre')
Skyrocket_CAFMyeloid_RefPre@meta.data <- droplevels(Skyrocket_CAFMyeloid_RefPre@meta.data)

Skyrocket_CAFMyeloid_RefPre$CellType1 <- factor(Skyrocket_CAFMyeloid_RefPre$CellType1, levels = c('Inflammatory CAF', 'CD86- Antigen-presenting CAF', 'CD86+ Antigen-presenting CAF',
                                                                                                  'Developmental CAF', 'ECM-CAFs', 'Myofibroblatic CAF', 'CXCL10+ Macrophage',
                                                                                                  'Monocyte', 'TREM2+ Macrophage', 'cDC2', 'FOLR2+ Macrophage', 'SPP1+ Macrophage',
                                                                                                  'ECM-remodeling TAMs', 'cDC1'))
Idents(Skyrocket_CAFMyeloid_RefPre) <- 'CellType1'
DotPlot(Skyrocket_CAFMyeloid_RefPre, features = c('CD274', 'PVR', 'NECTIN2')) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#Supplemenraty Figure 6B
selected_cell_type <- 'CXCL10+ Macrophage'
total_cells_per_sample <- as.data.frame(table(Skyroket_obj_Myeloid1s$sample_id))
filtered_data <- subset(Skyroket_obj_Myeloid1s, M_Celltype == selected_cell_type)

cell_counts_per_sample <- as.data.frame(table(filtered_data$sample_id))
dat <- merge(cell_counts_per_sample, total_cells_per_sample, by = "Var1", suffixes = c("_selected", "_total"))
dat$Freq <- dat$Freq_selected / dat$Freq_total

identical(cell_counts_per_sample$Var1, total_cells_per_sample$Var1)

colnames(dat) <- c("sample_id", "Selected_Count", "Total_Count", "Freq")

dat$PrePost <- 'none'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT12002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT39001'] <- 'Pre' 
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40002'] <- 'Post'

dat$PrePost <- factor(dat$PrePost, levels = c('Pre', 'Post'))

dat$Response <- 'none'
dat$Response[dat$sample_id == 'SKYSCRHBOT04001'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT04002'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT12002'] <- 'NCB'
dat$Response[dat$sample_id == 'SKYSCRHBOT25001'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT25002'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT39001'] <- 'NCB' 
dat$Response[dat$sample_id == 'SKYSCRHBOT40001'] <- 'NCB'
dat$Response[dat$sample_id == 'SKYSCRHBOT40002'] <- 'NCB'

dat$Response <- factor(dat$Response, levels = c('NCB', 'CB'))

table(dat$Response, dat$PrePost)
unique(dat$Response)
dat
ggplot(dat, aes(x = Response, y = Freq, fill = PrePost)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Pre" = "#3E9BCD", "Post" = "#E37A73")) +
  geom_signif(data = dat, test = "wilcox.test", comparisons = list(c("Pre", "Post")), 
              map_signif_level = F, step_increase = 0.1, y_position = 0.09, color = "black") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size =13),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "top")

selected_cell_type <- 'cDC2'
total_cells_per_sample <- as.data.frame(table(Skyroket_obj_Myeloid1s$sample_id))
filtered_data <- subset(Skyroket_obj_Myeloid1s, M_Celltype == selected_cell_type)

cell_counts_per_sample <- as.data.frame(table(filtered_data$sample_id))
dat <- merge(cell_counts_per_sample, total_cells_per_sample, by = "Var1", suffixes = c("_selected", "_total"))
dat$Freq <- dat$Freq_selected / dat$Freq_total

identical(cell_counts_per_sample$Var1, total_cells_per_sample$Var1)

colnames(dat) <- c("sample_id", "Selected_Count", "Total_Count", "Freq")

dat$PrePost <- 'none'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT12002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT39001'] <- 'Pre' 
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40002'] <- 'Post'

dat$PrePost <- factor(dat$PrePost, levels = c('Pre', 'Post'))

dat$Response <- 'none'
dat$Response[dat$sample_id == 'SKYSCRHBOT04001'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT04002'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT12002'] <- 'NCB'
dat$Response[dat$sample_id == 'SKYSCRHBOT25001'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT25002'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT39001'] <- 'NCB' 
dat$Response[dat$sample_id == 'SKYSCRHBOT40001'] <- 'NCB'
dat$Response[dat$sample_id == 'SKYSCRHBOT40002'] <- 'NCB'

dat$Response <- factor(dat$Response, levels = c('NCB', 'CB'))

table(dat$Response, dat$PrePost)
unique(dat$Response)
dat
ggplot(dat, aes(x = Response, y = Freq, fill = PrePost)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Pre" = "#3E9BCD", "Post" = "#E37A73")) +
  geom_signif(data = dat, test = "wilcox.test", comparisons = list(c("Pre", "Post")), 
              map_signif_level = F, step_increase = 0.1, y_position = 0.09, color = "black") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size =13),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "top")

selected_cell_type <- 'CD86+ Antigen-presenting CAF'
total_cells_per_sample <- as.data.frame(table(Skyroket_obj_CAF1s1$sample_id))
filtered_data <- subset(Skyroket_obj_CAF1s1, C_Celltype1s == selected_cell_type)
cell_counts_per_sample <- as.data.frame(table(filtered_data$sample_id))
dat <- merge(cell_counts_per_sample, total_cells_per_sample, by = "Var1", suffixes = c("_selected", "_total"))
dat$Freq <- dat$Freq_selected / dat$Freq_total
identical(cell_counts_per_sample$Var1, total_cells_per_sample$Var1)

colnames(dat) <- c("sample_id", "Selected_Count", "Total_Count", "Freq")

dat$PrePost <- 'none'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT04002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT12002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT25002'] <- 'Post'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40001'] <- 'Pre'
dat$PrePost[dat$sample_id == 'SKYSCRHBOT40002'] <- 'Post'

dat$PrePost <- factor(dat$PrePost, levels = c('Pre', 'Post'))

dat$Response <- 'none'
dat$Response[dat$sample_id == 'SKYSCRHBOT04001'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT04002'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT12002'] <- 'NCB'
dat$Response[dat$sample_id == 'SKYSCRHBOT25001'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT25002'] <- 'CB'
dat$Response[dat$sample_id == 'SKYSCRHBOT40001'] <- 'NCB'
dat$Response[dat$sample_id == 'SKYSCRHBOT40002'] <- 'NCB'

dat$Response <- factor(dat$Response, levels = c('NCB', 'CB'))

ggplot(dat, aes(x = Response, y = Freq, fill = PrePost)) +
  geom_boxplot() +
  geom_signif(data = dat, test = "wilcox.test", comparisons = list(c("Pre", "Post")), 
              map_signif_level = F, step_increase = 0.1, y_position = 0.09, color = "black") +
  scale_fill_manual(values = c("Pre" = "#3E9BCD", "Post" = "#E37A73")) +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 13),
    axis.text.y = element_text(color = "black", size = 13, face = "bold"),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "top")

#Supplemenraty Figure 6C
par(mfrow = c(1,2))
pathways.show <- c("NECTIN")
netVisual_aggregate(cellchat_NCB, signaling = pathways.show, layout = "circle", vertex.weight = groupSize_NCB,
                    color.use = color_use_NCB, vertex.label.cex = 1, weight.scale = TRUE, sources.use = NULL,  targets.use = NULL,
                    vertex.receiver = vertex.receiver, arrow.size = 0.15,  thresh = 0.05, edge.width.max = 5)
netVisual_aggregate(cellchat_CB, signaling = pathways.show, layout = "circle", vertex.weight = groupSize_CB,
                    color.use = color_use_CB, vertex.label.cex = 1, weight.scale = TRUE, sources.use = NULL,  targets.use = NULL,
                    vertex.receiver = vertex.receiver, arrow.size = 0.15, thresh = 0.05, edge.width.max = 5)

#Supplemenraty Figure 6D
pathway_name <- "NECTIN"

df_CB <- subsetCommunication(cellchat_CB, signaling = pathway_name) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(group = "CB")

df_NCB <- subsetCommunication(cellchat_NCB, signaling = pathway_name) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(group = "NCB")

df_combined <- dplyr::bind_rows(df_CB, df_NCB)

info_flow_summary <- df_combined %>%
  dplyr::group_by(interaction_name, group) %>%
  dplyr::summarise(sum_prob = sum(prob, na.rm = TRUE), .groups = "drop")

info_flow_wide <- info_flow_summary %>%
  tidyr::pivot_wider(
    names_from  = group,
    values_from = sum_prob,
    values_fill = 0
  )

info_flow_wide <- info_flow_wide %>%
  mutate(
    total = CB + NCB,
    CB_ratio = ifelse(total > 0, CB / total, 0),
    NCB_ratio = ifelse(total > 0, NCB / total, 0)
  )

info_flow_plot <- info_flow_wide %>%
  select(interaction_name, CB_ratio, NCB_ratio) %>%
  tidyr::pivot_longer(cols = c(CB_ratio, NCB_ratio),
                      names_to = "group", values_to = "ratio")

info_flow_plot$group <- factor(info_flow_plot$group,
                               levels = c("NCB_ratio", "CB_ratio"),
                               labels = c("NCB", "CB"))

ggplot(info_flow_plot, aes(x = reorder(interaction_name, -ratio), y = ratio, fill = group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste0(pathway_name, " signaling - Normalized info flow"),
       y = "Relative contribution (0–1 scale)",
       x = "Ligand -> Receptor") +
  scale_fill_manual(values = c("NCB" = "#2b3f8c", "CB" = "#e31a1c")) +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()
  )

df_for_test <- df_combined %>%
  mutate(
    LR_pair = case_when(
      "interaction_name" %in% names(.) ~ as.character(interaction_name_2),
      TRUE ~ as.character(interaction_name)  # fallback
    ),
    group = as.character(group)
  ) %>%
  filter(!is.na(LR_pair), group %in% c("CB", "NCB"))

split_list <- split(df_for_test, df_for_test$LR_pair)

pval_df <- do.call(rbind, lapply(names(split_list), function(lr){
  d <- split_list[[lr]]
  cb  <- d$prob[d$group == "CB"]
  ncb <- d$prob[d$group == "NCB"]
  
  data.frame(
    LR_pair = lr,
    n_CB  = length(cb),
    n_NCB = length(ncb),
    p_value = if (length(cb) >= 1 && length(ncb) >= 1) wilcox.test(cb, ncb)$p.value else NA_real_
  )
}))

pval_df

SkyrocketCellchatPost_Ref1$CellType1s <- SkyrocketCellchatPost_Ref1$CellType1
SkyrocketCellchatPost_Ref1$CellType1s <- as.character(SkyrocketCellchatPost_Ref1$CellType1s)
SkyrocketCellchatPost_Ref1$CellType1s[SkyrocketCellchatPost_Ref1$CellType1s == 'CAF'] <- 'ECM-CAFs'
unique(SkyrocketCellchatPost_Ref1$CellType1s)

#Supplemenraty Figure 6E
SkyrocketCellchatPost_Ref1$ind <- "FALSE"
SkyrocketCellchatPost_Ref1$ind[SkyrocketCellchatPost_Ref1$CellType1s %in%c('ECM-remodeling TAMs', 'CXCL10+ Macrophage', 'cDC2', 'CD86+ Antigen-presenting CAF', 'ECM-CAFs')] <- TRUE
SkyrocketCellchatPost_Ref_L <- subset(SkyrocketCellchatPost_Ref1, subset = ind == TRUE)
SkyrocketCellchatPost_Ref_L@meta.data <- droplevels(SkyrocketCellchatPost_Ref_L@meta.data)

SkyrocketCellchatPost_Ref1$ind <- "FALSE"
SkyrocketCellchatPost_Ref1$ind[SkyrocketCellchatPost_Ref1$CellType1s %in%c('CD8 Tpex', 'CD8 Tex', 'CD4 Tem', 'CD4 Treg', 'CD4 CXCL13')] <- TRUE
SkyrocketCellchatPost_Ref_R <- subset(SkyrocketCellchatPost_Ref1, subset = ind == TRUE)
SkyrocketCellchatPost_Ref_R@meta.data <- droplevels(SkyrocketCellchatPost_Ref_R@meta.data)

SkyrocketCellchatPost_Ref_L$CellType1s <- factor(SkyrocketCellchatPost_Ref_L$CellType1s, levels = c('CXCL10+ Macrophage', 'ECM-remodeling TAMs', 'cDC2', 'CD86+ Antigen-presenting CAF', 'ECM-CAFs'))
Idents(SkyrocketCellchatPost_Ref_L) <- 'CellType1s'
DotPlot(SkyrocketCellchatPost_Ref_L, features = c('NECTIN2', 'NECTIN3')) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

SkyrocketCellchatPost_Ref_R$CellType1s <- factor(SkyrocketCellchatPost_Ref_R$CellType1s, levels = c('CD8 Tpex', 'CD8 Tex', 'CD4 Tem', 'CD4 CXCL13', 'CD4 Treg'))
Idents(SkyrocketCellchatPost_Ref_R) <- 'CellType1s'
DotPlot(SkyrocketCellchatPost_Ref_R, features = c('TIGIT', 'CD226', 'PVR')) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
