# Bassem work individually in scRNA seq
library(Seurat)
library(ggplot2)
library(dplyr)
sc_data = Read10X(data.dir = "../scRNA GitHub Data/GSM5573504_sample39.csv")
sc_data2 = Read10X(data.dir = "../scRNA GitHub Data/GSM5573499_sample34.csv/")

seurat_object = (CreateSeuratObject(counts = sc_data, project = "sample39"))
seurat_object = (CreateSeuratObject(counts = sc_data, project = "sample34"))

# 2: QC
sample_d1 = CreateSeuratObject(counts = Read10X("../scRNA GitHub Data/GSM5573504_sample39.csv") ,
                               project = "SMPL_Donor1")
sample_d2 = CreateSeuratObject(counts = Read10X("../scRNA GitHub Data/GSM5573499_sample34.csv/") ,
                               project = "SMPL_Donor2")

sample_d1 [["percent.mt"]] = PercentageFeatureSet(sample_d1, pattern = "^MT-")
sample_d2 [["percent.mt"]] = PercentageFeatureSet(sample_d2, pattern = "^MT-")

# visualization before QC by #Violin Plot
VlnPlot(sample_d1,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
VlnPlot(sample_d2,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)

# Scatter Plot -> appear coorelation between 2 factors
FeatureScatter(sample_d1, feature1= "nCount_RNA", feature2= "nFeature_RNA")
FeatureScatter(sample_d1, feature1= "nCount_RNA", feature2= "percent.mt")
FeatureScatter(sample_d2, feature1= "nCount_RNA", feature2= "nFeature_RNA")
FeatureScatter(sample_d2, feature1= "nCount_RNA", feature2= "percent.mt")
# 1 Cell Filtering
# subset-> هي اللي بتعمل filtration
sample_d1 = subset(sample_d1, subset = nFeature_RNA >200&nFeature_RNA<6000&percent.mt<10)
sample_d2 = subset(sample_d2, subset = nFeature_RNA >200&nFeature_RNA<6000&percent.mt<10)


#بعد عمليه الfiltration دي هتلاقي الخلايا اتلفترت و التاكيد انك شوفها فيال violin plots

# 2 Doublet-> using: scDblFinder
library(scDblFinder)
library(SingleCellExperiment)
# A.Convert Seurat to SingleCellExperiment
sce_d1 = as.SingleCellExperiment(sample_d1)
sce_d2 = as.SingleCellExperiment(sample_d2)

# B.Run scDblFinder
set.seed(100)
sce_d1 = scDblFinder(sce_d1) 
sce_d2 = scDblFinder(sce_d2)        


# C.Seurat,so add result back to Seurat
sample_d1$doublet_score = colData(sce_d1)$scDblFinder.score
sample_d1$doublet_class = colData(sce_d1)$scDblFinder.class
table(sample_d1$doublet_class)
sample_d2$doublet_score = colData(sce_d2)$scDblFinder.score
sample_d2$doublet_class = colData(sce_d2)$scDblFinder.class
table(sample_d2$doublet_class)
#هيظهر تحت ان فيه 304 doublet المفروض يتشالوا
library(ggplot2)
ggplot(sample_d1@meta.data, aes(x= nCount_RNA, y= doublet_score, colour = doublet_class)) +
  geom_point(size=1, alpha = 0.6) + labs(x= "nCount_RNA", y="Doublet Score") + theme_classic()

VlnPlot(sample_d1, features = "doublet_score", group.by = "doublet_class")

ggplot(sample_d2@meta.data, aes(x= nCount_RNA, y= doublet_score, colour = doublet_class)) +
  geom_point(size=1, alpha = 0.6) + labs(x= "nCount_RNA", y="Doublet Score") + theme_classic()

VlnPlot(sample_d2, features = "doublet_score", group.by = "doublet_class")
# D. Doublet Removal

sample_d1 = subset(sample_d1, subset = doublet_class == "singlet")
sample_d2 = subset(sample_d2, subset = doublet_class == "singlet")


#3.Correct Ambient RNA contamination
# Best tools are in Python and bash like: CellBender
so will explain code, but not run in R
cellbender remove-background \
--input raw_feature_bc_matrix.h5 \
--output cellbender_output.h5 \
--expected-cells 5000 \
--total-droplets-include 15000
# 3: Normalization-Scaling-FeauterSelection
# A. Normalization:
sample_d1 = NormalizeData(sample_d1, normalization.method = "LogNormalize", scale.factor = 10000)
sample_d2 = NormalizeData(sample_d2, normalization.method = "LogNormalize", scale.factor = 10000)

# B. Scaling:
sample_d1 = ScaleData(sample_d1)
sample_d2 = ScaleData(sample_d2)

# C. Find Variable Genes
sample_d1 = FindVariableFeatures(sample_d1, selection.method = "vst", nfeatures = 2000)
sample_d2 = FindVariableFeatures(sample_d2, selection.method = "vst", nfeatures = 2000)

# Visaulization Highly Variable Genes
VariableFeaturePlot(sample_d1)
head(VariableFeatures(sample_d1), 10)
LabelPoints(plot = VariableFeaturePlot(sample_d1), points = head(VariableFeatures(sample_d1),10),
            repel = TRUE)

VariableFeaturePlot(sample_d2)
head(VariableFeatures(sample_d2), 10)
LabelPoints(plot = VariableFeaturePlot(sample_d2), points = head(VariableFeatures(sample_d2),10),
            repel = TRUE)
#Scaling on first 2000 gene only 
sample_d1 = ScaleData(sample_d1, features = VariableFeatures(sample_d1))
sample_d2 = ScaleData(sample_d2, features = VariableFeatures(sample_d2))

# 4: Dimentionally Reduction (PCA -Clustering- UMAP)
# 1.Run PCA
sample_d1 = RunPCA(sample_d1, features = VariableFeatures(object = sample_d1))
sample_d2 = RunPCA(sample_d2, features = VariableFeatures(object = sample_d2))

# 2.Visualize PCA
DimPlot(sample_d1, reduction = "pca")

DimPlot(sample_d2, reduction = "pca")
# 3. Elbow Plot
ElbowPlot(sample_d1)
ElbowPlot(sample_d2)
# 4. PCA Heatmap
DimHeatmap(sample_d1, dims = 1:10, cells = 500, balanced = TRUE)
DimHeatmap(sample_d2, dims = 1:10, cells = 500, balanced = TRUE)

# when you visualize pca heatmap plot will show yellow color and purple-> this indicate
#to YELLOW-> highly expressed, PURPLE-> Low Expressed

# Find Neibhor Graph 
sample_d1 = FindNeighbors(sample_d1, dims = 1:10)
sample_d2 = FindNeighbors(sample_d2, dims = 1:10)

# Clustering the cells
sample_d1 = FindClusters(sample_d1, resolution = 0.5)
sample_d2 = FindClusters(sample_d2, resolution = 0.5)

# RUn UMAP
sample_d1 = RunUMAP(sample_d1, dims = 1:10)
sample_d2 = RunUMAP(sample_d2, dims = 1:10)

# Visualize UMAP
DimPlot(sample_d1, reduction = "umap", label = TRUE) + ggtitle("SMPL_Donor 1 _UMAP")
DimPlot(sample_d2, reduction = "umap", label = TRUE) + ggtitle("SMPL_Donor 2 _UMAP")


# 5: (Cell Type Annotation)

#==========================================================
# B) Automated Cell-Type Annotation Using SingleR (Ref Based)

# # Install required packages if needed to work in this method 
library(SingleR)
library(celldex)
# celldex include (have) load human primary cell Atlas References
hpca_ref = celldex::HumanPrimaryCellAtlasData()
# Convert Seurat to SingleCellExperiment
SCE_d1 <- as.SingleCellExperiment(sample_d1)
SCE_d2 <- as.SingleCellExperiment(sample_d2)

# Run SingleR
singler_d1 <- SingleR(
  test = SCE_d1,
  ref = hpca_ref,
  labels = hpca_ref$label.main
)

singler_d2 <- SingleR(
  test = SCE_d2,
  ref = hpca_ref,
  labels = hpca_ref$label.main
)
# Add SingleR Results to Seurat
sample_d1$SingleR_annotation <- singler_d1$labels
sample_d2$SingleR_annotation <- singler_d2$labels

# Visualize SingleR UMAP
p1=DimPlot(
  sample_d1,
  reduction = "umap",
  group.by = "SingleR_annotation",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("SMPL_Donor1 — SingleR Automated Annotation")

p2=DimPlot(
  sample_d2,
  reduction = "umap",
  group.by = "SingleR_annotation",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("SMPL_Donor2 — SingleR Automated Annotation")

p1+p2
###############################################################################
# this step integrated betwwen 2 samples and pre-step before integrated and annotation 
#on new integrated results
library(Seurat)
library(harmony)

# 1) تحديد مصدر كل سامبل
sample_d1$sample_id <- "CIN_39"
sample_d2$sample_id <- "CIN_34"

# حفظ الـ annotation القديم (SingleR) عشان تراجعه بعد الدمج
sample_d1$old_annotation <- sample_d1$SingleR_annotation
sample_d2$old_annotation <- sample_d2$SingleR_annotation

# 2) الدمج
merged_obj <- merge(sample_d1, y = sample_d2,
                    add.cell.ids = c("D1_S39", "D2_S34"),
                    merge.data = TRUE)

# 3) إعادة Normalization / VariableFeatures / Scaling على الـ merged object
merged_obj <- NormalizeData(merged_obj, normalization.method = "LogNormalize", scale.factor = 10000)
merged_obj <- FindVariableFeatures(merged_obj, selection.method = "vst", nfeatures = 2000)
merged_obj <- ScaleData(merged_obj, features = VariableFeatures(merged_obj))

# 4) PCA على الـ merged data
merged_obj <- RunPCA(merged_obj, features = VariableFeatures(merged_obj))

# 5) Integration بـ Harmony
merged_obj <- RunHarmony(merged_obj, group.by.vars = "sample_id")

# 6) Clustering و UMAP على الـ corrected embedding (استخدمت dims 1:10 زي ما كنت عامل قبل كذا)
merged_obj <- FindNeighbors(merged_obj, reduction = "harmony", dims = 1:10)
merged_obj <- FindClusters(merged_obj, resolution = 0.5)
merged_obj <- RunUMAP(merged_obj, reduction = "harmony", dims = 1:10)

# 7) تقييم الدمج - قبل/بعد
DimPlot(merged_obj, reduction = "umap", group.by = "sample_id") + ggtitle("Sample Mixing after Harmony")
DimPlot(merged_obj, reduction = "umap", group.by = "old_annotation", label = TRUE, repel = TRUE) +
  ggtitle("Old SingleR Annotation on Integrated UMAP")
DimPlot(merged_obj, reduction = "umap", label = TRUE) + ggtitle("New Clusters on Integrated Data")

merged_obj = JoinLayers(merged_obj)
library(SingleR)
library(scrapper)
pred = SingleR(test = GetAssayData(merged_obj, layer = "data"),
               ref = ref, labels = ref$label.main,
               clusters = merged_obj$seurat_clusters)

merged_obj$new_annotation = pred$labels[match(merged_obj$seurat_clusters, rownames(pred))]

DimPlot(merged_obj, reduction = "umap", group.by = "new_annotation",
        label = TRUE, repel = TRUE)

D1=DimPlot(merged_obj, reduction = "umap", group.by = "old_annotation", label = TRUE, repel = TRUE) +
  ggtitle("Old SingleR Annotation on Integrated UMAP")

D2=DimPlot(merged_obj, reduction = "umap", group.by = "new_annotation",
           label = TRUE, repel = TRUE)

D1+D2
# DotPlot to show the Expression
markers <- FindAllMarkers(merged_obj, only.pos = TRUE)
top3 <- markers |> dplyr::group_by(cluster) |> dplyr::slice_max(avg_log2FC, n = 3)
DotPlot(merged_obj, features = unique(top3$gene)) + RotatedAxis()
# if you want to show expression of specific cell (ex: Macrophage)
macro <- c("CD68","CD14","CD163","CSF1R","C1QA","C1QB","MRC1","MARCO","AIF1","LYZ")

DotPlot(merged_obj, features = macro) + RotatedAxis()

# عشان تزود ال resolution of annotation
#1. احتفظ بالـ annotation القديم (اختياري)
merged_obj$annotation_old_res = merged_obj$new_annotation

#2. ابنِ الجيران (neighbors) على Harmony
merged_obj = FindNeighbors(merged_obj, reduction = "harmony", dims = 1:30)

#3. زوّد الـ resolution
merged_obj = FindClusters(merged_obj, resolution = 1.2)

#4. حدّث الـ UMAP
merged_obj = RunUMAP(merged_obj, reduction = "harmony", dims = 1:30)

#5. ارسم الـ clusters الجديدة (أرقام)
DimPlot(merged_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

#6. شغّل SingleR على الـ clusters الجديدة
pred = SingleR(test = GetAssayData(merged_obj, layer = "data"),
               ref = ref, labels = ref$label.main,
               clusters = merged_obj$seurat_clusters)

#7. حط الأسماء على الخلايا
merged_obj$new_annotation = pred$labels[match(merged_obj$seurat_clusters, rownames(pred))]

#8. ارسم الأسماء
DimPlot(merged_obj, reduction = "umap", group.by = "new_annotation",
        label = TRUE, repel = TRUE)



saveRDS(merged_obj, "merged_CIN_39_34_integrated.rds")

# Annotation (Ref Based)
library(SingleR)
library(celldex)

hpca_ref <- celldex::HumanPrimaryCellAtlasData()

# تحويل الـ merged object لـ SingleCellExperiment
SCE_merged <- as.SingleCellExperiment(merged_obj)

# تشغيل SingleR
singler_merged <- SingleR(
  test = SCE_merged,
  ref = hpca_ref,
  labels = hpca_ref$label.main
)

# إضافة النتيجة للـ Seurat object
merged_obj$SingleR_annotation_integrated <- singler_merged$labels

# Visualization
DimPlot(merged_obj,
        reduction = "umap",
        group.by = "SingleR_annotation_integrated",
        label = TRUE,
        repel = TRUE) +
  ggtitle("Integrated Object — SingleR Automated Annotation")

###############################################################################
# Downstream

# ==========================================================
# 0) Setup - تحويل merged_obj لـ SingleCellExperiment
# ==========================================================
library(miloR)
library(SingleCellExperiment)
library(dplyr)

gastric_sce <- as.SingleCellExperiment(merged_obj, assay = "RNA")

# شوف اسم الـ reduced dimension الصحيح
reducedDimNames(gastric_sce)
# غالبًا هيظهر "HARMONY" - لو الاسم مختلف، استخدم الاسم الحقيقي في كل اللي جاي

# فلترة أي NA من الـ corrected embedding
gastric_sce <- gastric_sce[, apply(reducedDim(gastric_sce, "HARMONY"), 1,
                                   function(x) !all(is.na(x)))]

# ==========================================================
# 1) Build Milo Object + Neighborhoods
# ==========================================================
gastric_milo <- Milo(gastric_sce) %>%
  buildGraph(k = 30, d = 30, reduced.dim = "HARMONY") %>%
  makeNhoods(prop = 0.1, k = 30, d = 30, refined = TRUE, reduced_dims = "HARMONY") %>%
  countCells(meta.data = as.data.frame(colData(.)), sample = "sample_id") %>%
  calcNhoodDistance(d = 30, reduced.dim = "HARMONY")

# ==========================================================
# 2) شوف المتغيرات المتاحة قبل الـ Design Matrix
# ==========================================================
colnames(colData(gastric_milo))
# ==========================================================
# 3) Construct the Experimental Design Matrix
# ==========================================================

# استخراج المتغيرات المطلوبة للتحليل الإحصائي
gastric_design <- distinct(data.frame(colData(gastric_milo))[, c("sample_id"), drop = FALSE])

rownames(gastric_design) <- gastric_design$sample_id

gastric_design$sample_id <- factor(
  gastric_design$sample_id,
  levels = c("CIN_39", "CIN_34")
)

gastric_design
# تحويلها لـ data.frame صحيح (لأن عندنا متغيّر واحد بس، بيرجع vector مش df)
gastric_design <- data.frame(sample_id = gastric_design)

# تحديد أسامي الصفوف = sample IDs (شرط أساسي في miloR)
rownames(gastric_design) <- gastric_design$sample_id

# تحديد CIN_39 كـ reference level و CIN_34 كـ comparison level
gastric_design$sample_id <- factor(
  gastric_design$sample_id,
  levels = c("CIN_39", "CIN_34")
)

################################################################################  
# ==========================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)

# ==========================================================
# Define Cell Type + Sample Groups
# ==========================================================
merged_obj$celltype.sample <- paste(
  merged_obj$SingleR_annotation_integrated,
  merged_obj$sample_id,
  sep = "_"
)

Idents(merged_obj) <- "celltype.sample"

# تأكد من الأسامي المتاحة قبل ما تكتب ident.1 / ident.2
table(merged_obj$celltype.sample)

# ==========================================================
# Differential Expression: Macrophage CIN_39 vs CIN_34
# ==========================================================
DefaultAssay(merged_obj) <- "RNA"

macro.de <- FindMarkers(
  merged_obj,
  ident.1 = "Macrophage_CIN_39",
  ident.2 = "Macrophage_CIN_34",
  verbose = FALSE
)

head(macro.de, n = 10)

# ==========================================================
# Prepare Results for Visualization
# ==========================================================
macro.de$gene <- rownames(macro.de)

macro.de$significance <- "Not Significant"

macro.de$significance[
  macro.de$p_val_adj < 0.05 &
    macro.de$avg_log2FC > 0.25
] <- "Upregulated"

macro.de$significance[
  macro.de$p_val_adj < 0.05 &
    macro.de$avg_log2FC < -0.25
] <- "Downregulated"

top_upregulated <- macro.de %>% 
  filter(significance == "Upregulated") %>% 
  top_n(n = 10, wt = avg_log2FC) %>% 
  pull(gene)

top_downregulated <- macro.de %>% 
  filter(significance == "Downregulated") %>% 
  top_n(n = -10, wt = avg_log2FC) %>% 
  pull(gene)

genes_to_label <- c(top_upregulated, top_downregulated)

macro.de$label <- ifelse(macro.de$gene %in% genes_to_label, macro.de$gene, NA)

head(macro.de)

# ==========================================================
# Volcano Plot
# ==========================================================
colnames(macro.de)
class(macro.de)
sum(is.na(macro.de$p_val_adj))
table(macro.de$significance, useNA = "ifany")

macro.de$label <- ifelse(
  macro.de$significance != "Not Significant" & macro.de$p_val_adj < 0.05,
  macro.de$gene,
  NA
)

ggplot(
  macro.de,
  aes(
    x = avg_log2FC,
    y = -log10(p_val_adj),
    color = significance
  )
) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(
    aes(label = label),
    size = 3.5,
    box.padding = 0.5,
    max.overlaps = 30,
    na.rm = TRUE,
    color = "black",
    show.legend = FALSE
  ) +
  geom_vline(
    xintercept = c(-0.25, 0.25),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Differential Expression: Macrophage CIN_39 vs CIN_34",
    x = "Average Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  scale_color_manual(
    values = c(
      "Downregulated" = "red",
      "Not Significant" = "grey",
      "Upregulated" = "blue"
    )
  ) +
  theme_classic()

# ==========================================================
# Save Complete DE Results
# ==========================================================
write.csv(
  macro.de,
  "Macrophage_CIN39_vs_CIN34_DEG_results.csv",
  row.names = FALSE
)

# ==========================================================
# Filter Significant DEGs
# ==========================================================
deg <- subset(
  macro.de,
  p_val_adj < 0.05 &
    abs(avg_log2FC) > 0.25
)

head(deg, 10)

# ==========================================================
# Separate Upregulated and Downregulated Genes
# ==========================================================
up_genes <- rownames(
  subset(
    macro.de,
    p_val_adj < 0.05 &
      avg_log2FC > 0.25
  )
)

down_genes <- rownames(
  subset(
    macro.de,
    p_val_adj < 0.05 &
      avg_log2FC < -0.25
  )
)

head(up_genes)
head(down_genes)

write.csv(
  data.frame(Gene = up_genes),
  "UPregulated_genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = down_genes),
  "DOWNregulated_genes.csv",
  row.names = FALSE
)
################################################################################
#Cell-Cell Communication
library(CellChat)
cellchat <- createCellChat(object = merged_obj, group.by = "SingleR_annotation_integrated")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways

# SPP1 Signaling Pathway - مرتبط مباشرة بفرضية الورقة البحثية
netVisual_aggregate(
  cellchat,
  signaling = "SPP1",
  layout = "circle"
)
df_net = subsetCommunication(cellchat)
spp1_lr = df_net[df_net$ligand == "SPP1", ]
head(spp1_lr)

spp1_lr <- subsetCommunication(cellchat, signaling = "SPP1")
head(spp1_lr)

# فلترة بس على Macrophage → T_cells
spp1_macro_tcell <- spp1_lr[spp1_lr$source == "Macrophage" & spp1_lr$target == "T_cells", ]
spp1_macro_tcell
################################################################################
# 3: Functional Analysis (GO and KEGG)
#(A) GO: 
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
up_df <- bitr(
  up_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

head(up_df)
# قبل ما تكمل: شوف كام جين اتحوّل من الإجمالي:
nrow(up_df)
length(up_genes)
# لو الرقمين قريبين، تمام.

#GO Biological Process — Upregulated Genes
ego_up_BP <- enrichGO(
  gene = up_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",# Benjamini-Hochberg (BH)
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(as.data.frame(ego_up_BP), 10)
#Visualize
dotplot(
  ego_up_BP,
  showCategory = 15
) +
  ggtitle("GO Biological Process - UPregulated Genes")
#GO Biological Process — Downregulated Genes
down_df <- bitr(
  down_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

head(down_df)

ego_down_BP <- enrichGO(
  gene = down_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(as.data.frame(ego_down_BP), 10)
#Visualize
dotplot(
  ego_down_BP,
  showCategory = 15
) +
  ggtitle("GO Biological Process - DOWNregulated Genes")
#GO Molecular Function — Upregulated Genes
ego_up_MF <- enrichGO(
  gene = up_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(as.data.frame(ego_up_MF), 10)

dotplot(
  ego_up_MF,
  showCategory = 15
) +
  ggtitle("GO Molecular Function - UPregulated Genes")
#GO Cellular Component — Upregulated Genes
ego_up_CC <- enrichGO(
  gene = up_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(as.data.frame(ego_up_CC), 10)

dotplot(
  ego_up_CC,
  showCategory = 15
) +
  ggtitle("GO Cellular Component - UPregulated Genes")
# (B) KEGG Pathway Enrichment — Upregulated Genes
ekegg_up <- enrichKEGG(
  gene = up_df$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)
# if this code don't work, try this will work
library(msigdbr)
kegg_sets = msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
kegg_t2g = kegg_sets[, c("gs_name", "gene_symbol")]
ekegg_up = enricher(gene = up_genes,
                    universe = rownames(macro.de),
                    TERM2GENE = kegg_t2g,
                    pvalueCutoff = 0.05)

head(as.data.frame(ekegg_up), 10)
# visualize 
dotplot(
  ekegg_up,
  showCategory = 15
) +
  ggtitle("KEGG Pathway Enrichment - UPregulated Genes")
#We can also convert the KEGG gene IDs back to readable gene symbols for easier interpretation.
ekegg_up <- setReadable(
  enrichKEGG(
    gene = up_df$ENTREZID,
    organism = "hsa",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  ),
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

head(as.data.frame(ekegg_up), 10)
dotplot(
  ekegg_up,
  showCategory = 15
) +
  ggtitle("KEGG Pathway Enrichment - UPregulated Genes")
#KEGG Pathway Enrichment — Downregulated Genes
ekegg_down <- enrichKEGG(
  gene = down_df$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)
down_genes = macro.de %>%
  filter(significance == "Downregulated") %>%
  pull(gene)
ekegg_down = enricher(gene = down_genes,
                      universe = rownames(macro.de),
                      TERM2GENE = kegg_t2g,
                      pvalueCutoff = 0.05)

head(as.data.frame(ekegg_down), 10)
dotplot(ekegg_down, showCategory = 10, font.size = 10, label_format = 40) +
  ggtitle("KEGG Pathways - Downregulated Genes")


# Convert Entrez IDs back to gene symbols for readability
if (!is.null(ekegg_down)) {
  ekegg_down <- setReadable(
    ekegg_down,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID"
  )
}

head(as.data.frame(ekegg_down), 10)

if (!is.null(ekegg_down) && nrow(as.data.frame(ekegg_down)) > 0) {
  dotplot(
    ekegg_down,
    showCategory = 15
  ) +
    ggtitle("KEGG Pathway Enrichment - DOWNregulated Genes")
}


# Save GO Results 
write.csv(
  as.data.frame(ego_up_BP),
  "GO_BP_UPregulated.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_down_BP),
  "GO_BP_DOWNregulated.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_up_MF),
  "GO_MF_UPregulated.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_up_CC),
  "GO_CC_UPregulated.csv",
  row.names = FALSE
)
# Save KEGG Results
write.csv(
  as.data.frame(ekegg_up),
  "KEGG_UPregulated.csv",
  row.names = FALSE
)

if (exists("ekegg_down") && !is.null(ekegg_down)) {
  write.csv(
    as.data.frame(ekegg_down),
    "KEGG_DOWNregulated.csv",
    row.names = FALSE
  )
} else {
  message("No significant KEGG pathways found for downregulated genes.")
}

