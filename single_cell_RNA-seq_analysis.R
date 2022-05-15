rm(list=ls())
gc()

library(Seurat)
library(Matrix)
library(dplyr) # Dataframe manipulation

library(monocle)
library(DDRTree)
library(data.table)
library(pheatmap)
library(corrplot) 

library(dplyr)
library(ggsci
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(cowplot)
library(patchwork)

library(ggpubr)
library(ggthemes)
library(RColorBrewer)
#############################
#######sample-object construct########
###### use CS10_YS as an example######
#YS10-Ery
CS10_YS_data <- Read10X("./CS10_YS/outs/filtered_feature_bc_matrix/")
CS10_YS_meta <- read.csv("./CS10_YS/CS10_YS_doublet_annot.csv",header = T,row.names = 1)
#####Creat FL-Ery project######
CS10_YS <- CreateSeuratObject(counts = CS10_YS_data, project = "CS10_YS", 
                                    meta.data = CS10_YS_meta)
#remove potential doubles
Idents(CS10_YS) <- CS10_YS@meta.data$man_doublets
CS10_YS_singlet <- subset( x = CS10_YS, idents = "FALSE")
#remove low quality cells
CS10_YS@meta.data$samples <- "CS10_YS"
CS10_YS@meta.data$Stage <- "YS"
CS10_YS[["percent.mt"]] <- PercentageFeatureSet(CS10_YS, pattern = "^MT-")
CS10_YS <- subset(CS10_YS, subset = nFeature_RNA > 500 & percent.mt < 10)   
CS10_YS <- NormalizeData(object = CS10_YS)
CS10_YS <- FindVariableFeatures(object = CS10_YS)

CS10_YS <- CellCycleScoring(
  object = CS10_YS,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)

#############################
#######stage-object construct########
###### use YS as an example######
#Use CCA in Seurat for data integrate
YS.anchors <- FindIntegrationAnchors(object.list = c(CS10_YS,CS11_YS,CS15_YS), dims = 1:30)
YS.combined <- IntegrateData(anchorset = YS.anchors, dims = 1:30)

DefaultAssay(YS.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
all.genes <- rownames(YS.combined)
YS.combined <- ScaleData(YS.combined, features = all.genes, vars.to.regress =c("nCount_RNA","percent.mt","S.Score","G2M.Score"))
YS.combined <- RunPCA(YS.combined, npcs = 50, features = VariableFeatures(object = YS.combined))

ElbowPlot(YS.combined)
dim(YS.combined@meta.data)
YS.combined <- RunUMAP(YS.combined, reduction = "pca", dims = 1:15)
YS.combined <- RunTSNE(YS.combined, reduction = "pca", dims = 1:15)
YS.combined <- FindNeighbors(YS.combined, reduction = "pca", dims = 1:15)
YS.combined <- FindClusters(YS.combined, resolution = 0.7)
table(YS.combined@active.ident)

DefaultAssay(YS.combined) <- "RNA"
YS.combined_markers <- FindAllMarkers( object = YS.combined,only.pos = T,
                                       test.use = "wilcox")
YS.combined_markers <- YS.combined_markers[YS.combined_markers$p_val_adj < 0.01,]
YS.combined_markers <- YS.combined_markers %>% group_by(gene) %>% top_n(n = 1, wt = avg_logFC)

#############################
#######Integrated-object construct########
YS_Ery <- subset(x=YS.combined, idents="Ery")
FL_Ery <- subset(x=FL.combined, idents="Ery")
UCB_Ery <- subset(x=UCB.combined, idents="Ery")
Ery_10X_combined <- merge( x = YS_Ery, y = c(FL_Ery,UCB_Ery),
                           add.cell.ids = c("YS","FL","UCB"))
Ery.list <- SplitObject(Ery_10X_combined, split.by = "samples")

#########################
for (i in names(Ery.list)) {
  Ery.list[[i]] <- SCTransform(Ery.list[[i]],variable.features.n = 3000,
                               verbose = FALSE)
}
Ery.features <- SelectIntegrationFeatures(object.list = Ery.list, nfeatures = 2000)
Ery.list <- PrepSCTIntegration(object.list = Ery.list, anchor.features = Ery.features)
Ery.anchors <- FindIntegrationAnchors(object.list = Ery.list, normalization.method = "SCT", 
                                      anchor.features = Ery.features)

Ery.integrated <- IntegrateData(anchorset = Ery.anchors, normalization.method = "SCT")

DefaultAssay(Ery.integrated) <- "integrated"
ElbowPlot(Ery.integrated)

Ery.integrated <- RunPCA(object = Ery.integrated, verbose = FALSE,npcs = 50)
Ery.integrated <- RunUMAP(object = Ery.integrated, dims = 1:6,min.dist = 0.16)
Ery.integrated <- RunTSNE(object = Ery.integrated, dims = 1:6)

Ery.integrated <- FindNeighbors(Ery.integrated, reduction = "pca", dims = 1:6)
Ery.integrated <- FindClusters(Ery.integrated, resolution = 0.1)
table(Ery.integrated@active.ident)

#############################
#######Figure1-Sfigure1######
load("./Ery.integrated.RData")
cluster_col <- c("#1F77B4","#FF7F0E","#8C564B","#D62728")
stage_col <- c("#461652","#FAE40E","#53C7DC")
phase_col <- c("#0E1674","#FF4A46","#4C9A25")

DimPlot( object = Ery.integrated,reduction = "umap",pt.size = 1,cols = cluster_col, 
         split.by = "samples",ncol = 3)
VlnPlot(Ery.integrated, features = c("nCount_RNA"), pt.size = 0.1,cols = cluster_col)		 

#####DEGs identification##########
DefaultAssay(Ery.integrated) <- "RNA"
Ery.integrated.markers <- FindAllMarkers( object = Ery.integrated,only.pos = T,logfc.threshold = 0.25,
                               test.use = "wilcox")
Ery.integrated.markers <- Ery.integrated.markers[Ery.integrated.markers$p_val_adj < 0.01,]

Ery.integrated.markers <- Ery.integrated.markers %>% group_by(gene) %>% top_n(n = 1, wt = avg_logFC)

genes_plot <- c("RPL23","RPS18","EIF3E","SRSF7","FAM178B","MKI67","CDK1","E2F4","CDC27",
                "CCNB1","BNIP3L","HBA1","ALAS2","GABARAPL2","MAP1LC3B","SRGN","NFKBIA","FOS","IFITM3","LYZ")
DotPlot( object = Ery.integrated, features = as.character(genes_plot),cols = c("lightgrey","red"))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ coord_flip()
###############
##########G0-Terms###########
gene_list <- read.csv("./Ery.integrated.markers.csv",header = T)
ids_changed <- bitr(gene_list$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
ids_changed_genes <- merge(gene_list,ids_changed,by.x="gene", by.y="SYMBOL",sort=F)
#######################
GO_Terms_Cluster <- compareCluster(gene~cluster, data=ids_changed_genes,ont = "BP", keyType="SYMBOL",
                                   fun='enrichGO', OrgDb='org.Hs.eg.db')

######correlation############
Ery.cor <- cor(Ery_average_genes_merge,method = "spearman")
col_pie <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                           "yellow", "#FF7F00", "red", "#7F0000"))
corrplot(Ery.cor, order = "original", method = "pie",addgrid.col = "black",
         addCoef.col = "black",tl.col = "black",col = col_pie)

############
genesets <- read.csv("./genesets_ery.csv",header = T)
head(genesets)
Erythrocyte_maturation  <- list(as.character(genesets$GO_ERYTHROCYTE_MATURATION )[1:15])

DefaultAssay(Ery.integrated) <- "RNA"
Ery.integrated <- AddModuleScore( object = Ery.integrated,
                          features = Erythrocyte_maturation,
                          name = "Erythrocyte_maturation")

#############
meta_info_ery <- as.data.frame(Ery.integrated@meta.data)

meta_info_ery$stage <- factor(meta_info_ery$stage,levels=c("YS","FL","UCB"),
                              ordered=T) 
ggplot(meta_info_ery, aes(x = stage, y = 100, fill = defined)) +
  geom_bar(stat = "identity", position = 'fill')+scale_fill_manual(values= cluster_col)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(meta_info_ery, aes(x = defined, y = Erythrocyte_maturation  ))+
	geom_boxplot(aes(fill=defined),notch= FALSE)+ scale_fill_manual(values = cluster_col)+
	theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.line = element_line(colour = "black"))

###########trajectory analysis of C4 cells########
Ery.integrated_C4 <- subset( x = Ery.integrated, idents = "C4")
Ery_Pathway_rawdata <- as.matrix(Ery.integrated_C4@assays$RNA@counts)
#######
gene_ann <- data.frame(gene_short_name = row.names(Ery_Pathway_rawdata), row.names = row.names(Ery_Pathway_rawdata))
head(gene_ann)
pd <- new("AnnotatedDataFrame",data=cell.meta.data)
fd <- new("AnnotatedDataFrame",data=gene_ann)
Ery_Pathway_cds <- newCellDataSet(Ery_Pathway_rawdata, 
                                  phenoData = pd,
                                  featureData =fd,
                                  expressionFamily = negbinomial.size(),
                                  lowerDetectionLimit=0.1)
Ery_Pathway_cds<- estimateSizeFactors(Ery_Pathway_cds)
Ery_Pathway_cds<- estimateDispersions(Ery_Pathway_cds)

Ery_Pathway_cds <- detectGenes(Ery_Pathway_cds, min_expr = 0.5)

expressed_genes <- row.names(subset(fData(Ery_Pathway_cds),
                                    num_cells_expressed > 1))
print(head(pData(Ery_Pathway_cds)))

####choose ordering genes based on PCA loading
Ery_Pathway_cds <- Ery_Pathway_cds[expressed_genes,]
exprs_filtered <- t(t(exprs(Ery_Pathway_cds)/pData(Ery_Pathway_cds)$Size_Factor))
#nz_genes <- which(exprs_filtered != 0, arr.ind = T)
#exprs_filtered@x <- log(exprs_filtered@x + 1)
# Calculate the variance across genes without converting to a dense
# matrix:
expression_means <- Matrix::rowMeans(exprs_filtered)
expression_vars <- Matrix::rowMeans((exprs_filtered - expression_means)^2)
# Filter out genes that are constant across all cells:
genes_to_keep <- expression_vars > 0
exprs_filtered <- exprs_filtered[genes_to_keep,]
expression_means <- expression_means[genes_to_keep]
expression_vars <- expression_vars[genes_to_keep]
# Here's how to take the top PCA loading genes, but using
# sparseMatrix operations the whole time, using irlba. Note
# that the v matrix from irlba is the loading matrix
set.seed(0)
irlba_pca_res <- irlba(t(exprs_filtered),
                       nu=0,
                       center=expression_means,
                       scale=sqrt(expression_vars),
                       right_only=TRUE)$v
row.names(irlba_pca_res) <- row.names(exprs_filtered)
# Here, we will just
# take the top 200 genes from components 2 and 3.
# Component 1 usually is driven by technical noise.
# We could also use a more principled approach,
# similar to what dpFeature does below
PC2_genes <- names(sort(abs(irlba_pca_res[, 2]), decreasing = T))[1:200]
PC3_genes <- names(sort(abs(irlba_pca_res[, 3]), decreasing = T))[1:200]

ordering_genes <- union(PC2_genes,PC3_genes)###used for C4 pseudotime analysis

##############################
Ery_Pathway_cds <- setOrderingFilter(Ery_Pathway_cds, ordering_genes)

Ery_Pathway_cds <- reduceDimension(Ery_Pathway_cds, max_components = 2,
                                   method = 'DDRTree')

Ery_Pathway_cds <- orderCells(Ery_Pathway_cds)

######pseudotime#######
P <- plot_cell_trajectory(Ery_Pathway_cds, color_by = "Pseudotime", show_tree = T,
                          cell_size = 1,show_branch_points=F, return=T,theta = 180)+theme(legend.position=c(.1,.8))
P1 <- P + scale_colour_gradientn(colours = hcl.colors(1000), values = NULL, space = "Lab",
                         na.value = "grey50", guide = "colourbar", aesthetics = "colour")
P1
##########
my_genes <- row.names(subset(fData(Ery_Pathway_cds),
                             gene_short_name %in% c("HBB","HBA2")))
subset <- Ery_Pathway_cds[my_genes,]
P <- plot_genes_in_pseudotime(subset,color_by="defined",ncol = 2,
                         cell_size = 1,panel_order = c("HBB","HBA2"))

#############################
#######Figure2-Sfigure2######						
C1 <- subset( x = Ery.integrated, idents = "C1")
DefaultAssay(C1) <- "RNA"
Idents(C1) <- C1@meta.data$stage
levels( x = C1) <- c("YS","FL","UCB")

C1.markers <- FindAllMarkers( object = C1,only.pos = T,logfc.threshold = 0.25,
                              test.use = "wilcox")
C1.markers <- C1.markers[C1.markers$p_val_adj < 0.01,]
C1.markers <- C1.markers %>% group_by(gene) %>% top_n(n = 1, wt = avg_logFC)
###############
MTs <- list(as.character(genesets$MTs[1:11]))
Ribosome <- list(as.character(genesets$GO_RIBOSOME[1:228]))
Glycolysis <- list(as.character(genesets$Glycolysis[1:200]))
G2M_Transition <- list(as.character(genesets$GO_CELL_CYCLE_G2_M_PHASE_TRANSITION[1:273]))

Ery.integrated <- AddModuleScore( object = Ery.integrated,
                          features = MTs,
                          name = "MTs")

##########Primitive vs Definitive############
DefaultAssay(Ery.integrated) <- "RNA"
Idents(Ery.integrated) <- Ery.integrated@meta.data$stage
YS_Ery <- subset( x= Ery.integrated, idents = "YS")
FL_Ery <- subset( x= Ery.integrated, idents = "FL")

#####Primitive Ery
YS_Ery_norm <- as.data.frame(YS_Ery@assays$RNA@data)
YS_Ery_norm <- as.data.frame(t(YS_Ery_norm))
YS_Ery_norm_HBE1_HBG1 <- YS_Ery_norm[,c("HBE1","HBG1","HBG2")]
YS_Ery_norm_HBE1_HBG1$G <- (YS_Ery_norm_HBE1_HBG1$HBG1+YS_Ery_norm_HBE1_HBG1$HBG2)/2
YS_Ery_norm_HBE1_HBG1$E <- YS_Ery_norm_HBE1_HBG1$HBE1
YS_Ery_norm_HBE1_HBG1$E_G <- YS_Ery_norm_HBE1_HBG1$HBE1 - YS_Ery_norm_HBE1_HBG1$G
Pri_YS_ery <- subset(YS_Ery_norm_HBE1_HBG1, subset =  E_G > 0 )
Pri_YS_ery_cells <- row.names(Pri_YS_ery)
Ery.integrated <- SetIdent(object = Ery.integrated, cells = Pri_YS_ery_cells, 
                           value = 'YS_Pri')
#####Definitive Ery
FL_Ery_norm <- as.data.frame(FL_Ery@assays$RNA@data)
FL_Ery_norm <- as.data.frame(t(FL_Ery_norm))
FL_Ery_norm_HBE1_HBG1 <- FL_Ery_norm[,c("HBE1","HBG1","HBG2")]
FL_Ery_norm_HBE1_HBG1$G <- (FL_Ery_norm_HBE1_HBG1$HBG1+FL_Ery_norm_HBE1_HBG1$HBG2)/2
FL_Ery_norm_HBE1_HBG1$E <- FL_Ery_norm_HBE1_HBG1$HBE1
FL_Ery_norm_HBE1_HBG1$E_G <- FL_Ery_norm_HBE1_HBG1$HBE1 - FL_Ery_norm_HBE1_HBG1$G
Def_FL_ery <- subset(FL_Ery_norm_HBE1_HBG1, subset =  E_G < 0 )
Def_FL_ery_cells <- row.names(Def_FL_ery)
Ery.integrated <- SetIdent(object = Ery.integrated, cells = Def_FL_ery_cells, 
                           value = 'FL_Def')
Ery_pri_def <- subset( x = Ery.integrated, idents = c("YS_Pri","FL_Def"))
levels(Ery_pri_def) <- c("YS_Pri","FL_Def")

#####construct object###
Ery_pri_def_norm <- as.data.frame(Ery_pri_def@assays$RNA@data)
Ery_pri_def_norm <- as.data.frame(t(Ery_pri_def_norm))
Ery_pri_def_norm_HBE1_HBG1 <- Ery_pri_def_norm[,c("HBE1","HBG1","HBG2")]
Ery_pri_def_norm_HBE1_HBG1$G <- (Ery_pri_def_norm$HBG1+Ery_pri_def_norm$HBG2)/2
Ery_pri_def_norm_HBE1_HBG1$E <- Ery_pri_def_norm$HBE1
Ery_pri_def_HBE_HBG <- Ery_pri_def_norm_HBE1_HBG1[,c("G","E")]
colnames(Ery_pri_def_HBE_HBG) <- c("UMAP_1","UMAP_2")
Ery_pri_def@reductions$umap@cell.embeddings <- as.matrix(Ery_pri_def_HBE_HBG) 
DimPlot(Ery_pri_def,pt.size = 1)

###DEGs###
DefaultAssay(object = Ery_pri_def) <- "RNA"
Ery_pri_def.markers <- FindAllMarkers( object = Ery_pri_def,only.pos = T,logfc.threshold = 0.25,
                                          test.use = "wilcox")
Ery_pri_def.markers <- Ery_pri_def.markers[Ery_pri_def.markers$p_val_adj < 0.01,]
Ery_pri_def.markers <- Ery_pri_def.markers %>% group_by(gene) %>% top_n(n = 1, wt = avg_logFC)
######TFs amd FLs##########
TF_list <- read.csv("./TF_database.csv",header = T)
(TF_list)
SF_list <- read.csv("./Surface_marker.csv",header = T)

Ery_pri_def_TFs <- merge( x = Ery_pri_def.markers, y = TF_list, 
                         by.x = "gene", by.y = "TF_genes", sort = F)

Ery_pri_def_SFs <- merge( x = Ery_pri_def.markers, y = SF_list, 
                         by.x = "gene", by.y = "gene_symbol", sort = F)

#############################
#######Figure3-Sfigure3######	
load("./Ery.integrated.RData")					
load("./ES_Ery.combined.RData")
Vivo_vitro_merge_object <- merge(Ery.integrated,ES_Ery.combined)
Vivo_vitro_merge_object@meta.data$defined_stage <- paste(Vivo_vitro_merge_object@meta.data$stage, Vivo_vitro_merge_object@meta.data$defined, sep = "_")
Idents(Vivo_vitro_merge_object) <- Vivo_vitro_merge_object@meta.data$defined_stage

C1_YS_Merge <- subset( x = Vivo_vitro_merge_object, idents = c("YS_C1","hESC_hESC-Ery1"))
C1_FL_Merge <- subset( x = Vivo_vitro_merge_object, idents = c("FL_C1","hESC_hESC-Ery1"))

C2_YS_Merge <- subset( x = Vivo_vitro_merge_object, idents = c("YS_C2","hESC_hESC-Ery2"))
C2_FL_Merge <- subset( x = Vivo_vitro_merge_object, idents = c("FL_C2","hESC_hESC-Ery2"))

C3_YS_Merge <- subset( x = Vivo_vitro_merge_object, idents = c("YS_C3","hESC_hESC-Ery3"))
C3_FL_Merge <- subset( x = Vivo_vitro_merge_object, idents = c("FL_C3","hESC_hESC-Ery3"))

###########DEGs between YS-hESC or FL-hESC#########
vivo_hESC_DEG <- FindAllMarkers(C1_YS_Merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                                 test.use = "wilcox")
vivo_hESC_DEG <-vivo_hESC_DEG[vivo_hESC_DEG$p_val_adj<0.01,]
vivo_hESC_DEG <- vivo_hESC_DEG %>% group_by(gene) %>% top_n(n = 1, wt = avg_logFC)

#############################
#######Figure4-Sfigure4######	
load("./Ery.integrated.RData")					
Immune_active <- list(as.character(genesets$GO_ACTIVATION_OF_IMMUNE_RESPONSE)[1:562])
Cytokine_pathway <- list(as.character(genesets$GO_CYTOKINE_MEDIATED_SIGNALING_PATHWAY)[1:791])
Ery.integrated <- AddModuleScore( object = Ery.integrated,
                          features = Immune_active,
                          name = "Immune_active")
pdata <- as.data.frame(Ery.integrated@meta.data)
C1.meta <- pdata[pdata$defined == "C1",]
C2.meta <- pdata[pdata$defined == "C2",]
C3.meta <- pdata[pdata$defined == "C3",]
C4.meta <- pdata[pdata$defined == "C4",]

wilcox.test(C4.meta$Immune_active, C1.meta$Immune_active)
wilcox.test(C4.meta$Immune_active, C2.meta$Immune_active)
wilcox.test(C4.meta$Immune_active, C3.meta$Immune_active)

#########SFs of C4 cells########						  
SF_list <- read.csv("./Surface_marker.csv",header = T)

Ery_integrated_SFs <- merge( x = Ery.integrated.markers, y = SF_list, 
                         by.x = "gene", by.y = "gene_symbol", sort = F)

average_genes <- AverageExpression( object = Ery.integrated)
average_genes_integrated <- average_genes$RNA

DEG_average <-average_genes[as.character(unique(genes_plot)),]
pheatmap(DEG_average,scale = "row",cluster_rows = F,cluster_cols = F,
         color =  rev(RColorBrewer::brewer.pal(n = 8, name = "RdYlBu")))

####KEGG of C4 TFs target genes####### 		 
gene_list <- read.csv("./scenic/TFtarget_C4.csv",header = T)
ids_changed <- bitr(gene_list$target, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

ids_changed_genes <- merge(gene_list,ids_changed,by.x="target", by.y="SYMBOL",sort=F)
#######################
genes_target <- ids_changed_genes$ENTREZID
KEGG_target <- enrichKEGG(genes_target, organism = "hsa", keyType = "kegg",
           pvalueCutoff = 0.05, pAdjustMethod = "BH")

#####scRT-qPCR########
rna_data <- read.csv("./scRT_qPCR_dataset.csv",header = T)
t_data <- melt(rna_data)
ggplot(t_data, aes(x = variable, y = value  ))+
  geom_boxplot(aes(fill=variable),
               notch= FALSE)+   
  theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.line = element_line(colour = "black"))

wilcox.test( t_data$C1.C3.LPS,t_data$C4.LPS)
wilcox.test( t_data$PB.C1.C3.LPS,t_data$PB.C4.LPS)

#############################
#######Sfigure5##############			   
load("./pre_UCB_Ery.Rdata")
pre_UCB_Ery.markers <- FindAllMarkers( object = pre_UCB_Ery,only.pos = T,logfc.threshold = 0.25,
                                          test.use = "wilcox")
pre_UCB_Ery.markers <- pre_UCB_Ery.markers[pre_UCB_Ery.markers$p_val_adj < 0.01,]
pre_UCB_Ery.markers <- pre_UCB_Ery.markers %>% group_by(gene) %>% top_n(n = 1, wt = avg_logFC)

pre_UCB_Ery_average <- AverageExpression( object = pre_UCB_Ery)
pre_UCB_Ery_average_RNA <- pre_UCB_Ery_average$RNA
########
genes_plot <- c("FAM178B","RPS26","RPL4","EIF5A","SRSF7","MKI67","CDK1",
                "CCNB1","E2F4","CDC27","BNIP3L","HBA1","ALAS2",
                "GABARAPL2","MAP1LC3B","SRGN","FOS","IFITM3","LYZ","NFKBIA")
DEG_average <-pre_UCB_Ery_average_RNA[as.character(genes_plot),]
pheatmap(DEG_average,scale = "row",cluster_rows = F,cluster_cols = F,show_rownames = T,
         color = colorRampPalette(c("#0000CD","#F8F8FF","#FF0000"))(100))

TFs_plot <- c("GATA2","MEF2C","SPI1","JUN","FOS")
P <- DotPlot( object = pre_UCB_Ery, features = TFs_plot,
              cols = c("#F5F5F5","red"),scale.by = "size")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  coord_flip()
P+scale_colour_gradient(low = "white", high = "red",
                        na.value = NA,space = "Lab",guide = "colourbar",
                        aesthetics = "colour")
#####Jaccard similarity coefficient based on DEGs$$$$$$$$
YS_FL_pre.UCB_UCB_object <- merge(Ery.integrated,pre_UCB_Ery,
                                  add.cell.ids = c("int","pre"))

C4_stage_object <- subset( x = YS_FL_pre.UCB_UCB_object, idents = c("YS_C4","FL_C4",
                                                                    "pre_UCB_C4","UCB_C4"))

C4_stage.markers <- FindAllMarkers( object = C4_stage_object,only.pos = T,logfc.threshold = 0.25,
                                  test.use = "wilcox")
C4_stage.markers <- C4_stage.markers[C4_stage.markers$p_val_adj < 0.01,]
C4_stage.markers <- C4_stage.markers %>% group_by(gene) %>% top_n(n = 1, wt = avg_logFC)

jac_matrix <- matrix(nrow = 5, ncol = 5)
rownames(jac_matrix) <- sort(unique(C4_stage.markers$cluster))
colnames(jac_matrix) <- sort(unique(C4_stage.markers$cluster))
for (i in sort(unique(C4_stage.markers$cluster))) {
  for (j in sort(unique(C4_stage.markers$cluster))) {
    x = subset(C4_stage.markers, cluster == i)$gene
    y = subset(C4_stage.markers, cluster == j)$gene
    jac_tmp <- length(intersect(x, y))/length(unique(union(x, y)))
    jac_matrix[i,j] = jac_tmp
  }
}
pheatmap(t(jac_matrix),color = rev(RColorBrewer::brewer.pal(n = 8, name = "RdYlBu")),
         clustering_method = "mcquitty")
###################pre_UCB_C4 vs UCB_C4###############
############Differences between stages##########

C4_UCB <- FindMarkers( object = YS_FL_pre.UCB_UCB_object, ident.1 = "pre_UCB_C4",ident.2 = "UCB_C4",
                       logfc.threshold = 0,min.pct = 0)

deg.data <- as.data.frame(C4_UCB)
deg.data$gene <- rownames(deg.data)
deg.data$adjP <- -log10(deg.data$p_val_adj)
deg.data$Group = "Non-significant"

deg.data$Group[which((deg.data$p_val_adj<0.05) & (deg.data$avg_logFC > 0.25 ))] = "Up in pre_UCB_C4"
deg.data$Group[which((deg.data$p_val_adj<0.05) & (deg.data$avg_logFC < -0.25 ))] = "Up in UCB_C4"
table(deg.data$Group)

deg.data$Label = ""
deg.data <- deg.data[order(deg.data$adjP),]
pre_UCB_up_genes <- tail(deg.data$gene[which(deg.data$Group == "Up in pre_UCB_C4")],10)
UCB_up_genes <- tail(deg.data$gene[which(deg.data$Group == "Up in UCB_C4")],10)

deg.top10.genes <- c(as.character(pre_UCB_up_genes),as.character(UCB_up_genes))
deg.data$Label[match(deg.top10.genes,deg.data$gene)] <- deg.top10.genes
deg.data$Group <- factor(deg.data$Group,levels=c("Up in pre_UCB_C4","Non-significant","Up in UCB_C4"
),ordered=T) 

ggscatter(deg.data,x = "avg_logFC", y = "adjP",
          color = "Group",palette = c("#d5b60a","#BBBBBB","#984EA3"),size = 1,
          
          font.label = 15,
          repel = T,
          xlab = "log2FoldChange",
          ylab = "-log10(P-adjust)") +theme_base()+
  geom_hline(yintercept = 1.30103, linetype = "dashed")+
  geom_vline(xintercept = c(-0.25,0.25), linetype = "dashed")

#############################
#######Figure5-Sfigure6##############
load("./human_BM_Ery_smart_seq2.Rdata")
Immune_active <- list(as.character(genesets$GO_ACTIVATION_OF_IMMUNE_RESPONSE)[1:562])
Cytokine_pathway <- list(as.character(genesets$GO_CYTOKINE_MEDIATED_SIGNALING_PATHWAY)[1:791])
Erythrocyte_maturation  <- list(as.character(genesets$GO_ERYTHROCYTE_MATURATION )[1:15])

human_BM_Ery_smart_seq2 <- AddModuleScore( object = human_BM_Ery_smart_seq2,
                          features = Immune_active,
                          name = "Immune_active")
pdata <- as.data.frame(human_BM_Ery_smart_seq2@meta.data)
C1.meta <- pdata[pdata$defined == "C1",]
C2.meta <- pdata[pdata$defined == "C2",]
C3.meta <- pdata[pdata$defined == "C3",]
C4.meta <- pdata[pdata$defined == "C4",]

wilcox.test(C4.meta$Immune_active, C1.meta$Immune_active)
wilcox.test(C4.meta$Immune_active, C2.meta$Immune_active)
wilcox.test(C4.meta$Immune_active, C3.meta$Immune_active)

###########DEGs
DefaultAssay(human_BM_Ery_smart_seq2) <- "RNA"

BM.DEGs <- FindAllMarkers( object = Integrated_TDE, only.pos = T)
BM.DEGs <- BM.DEGs[BM.DEGs$p_val_adj < 0.01,]
BM.DEGs <- BM.DEGs %>% group_by(gene) %>% top_n(n = 1, wt = avg_logFC)

top10 <- BM.DEGs %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

DoHeatmap(human_BM_Ery_smart_seq2, features = as.character(top10$gene), disp.min = -1.5,disp.max = 1.5) + 
	scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 4, name = "RdYlBu")))

#########GO-Terms##########
gene_list <- BM.DEGs
ids_changed <- bitr(gene_list$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
ids_changed_genes <- merge(gene_list,ids_changed,by.x="gene", by.y="SYMBOL",sort=F)

#######################
GO_Terms_Cluster <- compareCluster(gene~cluster, data=ids_changed_genes,ont = "BP", keyType="SYMBOL",
                                   fun='enrichGO', OrgDb='org.Hs.eg.db')	

#############################
#######Figure6-Sfigure7##############
#####use FL as example########
data_cpDB_sig_mean <- read.csv(file = "./FL/significant_means_filtered.csv", header = T)
data_cpDB_pval <- read.csv(file = "./FL/pvalues_filtered.csv", header = T)

rownames(data_cpDB_sig_mean) <- data_cpDB_sig_mean$interacting_pair
data_cpDB_sig_mean <- data_cpDB_sig_mean[, -1]
data_cpDB_sig_mean[is.na(data_cpDB_sig_mean)] <- 0
range(data_cpDB_sig_mean)
data_cpDB_sig_mean_log <- log2(data_cpDB_sig_mean+1)
data_cpDB_sig_mean_filtered <- data_cpDB_sig_mean_log[rowSums(data_cpDB_sig_mean > 0) > 0, ]
ggdata <- reshape2::melt(t(data_cpDB_sig_mean_filtered))

#data_cpDB_pval_filtered <- data_cpDB_pval[,c("interacting_pair", pairs)]
data_cpDB_pval_filtered <- data_cpDB_pval
rownames(data_cpDB_pval_filtered) <- data_cpDB_pval_filtered$interacting_pair
data_cpDB_pval_filtered <- data_cpDB_pval_filtered[, -1]
data_cpDB_pval_filtered <- data_cpDB_pval_filtered[rownames(data_cpDB_sig_mean_filtered), ]
data_cpDB_pval_filtered <- -log10(data_cpDB_pval_filtered + 1e-5)
ggdata2 <- reshape2::melt(t(data_cpDB_pval_filtered))
ggdata$p_val <- ggdata2$value
ggdata$p_val[ggdata$p_val > 3] <- 3
colnames(ggdata) <- c("cluster_pair", "LR_pair", "mean", "p_val")


p <- ggplot(ggdata, mapping = aes(x = cluster_pair, y = LR_pair)) + 
  geom_point(aes(color = mean, size =  p_val)) + 
  scale_size_continuous(range=c(0,5)) +
  scale_color_gradientn(colors = c("#010103","#2b378b","#e5e132","#ea4e1a","#e1311f")) +
  theme(axis.text = element_text(face = "bold"), axis.title = element_blank(), 
        axis.text.x = element_text(angle=90, size= 12,color="darkred",hjust = 1),
        axis.text.y = element_text(angle=0, size= 12,color="darkred",hjust = 1))
p +  coord_flip()
		