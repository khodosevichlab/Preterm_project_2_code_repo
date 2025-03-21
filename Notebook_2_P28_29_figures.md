---
title: "Notebook_2_P28_29_figures"
author: "Laura Wolbeck"
date: "2023-12-03"
output: html_document
---

This script generates the Fig. and Suppl. Fig. 3 D - I.

# Setup
```{r setup}
source("../scRNA_helper.R")
library(CRMetrics)
library(magrittr)
library(dplyr)
library(qs)
library(pagoda2)
library(conos)
library(ggplot2)
library(ggrastr)
library(cowplot)
library(Seurat)
```

# 1. Read in data

```{r}
anno <- qread("anno_rough_1.qs") #annotation file
con <- qread("con.qs")
```

# Fig. 4I  
defining colours of clusters
```{r}
colours <- c("#FF3232","#3232FF", "#32FF32", "#FFFF32", "#FFC8FF", "#FFAF32","#E632E6","#AF0F0F","#320faf","#AF32FF","#ff6432", "#9cff00", "#00ffcc") %>% setNames(c( "A_cells","Astrocytes", "C_cells", "Dividing_cells", "Endothelial_cells","Ependymal_cells","Microglia" ,"Neurons","Oligodendrocytes", "OPCs", "Pericytes_VSMCs" ,"Erythrocytes" ,"ChP_epithelial_cells" ))
```

exporting final UMAP 
```{r}
plot <- con$plotGraph(groups=anno, label="") + scale_colour_manual(values=colours)+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
plot
rasterize(plot, layers='Point', dpi=1000)
ggsave("UMAP_all_cells.pdf", width=14, height=8.6)
```

# Fig. S4D
plot UMAP with gene expression (Astrocytes)
```{r}
genes <- c("Gfap", "Aqp4", "Thbs4", "Cdh5")
```
```{r, fig.height=8}
plots <- lapply(genes, function(gene) {
 con$plotGraph(gene=gene, title=gene, size=0.2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border=element_blank())
})

plot <- cowplot::plot_grid(plotlist = plots)
rasterize(plot, layers='Point', dpi=1000)
ggsave("UMAP_gene_expression.pdf", width=16, height = 12, units = "cm")
```

# Suppl. Fig. 3H

```{r}
plotClusterBarplots(con, groups = anno , show.entropy = F, show.size = F  , sample.factor = conditions) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_discrete(name = "condition")+  scale_x_discrete(limits = level_order)+
  theme(plot.margin = margin(0.5,0,0,2, "cm")) +xlab(" ")
ggsave("Barplot_condition.pdf", width=14, height = 8, units = "cm")
```

# Suppl. Fig. 3I
```{r}
plotClusterBarplots(con, groups = anno , show.entropy = F, show.size = F) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_discrete(name = "condition")+  scale_x_discrete(limits = level_order)+
  theme(plot.margin = margin(0.5,0,0,2, "cm")) +xlab(" ")
ggsave("Barplot_samples.pdf", width=15, height = 8, units = "cm")
```


# 2. Generate seurat object
generate seurat object to use seurats plotting functions

get raw counts from filtered cms
```{r}
cm <- con$getJointCountMatrix(raw=T) %>% 
  Matrix::t() %>% 
  as.sparse()
```

create Seurat object
```{r}
seurat <- CreateSeuratObject(counts = cm, project = "P_28_29", min.cells = 0, min.features = 0)
```

add normalized cm
```{r}
cm_norm <- con$getJointCountMatrix(raw=F) %>% 
  Matrix::t() %>% 
  as.sparse()
```

```{r}
seurat <- SetAssayData(object= seurat, slot= "data", cm_norm )
```

to make sure that cells in annotations are in same order as in cm
read in anno 
```{r}
# order cells in "anno" as in "cm"
anno %<>% .[match(colnames(cm), names(.))] 
```

add annotations as metadata
```{r}
seurat$annotation <- anno
```

create condition factor
```{r}
sample <- con$getDatasetPerCell()

conditions <- gsub("full_P28_1", "fullterm", sample)
conditions <- gsub("full_P28_2", "fullterm", conditions)
conditions <- gsub("pre_P29_1", "preterm", conditions)
conditions <- gsub("pre_P29_2", "preterm", conditions)
```

```{r}
conditions %<>% as.factor()
names(conditions) <- names(sample)
```

to make sure that cells in conditions are in same order as in cm
```{r}
# order cells in "conditions" as in "cm"
conditions %<>% .[match(colnames(cm), names(.))] 
```

add condition metadata to seurat object
```{r}
seurat$condition <- conditions
```
```{r}
qsave(seurat, "seurat.qs")
seurat <- qread("seurat.qs")
```

# Suppl. Fig. 3E
```{r}
genes <- c("Gfap",
"Aqp4",
"Thbs4",
"Notum",
"Ascl1",
"Egfr",
"Mki67",
"Dcx",
"Gad1",
"Rbfox3",
"Pdgfa",
"Cspg4",
"Mog",
"Hbb-bs" ,
"Hba-a2",
"Cx3cr1",
"Mpeg1",
"Flt1",
"Cdh5",
"Pdgfrb",
"Foxj1",
"Enkur",
"Ttr")
```

```{r}
Idents(seurat) <- "annotation"
```

```{r}
level_order <- c("Astrocytes",
"C_cells",
"Dividing_cells",
"A_cells",
"Neurons",
"OPCs",
"Oligodendrocytes",
"Erythrocytes",
"Microglia",
"Endothelial_cells",
"Pericytes_VSMCs",
"Ependymal_cells",
"ChP_epithelial_cells")
```

```{r}
seurat$annotation <- factor(x = seurat$annotation, levels = level_order)
DotPlot(seurat, features = rev(genes), cols = c("Grey", "Blue")) + coord_flip()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab(" ") + ylab(" ")
ggsave("Dotplot.pdf", width=16, height = 16, units = "cm")
```

# Suppl. Fig. 3F
```{r}
VlnPlot(seurat, features= "nFeature_RNA", pt.size = 0, log = F) +  ggtitle(" ") + ylab("genes per cell") + xlab(" ")+ theme_bw(base_size = 16) + theme( legend.position = "none") +  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.margin = margin(0.5,0,0,1, "cm")) 
ggsave("VlnPlot_genes_clusterwise.pdf")
```

# Suppl. Fig. 3G
```{r}
VlnPlot(seurat, features= "nCount_RNA", pt.size = 0, log = T) +  ggtitle(" ") + ylab("UMIs per cell") + xlab(" ")+ theme_bw(base_size = 16) + theme( legend.position = "none") +  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.margin = margin(0.5,0,0,1, "cm")) 
ggsave("VlnPlot_UMIs_clusterwise.pdf")
```
