---
title: "Notebook_5_cell_adhesion_RG_Endo"
author: "Laura Wolbeck"
date: "2023-12-03"
output: html_document
---

This script uses the fullterm P2 data from the previous paper: Kawase K, Nakamura Y, Wolbeck L, Takemura S, et al. Significance of birth in the maintenance of quiescent neural stem cells. Sci Adv. 2025 Jan 24;11(4):eadn6377. doi: 10.1126/sciadv.adn6377. Epub 2025 Jan 22. PMID: 39841848; PMCID: PMC11753423.

# Setup
```{r setup, results='hide', warning=FALSE, message = FALSE}
library(dplyr)
library(magrittr)
library(openxlsx)
```

# 1. Read in data
Read in conos object and annotation from "late timepoints" dataset of the previous paper.
```{r}
con <- readRDS("con2_m5_2.rds")
anno <- readRDS("anno2_final6_rough.rds")
```

# 2. Get list of genes of the GO term "cell adhesion"
```{r}
go_term <- "GO:0007155"  # cell adhesion
```

Extract all genes associated with the GO term GO:0007155 cell adhesion
```{r}
gene_list <- AnnotationDbi::select(org.Mm.eg.db, 
                                   keys = "GO:0007155", 
                                   columns = c("SYMBOL", "GENENAME"), 
                                   keytype = "GO")
```

Extract the gene symbols as a vector
```{r}
gene_vector <- unique(gene_list$SYMBOL)
```

# Table S2
Expression of cell adhesion-related genes in RG and Endothelial cell/Pericytes/ Vascular smooth muscle cells clusters of term P2.5 samples.

Get joint count matrix (normalized) from conos object
```{r}
cm <- con$getJointCountMatrix(raw=FALSE)
```
```{r}
cells <- names(anno[anno %in% c("RG", "Endothelial/Pericytes/VSMC")])
cells_fullterm <- grep("full_P2", cells, value=T) # 1262 cells
```

Filter count matrix for selected cells 
```{r}
cm_RG <- cm[cells_fullterm, ]
```

Check detected genes
```{r}
genes <- colnames(cm)
```

Identify missing genes from GO term list
```{r}
missing_genes <- setdiff(gene_vector, genes)
cat("Number of missing genes:", length(missing_genes), "\n") # 69
```

Retain only detected GO term genes
```{r}
genes_final <- intersect(gene_vector, genes)
```

Filter count matrix for selected genes
```{r}
cm_RG <- cm_RG[, genes_final]
expression_data <- as.data.frame(cm_RG)
```

Extract annotations
```{r}
anno_RG <- anno[anno %in% c("RG", "Endothelial/Pericytes/VSMC")]
anno_RG %<>% droplevels()
anno_RG <- as.data.frame(anno_RG)
```

Merge expression data with annotations
```{r}
data_merged <- merge(expression_data, anno_RG, by = "row.names")
```
Summarize expression per cell type
```{r}
expression_celltype <- group_by(data_merged, anno_RG)
summary_df <- summarise(expression_celltype, across(-1, mean))
summary_df %<>% t() %>% as.data.frame()
colnames(summary_df) <- c("RG", "Endothelial/Pericytes/VSMC")
summary_df <- summary_df[-1, ]
```

Remove genes with zero expression in both cell types
```{r}
summary_df <- summary_df[rowSums(summary_df != 0) > 0, ]
```

Add gene names and select relevant columns
```{r}
summary_df %<>% mutate(gene = rownames(summary_df)) %>% select("gene", "RG", "Endothelial/Pericytes/VSMC")
```

Save results as Excel file
```{r}
write.xlsx(summary_df, file = "cell_adhesion_genes_termP2.xlsx")
```

# Table S3
Differentially expressed genes of cell adhesion-related genes in RG and Endothelial cells/Pericytes/ Vascular smooth muscle cells clusters of E18.5 and term P2.5 samples.

Read in cacoa object of fullterm samples which contains DEG analysis results.
```{r}
cao <- readRDS("cao_fullterm_DEG2.rds")
```

Get DEG of RG that are also in GO term cell adhesion. 
```{r}
RG_df <- cao$test.results$de$RG$res %>% arrange(padj)
RG_DEGs <- RG_df %>% filter(padj < 0.05)  # 1122 genes
RG_DEGs_vector <- rownames(RG_DEGs)

RG_adhesion_DEGs <- RG_DEGs[intersect(gene_vector, RG_DEGs_vector), ] # 29 genes
write.xlsx(RG_adhesion_DEGs, file = "RG_adhesion_DEGs.xlsx")
```

Get DEG of Endo cells that are also in GO term cell adhesion. 
```{r}
Endo_df <- cao$test.results$de$'Endo/Peri'$res %>% arrange(padj)
Endo_DEGs <- Endo_df %>% filter(padj < 0.05)  # 51 genes
Endo_DEGs_vector <- rownames(Endo_DEGs)

Endo_adhesion_DEGs <- Endo_DEGs[intersect(gene_vector, Endo_DEGs_vector), ] # 2 genes
write.xlsx(Endo_adhesion_DEGs, file = "Endo_adhesion_DEGs.xlsx")
```

# Table S6
Differentially expressed genes of cell adhesion-related genes in RG and Endothelial cell/Pericytes/ Vascular smooth muscle cells clusters of term P2.5 and preterm P3.5 samples. 

Read in cacoa object of fullterm samples which contains DEG analysis results.
```{r}
cao <- readRDS("cao_fullterm_DEG2.rds")
```

Get DEG of RG that are also in GO term cell adhesion. 
```{r}
RG_df <- cao$test.results$de$RG$res %>% arrange(padj)
RG_DEGs <- RG_df %>% filter(padj < 0.05)  # 109 genes
RG_DEGs_vector <- rownames(RG_DEGs)

RG_adhesion_DEGs <- RG_DEGs[intersect(gene_vector, RG_DEGs_vector), ] # 5 genes
write.xlsx(RG_adhesion_DEGs, file = "RG_late_adhesion_DEGs.xlsx")
```

Get DEG of Endo cells that are also in GO term cell adhesion. 
```{r}
Endo_df <- cao$test.results$de$'Endo/Peri'$res %>% arrange(padj)
Endo_DEGs <- Endo_df %>% filter(padj < 0.05)  # 12 genes
Endo_DEGs_vector <- rownames(Endo_DEGs)

Endo_adhesion_DEGs <- Endo_DEGs[intersect(gene_vector, Endo_DEGs_vector), ] # 2 genes
write.xlsx(Endo_adhesion_DEGs, file = "Endo_late_adhesion_DEGs.xlsx")
```
