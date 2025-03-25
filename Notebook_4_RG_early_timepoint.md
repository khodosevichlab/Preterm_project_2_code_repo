---
title: "Notebook_4_RG_early_timepoint"
author: "Laura Wolbeck"
date: "2023-12-03"
output: html_document
---

This script extracts the Radial glia clusters from the dataset "Early timepoints" (E 18.5 and preterm P0) from the previous paper: Kawase K, Nakamura Y, Wolbeck L, Takemura S, et al. Significance of birth in the maintenance of quiescent neural stem cells. Sci Adv. 2025 Jan 24;11(4):eadn6377. doi: 10.1126/sciadv.adn6377. Epub 2025 Jan 22. PMID: 39841848; PMCID: PMC11753423. \
Clustering and annotation of the Radial glia cells is performed.
Then, with the R package cacoa, the differentially expressed genes are estimated and Gene set enrichment analysis (GSEA) is performed.

# Setup
```{r setup, results='hide', warning=FALSE, message = FALSE}
#Load helper functions
source("../scRNA_helper.R")

#Load libraries
library(pagoda2)
library(magrittr)
library(ggplot2)
library(cacoa) #version 0.4.0
library(qs)
library(org.Mm.eg.db)
library(openxlsx)
library(ggrastr)
library(enrichplot)
```

# 1. Read in data
Read in the conos object containing early timepoints data (E 18.5 and preterm P0) and annotation from the previous paper "Significance of birth in the maintenance of quiescent neural stem cells"
```{r}
anno <- readRDS("anno.rds")
con <- readRDS("con.rds")
cms <- readRDS("cms.rds")
```

```{r}
#get RG cluster names
RG_clusters <- c("RG", "Mki67+_RG", "Intermediate_progenitor_cells")

#get RG cells
RG_cells <- names(anno[anno %in% RG_clusters])
```

```{r}
# keep only cells in cms that are in the RG cell vector
cms_RG <- cms %>% lapply(function(x) x[,colnames(x) %in% RG_cells])
saveRDS(cms_RG, "cms_RG.rds")
```

# 2. Initial conos object
```{r}
names = c("E18_5_1_",
           "E18_5_2_",
           "pre_P0_1_",
           "pre_P0_2_")
```
```{r}
con_RG <- quickConos(cms_RG,
                  names,
                  n.cores.p2=10,
                  n.cores.con=20, get.tsne = TRUE, alignment.strength=0.3)

con_RG <- con_RG$con
```
Rerun leiden clustering to find additional clusters

```{r, eval=FALSE}
con_RG$findCommunities(method = leiden.community, resolution = 3, min.group.size = 1)
```

```{r, eval=F}
leiden <- con_RG$clusters$leiden$groups %>% factor
```

```{r, eval=F}
con_RG$clusters$leiden$groups <- leiden
```

```{r}
saveRDS(con_RG, "con_RG.rds")
```

# 3. Remove cluster 7

```{r}
anno_RG <- getConosCluster(con_RG)
```
Check cells per sample
```{r}
sapply(cms_RG, function(d) dim(d)[2])
sum(sapply(cms_RG, function(d) dim(d)[2]))
```

remove "cluster 7" cells
```{r}
seven <- names(anno_RG[anno_RG %in% "7"])
cms_RG %<>% lapply(function(x) x[,!colnames(x) %in% seven])
```

Check cells per sample after removal 
```{r}
sapply(cms_RG, function(d) dim(d)[2])
sum(sapply(cms_RG, function(d) dim(d)[2]))
```
```{r}
# cms of RG clusters with cluster 7 cells removed
saveRDS(cms_RG, "cms_RG2.rds")
```

# 4. UMAP embedding and clustering

```{r}
con_RG <- quickConos(cms_RG2,
                  names,
                  n.cores.p2=10,
                  n.cores.con=20, get.tsne = F, alignment.strength= 0.3)

con_RG <- con_RG$con
```

```{r}
saveRDS(con_RG, "con_RG2.rds")
```

Final clusters were partly formed by manual selection of cells based on marker gene expression.

# 5. Create cacoa object

```{r,eval=FALSE }
condition <- setNames(c("fullterm", "fullterm", "preterm", "preterm"), names(con_RG$samples))
```
Create Cacoa object
```{r,eval=FALSE}
cao <- Cacoa$new(con_RG, cell.groups = anno_RG, sample.groups=condition, n.cores = 20, ref.level = "fullterm", target.level = "preterm")
```

# 6. DEG analysis
```{r,eval=FALSE}
cao$estimateDEPerCellType(n.cores=50, min.cell.count = 10, min.cell.frac = 0.05)
#DEs not calculated for 2 cell group(s): Immature ependymal cells-2, RG_5
```

```{r}
DEG_early_RG <- cao$test.results$de %>% lapply("[[", "res") 
```

sort ascending by padj
```{r}
DEG_early_RG %<>% lapply( function(x) x[order(x$padj),])
```

```{r}
DEG_early_RG_top100 <- lapply(DEG_early_RG, function(x) x[c(1:100),])
```

write top 100 (by padj) to excel file
```{r}
write.xlsx(x = DEG_early_RG_top100, file= "DEG_early_RG_top100.xlsx")
```

write to excel file, one sheet per cell type
```{r}
require(openxlsx)
write.xlsx(DEG_early_RG, file = "DEG_early_RG.xlsx")
```

# 7. GSEA (Table S1)
```{r}
cao$estimateOntology(type="GSEA", org.db=org.Mm.eg.db)
```

```{r}
cao$saveOntologyAsTable("early_RG_GSEA.tsv", name="GSEA")
write.xlsx(x = early_RG_GSEA, file= "early_RG_GSEA.xlsx")
```
```{r, eval=FALSE, echo=FALSE}
qsave(cao, "cao_RG_early.qs")
```

# Fig. S1D 
```{r}
#defining colours of clusters
colours <- c("#074770","#3232ff", "#AFB4DB", "#659AD2", "#cc3216","#0067C0", "#32ff32", "#ffff32", "#ffaf32") %>% setNames(c( "Lateral-ventral_RG", "Dorsal_RG", "Dorsal-lateral_RG", "Primed_active_RG", "Immature_E_cell", "Primed_active_RG_dorsal", "NPC", "Mki67_RG", "Pre_E"))

plot <- con$plotGraph(groups=anno, label="") + scale_colour_manual(values=colours)+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
plot
rasterize(plot, layers='Point', dpi=1000)
ggsave("UMAP_RG_new.pdf", width=7, height=4.3)
```

# Fig. S1E
create condition factor
```{r}
sample <- con$getDatasetPerCell()

conditions <- gsub("E18_5_1_", "fullterm", sample)
conditions <- gsub("E18_5_2_", "fullterm", conditions)
conditions <- gsub("pre_P0_1_", "preterm", conditions)
conditions <- gsub("pre_P0_2_", "preterm", conditions)

conditions %<>% as.factor()
names(conditions) <- names(sample)
```
```{r}
plotClusterBarplots(con, show.entropy = F, show.size = F  , sample.factor = conditions,groups = anno) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_discrete(name = "condition")+
 theme(plot.margin = margin(0.5,0,0,2, "cm")) +xlab(" ")
ggsave("Barplot_condition_RG.pdf", width=14, height = 8, units = "cm")
```

# Fig. S1F
```{r}
plotClusterBarplots(con, show.entropy = F, show.size = F,groups = anno) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_discrete(name = "samples")+ 
 theme(plot.margin = margin(0.5,0,0,2, "cm")) +xlab(" ")
ggsave("Barplot_samples_RG.pdf", width=14, height = 8, units = "cm")
```

# Fig. S1G
```{r}
genes <- c("Rlbp1", "Crym", "Emx1", "Ascl1")
```
```{r, fig.height=8}
plots <- lapply(genes, function(gene) {
 con$plotGraph(gene=gene, title=gene, size=0.6) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border=element_blank())
})

plot <- cowplot::plot_grid(plotlist = plots)
rasterize(plot, layers='Point', dpi=1000)
ggsave("UMAP_RG_gene_expression.pdf", width=12, height = 9, units = "cm")
```

# Fig. S1H
```{r}
enrichplot::dotplot(cao$test.results$GSEA$res$'Dorsal_RG'$BP,showCategory=c("regulation of gliogenesis",
"renal tubule morphogenesis",
"camera-type eye morphogenesis",
"eye morphogenesis",
"ureteric bud morphogenesis",
"mesonephric tubule morphogenesis",
"branching involved in ureteric bud morphogenesis",
"cardiac septum morphogenesis",
"embryonic morphogenesis",
"ear morphogenesis",
"pulmonary valve morphogenesis",
"embryonic organ morphogenesis",
"retina morphogenesis in camera-type eye",
"lens morphogenesis in camera-type eye",
"negative regulation of blood vessel morphogenesis"))+ theme(axis.text.y = element_text(size = 10))
ggsave("GSEA_Dorsal_RG.pdf", width=5.5)
```
