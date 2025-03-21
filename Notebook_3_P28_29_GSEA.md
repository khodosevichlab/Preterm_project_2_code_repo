---
title: "Notebook_3_P28_29_GSEA"
author: "Laura Wolbeck"
date: "2023-12-03"
output: html_document
---

This script covers the gene set enrichment analysis (GSEA) and uses the R package cacoa.

# Setup
```{r}
library(magrittr)
library(ggplot2)
library(cacoa) #version 0.4.0
library(qs)
library(openxlsx)
library(org.Mm.eg.db)
library(dplyr)
library(enrichplot)
```

# 1. Read in data
```{r}
con <- qread("con.qs")
cms <- qread("cms_P28_29_filtered.qs") 
anno <- qread("anno_rough_1.qs")
```

# 2. Create cacoa object
Create condition factor
```{r,eval=FALSE }
condition <- setNames(c("fullterm", "fullterm", "preterm", "preterm"), names(con$samples))
```

Create Cacoa object
```{r,eval=FALSE}
cao <- Cacoa$new(con, cell.groups = anno, sample.groups=condition, n.cores = 10, ref.level = "fullterm", target.level = "preterm")
```
```{r}
qsave(cao_1, "cao_1.qs")
```

# 3. DEG analysis

estimate DE genes
```{r,eval=FALSE}
cao_1$estimateDEPerCellType(n.cores=50) 
#DEs not calculated for 1 cell group(s): ChP_epithelial_cells
```

save DEG as excel file 
```{r}
DEG_1 <- cao_1$test.results$de %>% lapply("[[", "res") 
write.xlsx(x = DEG_1, file= "DEG_1.xlsx")
```

# 4. GSEA
```{r}
cao_1$estimateOntology(type="GSEA", org.db=org.Mm.eg.db)
```

```{r}
qsave(cao_1, "cao_1.qs")
```

save GSEA results of RG as excel file
```{r}
cao_1$saveOntologyAsTable("GSEA_1.tsv", name="GSEA")
GSEA_1 <- read.table("GSEA_1.tsv", sep="\t")
write.xlsx(x = GSEA_1, file= "GSEA_1.xlsx")
```

# Fig. 4J
```{r, fig.width=8, fig.heigth=20}
enrichplot::dotplot(cao_1$test.results$GSEA$res$'B/As_cells'$BP,showCategory=c("synaptic membrane adhesion","homotypic cell-cell adhesion",
"homophilic cell adhesion via plasma membrane adhesion molecules",
"positive regulation of cell-matrix adhesion",
"cell-cell adhesion via plasma-membrane adhesion molecules",
"regulation of cell-matrix adhesion",
"positive regulation of cell junction assembly",
"cell-matrix adhesion",
"regulation of cell-substrate adhesion",
"cell-substrate adhesion",
"regulation of cell junction assembly",
"cell adhesion",
"cell-cell adhesion",
"cell junction assembly",
"cell junction organization",
"regulation of cell adhesion",
"positive regulation of endothelial cell migration",
"regulation of endothelial cell migration",
"angiogenesis")) + theme(axis.text.y = element_text(size = 8))
ggsave("GSEA_Astrocytes.pdf", width = 17, height = 17, units = "cm")
```

# Fig. 4K
```{r}
enrichplot::dotplot(cao_1$test.results$GSEA$res$Endothelial_cells$BP,showCategory=c("homotypic cell-cell adhesion" , "heterophilic cell-cell adhesion via plasma membrane cell adhesion molecules"))
ggsave("GSEA_Endothelial_cells.pdf", width = 17, height = 11, units = "cm")
```

# Fig. S4J
```{r}
#define words to exclude from collapsed GO terms
ex_words <- c('regulation', 'process', 'cell')
As <- cao_1$plotOntologyHeatmapCollapsed(
  name="GSEA", genes="down", n=59, size.range=c(1, 4), subtype="BP", exclude.words=ex_words, clust.method="ward.D")

#filter for Astrocytes
As$data %<>% filter(G2=='B/As_cells', value > 0)

df <- As$data
df %<>% arrange(value) %>% mutate(G1 = G1 %>% factor(., levels = .))
```

```{r}
ggplot(df, aes(x=G2, y = G1, fill = value)) + geom_tile( colour="white", size=0.2) +  scale_fill_distiller(direction = 1 , palette = "Blues") + guides(fill = color.guide) + labs(title="Top 30 down (Astrocytes)", y = " ", x= " ") + theme_grey(base_size = 11) + theme(axis.ticks.x = element_blank(),  axis.text.x = element_blank())
ggsave("GO_As_down.pdf", width = 16, height = 11, units = "cm")
```

# Fig. S4K
```{r}
#define color legend
color.guide <- guide_colorbar(title = "-log10(p-value)", title.position = "left", 
        title.theme = element_text(angle = 90, hjust = 0.5))
```

```{r}
#define words to exclude from collapsed GO terms
ex_words <- c('regulation', 'process', 'cell')
As <- cao_1$plotOntologyHeatmapCollapsed(
  name="GSEA", genes="up", n=50, size.range=c(1, 4), subtype="BP", exclude.words=ex_words, clust.method="ward.D")

#filter for Astrocytes
As$data %<>% filter(G2=='B/As_cells', value > 0)

df <- As$data
df %<>% arrange(value) %>% mutate(G1 = G1 %>% factor(., levels = .))
```

```{r}
ggplot(df, aes(x=G2, y = G1, fill = value)) + geom_tile(colour="white", size=0.2) +  scale_fill_distiller(direction = 1 , palette = "Reds") + guides(fill = color.guide) + labs(title="Top 30 up (Astrocytes)", y = " ", x= " ") + theme_grey(base_size = 11) + theme(axis.ticks.x = element_blank(),  axis.text.x = element_blank())
ggsave("GO_As_up.pdf", width = 15, height = 11, units = "cm")
```

# Fig. S4L
```{r}
#define words to exclude from collapsed GO terms
ex_words <- c('regulation', 'process', 'cell')
Endo <- cao_1$plotOntologyHeatmapCollapsed(
  name="GSEA", genes="down", n=140, size.range=c(1, 4), subtype="BP", exclude.words=ex_words, clust.method="ward.D")

#filter for Endothelial cells
Endo$data %<>% filter(G2=='Endothelial_cells', value > 0)

df <- Endo$data
df %<>% arrange(value) %>% mutate(G1 = G1 %>% factor(., levels = .))
```

```{r}
ggplot(df, aes(x=G2, y = G1, fill = value)) + geom_tile( colour="white", size=0.2) +  scale_fill_distiller(direction = 1 , palette = "Blues") + guides(fill = color.guide) + labs(title="Top 30 down (Endothelial cells)", y = " ", x= " ") + theme_grey(base_size = 11) + theme(axis.ticks.x = element_blank(),  axis.text.x = element_blank())
ggsave("GO_Endo_down.pdf", width = 15, height = 11, units = "cm")
```

# Fig. S4M
```{r}
Endo <- cao_1$plotOntology(
  name="GSEA", genes="up", subtype="BP",  cell.type = 'Endothelial_cells' )

df <- Endo$data
df <- df[,c("Description", "p.adjust")]
df %<>% mutate(value = -log10(p.adjust), G2=c("Endothelial_cell", "Endothelial_cell", "Endothelial_cell", "Endothelial_cell", "Endothelial_cell", "Endothelial_cell", "Endothelial_cell", "Endothelial_cell", "Endothelial_cell"))
df %<>% arrange(value) %>% mutate(Description = Description %>% factor(., levels = .))
```

```{r}
ggplot(df, aes(x=G2, y = Description, fill = value)) + geom_tile( colour="white", size=0.2) +  scale_fill_distiller(direction = 1 , palette = "Reds") + guides(fill = color.guide) + labs(title="9 up (Endothelial cells)", y = " ", x= " ") + theme_grey(base_size = 11) + theme(axis.ticks.x = element_blank(),  axis.text.x = element_blank())
ggsave("GO_Endo_up.pdf", width = 14, height = 5, units = "cm") 
```
