---
title: "Notebook_1_P28_29_preprocessing"
author: "Laura Wolbeck"
date: "2023-11-21"
output: html_document
---
This script performs the preprocessing, filtering and initial visualizations of the data using the R package CRMetrics.
Finally, the UMAP embedding and clusters are estimated also within the CRMetrics package which uses pagoda2 and conos.  

# Setup
```{r setup}
source("scRNA_helper_functions_preterm_2.R")
library(CRMetrics)
library(magrittr)
library(dplyr)
library(qs)
library(pagoda2)
library(conos)
library(ggplot2)
```

# 1. Read in data
```{r}
crm <- CRMetrics$new(data.path = "/data/Japan_preterm_counts/preterm_P28_29/HN00143770_result_10X/", 
                    metadata = "/people/nrq364/preterm-2-NCU-UCPH_figures/metadata.csv", 
                    sep.meta = ",",
                     n.cores = 50,
                     verbose = T)
```

add count matrices
```{r}
crm$addDetailedMetrics()
```

median genes per cell across all samples
```{r}
median(crm$summary.metrics[crm$summary.metrics$metric=="median genes per cell",]$value)
```

cells per samples and sum of all cells
```{r}
sapply(crm$cms, function(d) dim(d)[2])
sum(sapply(crm$cms, function(d) dim(d)[2]))
```

## Suppl. Fig. 3A
plot sequenced cells per sample
```{r}
crm$plotSummaryMetrics(comp.group = "sample", second.comp.group = "condition", metrics = "estimated number of cells",plot.geom = "bar") + xlab("") + ylab("number of cells") + theme(legend.position = "right", legend.title = element_blank())
ggsave("number_of_cells.pdf", width=10, height = 8, units = "cm")
```

## Suppl. Fig. 3B
plot genes per cell
```{r}
crm$plotDetailedMetrics(comp.group = "condition",
                        metrics = "gene_count", 
                        plot.geom = "violin",hline = F)  + ylab("genes per cell") + theme(legend.position = "right", legend.title = element_blank()) +  stat_summary(fun.y=median, size=8, geom = "text", label = "-")  #+ theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))
ggsave("VlnPlot_genes_per_cell.pdf", width=10, height = 8, units = "cm")
```

# 2. Intitial conos object
before filtering
```{r}
crm$doPreprocessing(preprocess = "pagoda")
crm$createEmbedding()
```

# 3. Filtering
## 3.1 Scrublet to estimate doublets
in terminal run (or in terminal tab in R studio): 
module load miniconda/4.12.0 
conda activate doublets
```{r}
crm$detectDoublets(env = "doublets", conda.path = "/opt/software/miniconda/4.12.0/condabin/conda", method = "scrublet") 
#Detected 925 possible doublets out of 29127 cells 
```

## Suppl. Fig. 3C
```{r}
plot <- crm$plotEmbedding(doublet.method = "scrublet", size = 0.2)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
rasterize(plot, layers='Point', dpi=1000)
ggsave("Doublets_UMAP.pdf", width=7, height=4.3 )
```

## 3.2 Filter count matrices for depth, doublets and mito fraction

```{r}
crm$filterCms(depth.cutoff = 750, 
              mito.cutoff = 0.07, 
              doublets = "scrublet",
              species = "mouse")
#Removing 3231 of 29127 cells (11.1%)
```

save filtered count matrix
```{r}
qsave(crm$cms.filtered, "cms_P28_29_filtered.qs", 
      nthreads = 10)
```

# 4. UMAP embedding and clustering
pre-processing of each filtered dataset with pagoda2, followed by using conos to build the joint UMAP graph with forced alignment (alignment.strength = 0.3) to integrate samples well
```{r}
cms <- qread("cms_P28_29_filtered.qs") 
```

```{r}
#vector with sample names 
names = c("full_P28_1",
           "full_P28_2",
           "pre_P29_1",
           "pre_P29_2")
```
```{r, eval=FALSE}
con <- quickConos(cms,
                  names,
                  n.cores.p2=10,
                  n.cores.con=20, get.tsne = F, alignment.strength = 0.3) # 0 shows very strong batch by sample, 0.2 also

con <- con$con
```

Rerun leiden clustering to find additional clusters
```{r, eval=FALSE}
con$findCommunities(method = leiden.community, resolution = 5, min.group.size = 15)
```

save clusters and conos object
```{r}
leiden22<- con$clusters$leiden$groups %>% factor
qsave(leiden22, "leiden22.qs")
qsave(con, "con.qs")
```
For annotation, some of the clusters were collapsed or manually separated.



