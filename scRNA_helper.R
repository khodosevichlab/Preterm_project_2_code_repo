#Rename column names, making them unique
renameCols <- function(cms, names) {
  if(!length(cms)==length(names)) stop("Names must match number of count matrices.")
  
  mapply(function(c, n) {
    colnames(c) <- lapply(colnames(c), function(cn) paste0(n,cn)) %>% unlist
    return(c)
  }, c=cms, n=names, SIMPLIFY = F)
}

#Rename column/cell names in velocyto objects
renameVeloCells <- function(cms.list, pattern, replacement) {
  cnames.corr <- lapply(cms.list, colnames) %>% lapply(function(c) gsub(pattern, replacement, c))
  cms.list %<>% mapply(function(corr,cms) {
    colnames(cms) <- corr
    return(cms)
  }, corr = cnames.corr,
  cms = .)
  return(cms.list)
}

#Prepare velocyto data with Conos object
veloOnConos <- function(list.velo, filter.limit, con=con, clustering=NULL, groups=NULL, n.odgenes=2e3, verbose=TRUE, min.max.cluster.average.emat=0.5, min.max.cluster.average.nmat=0.05, min.max.cluster.average.smat=0.01, ncomps=100) {
  if(verbose) cat("Filter Velocyto object ")
  genes <- sapply(con$samples, conos:::getGenes) %>% unlist %>% unique
  cells <- sapply(con$samples, conos:::getCellNames) %>% unlist %>% unique
  
  list.velo %<>% lapply(conos:::prepareVelocity, genes=genes, cells=cells)
  if(verbose) cat(".")
  
  emat <- do.call(cbind, lapply(list.velo, function(x) {x[[1]]}))
  emat %<>% .[,colSums(.)>=filter.limit]
  rownames(emat) <- make.unique(rownames(emat))
  if(verbose) cat(".")
  
  nmat <- do.call(cbind, lapply(list.velo, function(x) {x[[2]]}))
  nmat %<>% .[,colnames(.) %in% colnames(emat)]
  rownames(nmat) <- rownames(emat)
  if(verbose) cat(".")
  
  smat <- do.call(cbind, lapply(list.velo, function(x) {x[[3]]}))
  smat %<>% .[,colnames(.) %in% colnames(emat)]
  rownames(smat) <- rownames(emat)
  if(verbose) cat(".")
  
  groups <- conos:::parseCellGroups(con, clustering, groups) %>% .[names(.) %in% colnames(emat)]
  if(any(table(groups)<2)) warning(" groups with less than 2 cells detected in 'groups'/'clustering', these are excluded ...")
  groups %<>% .[. %in% names(table(.)[table(.)>1])]
  
  emat %<>% .[,colnames(.) %in% names(groups)] %>% velocyto.R::filter.genes.by.cluster.expression(groups, min.max.cluster.average=min.max.cluster.average.emat)
  nmat %<>% .[,colnames(.) %in% names(groups)] %>% velocyto.R::filter.genes.by.cluster.expression(groups, min.max.cluster.average=min.max.cluster.average.nmat)
  smat %<>% .[,colnames(.) %in% names(groups)] %>% velocyto.R::filter.genes.by.cluster.expression(groups, min.max.cluster.average=min.max.cluster.average.smat)
  if(verbose) cat(" done!\n")
  
  cell.colors <- pagoda2:::fac2col(groups)
  emb <- con$embedding %>% .[rownames(.) %in% colnames(emat),]
  
  emat <- emat[,order(match(colnames(emat), rownames(emb)))]
  nmat <- nmat[,order(match(colnames(nmat), rownames(emb)))]
  smat <- smat[,order(match(colnames(smat), rownames(emb)))]
  
  pcs <- conos:::pcaFromConos(con$samples, n.odgenes=n.odgenes, ncomps = ncomps) %>% .[rownames(.) %in% colnames(emat),]
  pcs <- pcs[order(match(rownames(pcs), rownames(emb))),]
  
  if (verbose) cat("Calculating cell distances...\n")
  cell.dist <- as.dist(1 - velocyto.R::armaCor(t(pcs)))
  
  if (verbose) cat("All Done!")
  return(list(cell.dist=cell.dist, emat=emat, nmat=nmat, smat=smat, cell.colors=cell.colors, emb=emb))
}

#Mitochondrial fraction
mitoFraction <- function(con, species="human") {
  if(species=="human") lapply(con$samples, function(d) Matrix::rowSums(d$counts[,grep("MT-", colnames(d$counts))]) / Matrix::rowSums(d$counts)) %>% Reduce(c, .)
  else if(species=="mouse") lapply(con$samples, function(d) Matrix::rowSums(d$counts[,grep("mt-", colnames(d$counts))]) / Matrix::rowSums(d$counts)) %>% Reduce(c, .)
  else stop("Species must either be 'human' or 'mouse'.")
}

#Shannon entropy
ShannonEntropy <- function(con, levels=20, verbose=T) {
  shannon.entropy <- function(p)  
    {if (min(p) < 0 || sum(p) <= 0)
      return(NA)
    p.norm <- p[p>0]/sum(p)
    -sum(log2(p.norm)*p.norm)
  }
  
  #Calculate
  if(verbose) message("Merging matrices...")
  rc <- conos:::mergeCountMatrices(lapply(con$samples, conos:::getRawCountMatrix))
  
  if(verbose) {
    message("Calculating entropy...")
    entropy <- unlist(pbmcapply::pbmclapply(1:dim(rc)[2], function(x) shannon.entropy(rc[,x]), mc.cores=con$n.cores))
    } else {
      entropy <- unlist(parallel::mclapply(1:dim(rc)[2], function(x) shannon.entropy(rc[,x]), mc.cores=con$n.cores))
    }
  
  
  if(verbose) message("Normalizing based on levels...")
  names(entropy) <- colnames(rc)
  
  min <- min(entropy)
  entropynorm <- floor((entropy - min)/(max(entropy) - min) * levels)
  return(entropynorm)
}

#Proportions of Conos object
conProp <- function(con, clusters, ctrl, disease) {
  dpc <- con$getDatasetPerCell()
  cond <- as.factor(setNames(ifelse(grepl(ctrl,dpc),ctrl,disease), names(dpc)))
  
  samplenames <- unique(dpc)
  
  c <- clusters[names(clusters) %in% names(dpc[cond==ctrl])]
  c_samples <- lapply(samplenames[grep(ctrl,samplenames)], function(s) {
    sprop <- c[grep(s,names(c))]
    sprop <- table(sprop)/length(sprop)
  })
  names(c_samples) <- samplenames[grep(ctrl,samplenames)]
  c <- table(c)/length(c)
  
  d <- clusters[names(clusters) %in% names(dpc[cond==disease])]
  d_samples <- lapply(samplenames[grep(disease,samplenames)], function(s) {
    sprop <- d[grep(s,names(d))]
    sprop <- table(sprop)/length(sprop)
  })
  names(d_samples) <- samplenames[grep(disease,samplenames)]
  d <- table(d)/length(d)
  
  res <- c(list(c),c_samples,list(d),d_samples) %>% setNames(c(ctrl,lapply(samplenames[grep(ctrl,samplenames)], function(n) paste0(ctrl,"_",n)) %>% unlist,disease,lapply(samplenames[grep(disease,samplenames)], function(n) paste0(disease,"_",n))))
  res <- do.call("rbind", res)
  
  return(res)
}

addEmbeddingP2Web <- function(p2, con, embedding=NULL, name="UMAP") {
  if(is.null(embedding)) embedding <- con$embedding
  
  if(identical(dim(p2$originalP2object$embeddings$PCA[[1]]),dim(embedding))) {
    p2$originalP2object$embeddings$PCA[[name]] <- embedding
    return(p2)
  } else {
    stop("The embedding dimensions of the p2.web object and the input object are not identical.")
  }
}

embedUMAP <- function(con,
                      min.dist=0.01,
                      spread=15,
                      min.prob.lower=1e-7,
                      method=leiden.community,
                      resolution=1,
                      min.group.size=50) {
  message("Creating UMAP embedding...")
  con$embedGraph(method="UMAP", 
                 min.dist=min.dist, 
                 spread=spread,
                 min.prob.lower=min.prob.lower)
  
  message("Estimating clusters...")
  con$findCommunities(method=leiden.community, resolution=resolution, min.group.size=min.group.size)
  
  return(con)
}

buildConosGraph <- function(con,
                          k.conos=15, 
                          k.self=15, 
                          space='PCA', 
                          ncomps=40,
                          n.odgenes=2e3,
                          matching.method='mNN', 
                          metric='angular', 
                          score.component.variance=T,
                          alignment.strength=0,
                          min.dist=0.01, 
                          spread=15,
                          min.prob.lower=1e-3,
                          resolution=1,
                          min.group.size=50) {
  message("Building graph...")
  con$buildGraph(k=k.conos, 
                 k.self=k.self, 
                 space=space, 
                 ncomps=ncomps, 
                 n.odgenes=n.odgenes, 
                 matching.method=matching.method, 
                 metric=metric, 
                 verbose=T, 
                 score.component.variance=score.component.variance,
                 alignment.strength=alignment.strength)
  
  embedUMAP(con=con,
            min.dist=min.dist,
            spread=spread,
            min.prob.lower=min.prob.lower,
            method=leiden.community,
            resolution=resolution,
            min.group.size=min.group.size)
  
  return(con)
}

quickConos <- function(cms, 
                       sample.names,
                       n.cores.p2,
                       n.cores.con,
                       n.odgenes=3e3, 
                       nPcs = 50, 
                       k.p2 = 30, 
                       perplexity = 50, 
                       log.scale = TRUE, 
                       trim = 10, 
                       keep.genes = NULL, 
                       min.cells.per.gene = 3, 
                       min.transcripts.per.cell = 200, 
                       get.largevis = F, 
                       get.tsne = F, 
                       make.geneknn = TRUE,
                       k.conos=15, 
                       k.self=30, 
                       space='PCA', 
                       ncomps=40, 
                       matching.method='mNN', 
                       metric='angular', 
                       score.component.variance=T,
                       alignment.strength=0,
                       min.dist=0.01, 
                       spread=15) {
  if(length(cms)==length(sample.names)) {
    message("Performing P2 processing...")
    panel.preprocessed <- lapply(cms, function(x) basicP2proc(x, n.cores = n.cores.p2,
                                                              n.odgenes = n.odgenes, 
                                                              nPcs = nPcs,
                                                              k = k.p2, 
                                                              perplexity = perplexity, 
                                                              log.scale = log.scale, 
                                                              trim = trim, 
                                                              keep.genes = keep.genes, 
                                                              min.cells.per.gene = min.cells.per.gene, 
                                                              min.transcripts.per.cell = min.transcripts.per.cell, 
                                                              get.largevis = get.largevis, 
                                                              get.tsne = get.tsne, 
                                                              make.geneknn = make.geneknn))
    
    names(panel.preprocessed) = sample.names
    con <- Conos$new(panel.preprocessed, n.cores=n.cores.con)
    
    con <- buildConosGraph(con=con,
                           k.conos=k.conos, 
                           k.self=k.self, 
                           space=space, 
                           ncomps=ncomps, 
                           n.odgenes=n.odgenes, 
                           matching.method=matching.method, 
                           metric=metric, 
                           score.component.variance=score.component.variance,
                           alignment.strength=alignment.strength,
                           min.dist=min.dist, 
                           spread=spread)
    
    return(list(con=con, panel.preprocessed=panel.preprocessed))
  } else {
    stop("Sample names must match number of count matrices.")
  }
  
}

collapseAnnotation <- function(anno, label) {
  anno %<>% factor
  idx <- grepl(label,levels(anno))
  cat(paste0("Collapsing ",sum(idx)," labels containing '",label,"' in their name into one label.\n"))
  levels(anno)[idx] <- c(label)
  anno %<>% factor
  return(anno)
}

getConosDepth <- function(con) {
  lapply(con$samples, function(d) d$depth) %>% unlist %>% setNames(.,(strsplit(names(.), ".", T) %>% 
                                                                                 sapply(function(d) d[2])))
}

getConosCluster <- function(con, name="leiden") {
  con$clusters[[name]]$groups
}

plotDotMap <- function (markers, 
                        count.matrix, 
                        annotation, 
                        marker.colour="black",
                        cluster.colour="black",
                        text.angle = 45, 
                        gene.order = NULL, 
                        cols = c("blue", "red"),
                        col.min = -2.5,
                        col.max = 2.5,
                        dot.min = 0,
                        dot.scale = 6,
                        scale.by = "radius",
                        scale.min = NA,
                        scale.max = NA,
                        verbose=T) {
  scale.func <- switch(scale.by, 'size' = scale_size, 'radius' = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  if(verbose) cat("Plotting .")
  
  if(!is.character(markers)) stop("'markers' must be a character vector.")
  
  missing.markers <- setdiff(markers, colnames(count.matrix))
  if(length(missing.markers)>0) {
  cat("Not all markers are in 'count.matrix'. The following are missing:\n",paste(missing.markers, collapse=" "),"\n")
    stop("Please update 'markers'.")
  }
  
  # From CellAnnotatoR:::plotExpressionViolinMap, should be exchanged with generic function
  p.df <- lapply(markers, function(g) data.frame(Expr = count.matrix[names(annotation), g], Type = annotation, Gene = g)) %>% Reduce(rbind, .)
  if (is.logical(gene.order) && gene.order) {
    gene.order <- unique(markers)
  } else {
    gene.order <- NULL
  }
  
  if (!is.null(gene.order)) {
    p.df %<>% dplyr::mutate(Gene = factor(as.character(Gene), 
                                          levels = gene.order))
  }
  
  # Adapted from Seurat:::DotPlot
  if(verbose) cat(".")
  data.plot <- levels(annotation) %>% lapply(function(t) {
    markers %>% lapply(function(g) {
      df <- p.df %>% filter(Type==t, Gene==g)
      pct.exp <- sum(df$Expr>0)/dim(df)[1]*100
      avg.exp <- mean(df$Expr[df$Expr>0])
      res <- data.frame(gene=g,
                        pct.exp=pct.exp,
                        avg.exp=avg.exp)
      return(res)
    }) %>% Reduce(rbind, .)
  }) %>% 
    setNames(levels(annotation)) %>%
    bind_rows(., .id="cluster")
  
  data.plot$cluster %<>% factor(., levels=rev(unique(.)))
  
  data.plot %<>% arrange(gene)
  
  data.plot$avg.exp.scaled <- data.plot$gene %>% unique %>% sapply(function(g) {
    data.plot %>% .[.$gene == g, 'avg.exp'] %>% 
      scale %>% 
      MinMax(min = col.min, max = col.max)
  }) %>% unlist %>% as.numeric
  
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  
  cluster.colour %<>% rev
  
  plot <- ggplot(data.plot, aes_string("gene", "cluster")) +
    geom_point(aes_string(size = "pct.exp", color = "avg.exp.scaled")) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.text.x = element_text(angle=text.angle, hjust = 1, colour=marker.colour),
          axis.text.y = element_text(colour=cluster.colour),
          panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(size = guide_legend(title = 'Percent expressed'), color = guide_colorbar(title = 'Average expression')) +
    labs(x = 'Marker', y = 'Cluster') +
    scale_color_gradient(low = cols[1], high = cols[2])
  if(verbose) cat(" done!")
  return(plot)
}

clusterPlots <- function(genes, cluster, annotation, limits) {
  plot1 <- list(con$plotGraph(groups=annotation, subgroups=cluster, plot.na=F, size=0.1, alpha=0.3, shuffle.colors=F) + limits) %>%
    append(list(con$plotGraph(groups=ifelse(grepl(cluster,annotation),cluster,"Other") %>% setNames(names(annotation)), plot.na=F, size=0.01, alpha=0.1, shuffle.colors=F) + limits)) %>%
    append(list(con$plotGraph(groups=group.per.cell[names(group.per.cell) %in% names(annotation[annotation==cluster])], plot.na=F, size=0.01, alpha=0.1, shuffle.colors=F) + limits)) %>%
    cowplot::plot_grid(plotlist=., ncol=2)
  
  if(length(genes)>0) {
    plot2 <- genes %>%
      lapply(function(g) con$plotGraph(groups=annotation, gene=g, alpha=0.1, size=0.1, title=g, plot.na=F) + limits) %>%
      cowplot::plot_grid(plotlist=., ncol=2)
    
    plot3 <- genes %>% lapply(function(g) {
      con$plotGraph(groups=annotation, subgroups=cluster, gene=g, plot.na=F, size=1, alpha=0.3, shuffle.colors=T, title=g) + limits
    }) %>%
      cowplot::plot_grid(plotlist=., ncol=2)
    
    plot4 <- sccore:::dotPlot(genes, cluster.cms, cluster.per.cell, n.cores=1)
    
    plot5 <- sccore:::dotPlot(genes, cluster.cms, annotation, n.cores=1)
    
    plot6 <- sccore:::dotPlot(genes, cluster.cms, group.per.cell[names(group.per.cell) %in% names(annotation[annotation==cluster])], n.cores=1)
  }
  
  df <- table(sample.per.cell[names(sample.per.cell) %in% names(annotation[annotation==cluster])]) %>%
    data.frame %>% 
    setNames(c("sample","percent"))
  df$percent <- df$percent/table(sample.per.cell[names(sample.per.cell) %in% names(annotation)]) %>% as.numeric*100
  df$group <- c(rep("CTRL",10),
                rep("MSA",9),
                rep("PD",13))
  
  plot7 <- ggplot(df, aes(sample, percent, fill=group)) + geom_col() + theme(legend.position="none") + xlab("")+ ylab("% cells of cell type per sample") + theme(axis.text.x = element_text(angle=90))
  
  plot8 <- ggplot(df, aes(group, percent, fill=group)) + geom_boxplot() + theme(legend.position="none") + xlab("")+ ylab("% cells of cell type per sample") + theme(axis.text.x = element_text(angle=90))
  
  if(length(genes)>0) return(list(plot1,plot2,plot3,plot4,plot5,plot6,plot7,plot8)) else return(list(plot1,plot7,plot8))
}

typePlots <- function(markers, annotation, limits) {
  plot1 <- con$plotGraph(groups=annotation, plot.na=F, size=0.5, alpha=0.1, shuffle.colors=T, mark.groups=F, show.legend=T, legend.position="bottom") + limits
  plot2 <- con$plotGraph(groups=annotation, plot.na=F, size=0.5, alpha=0.1, shuffle.colors=T) + limits
  
  plot3 <- cowplot::plot_grid(plotlist=lapply(c("CTRL","PD","MSA"), function(p) {
    con$plotGraph(groups=group.per.cell[names(group.per.cell) %in% names(annotation)], subgroups=p, plot.na=F, size=0.1, alpha=0.3, shuffle.colors=T, title=p) + limits
  }), ncol=2)
  
  
  plot4 <- markers %>%
    lapply(function(g) gg <- con$plotGraph(groups=group.per.cell[names(group.per.cell) %in% names(annotation)], gene=g, alpha=0.1, size=0.1, title=g, plot.na=F) + limits) %>%
    cowplot::plot_grid(plotlist=., ncol=2)
  
  plot5 <- markers %>%
    sccore:::dotPlot(cluster.cms, annotation, n.cores=1)
  
  df <- table(sample.per.cell[names(group.per.cell)%in% names(annotation)]) %>%
    data.frame %>%
    setNames(c("sample","percent"))
  df$percent <- df$percent/table(sample.per.cell) %>% as.numeric*100
  df$group <- c(rep("CTRL",10),
                rep("MSA",9),
                rep("PD",13))
  
  plot6 <- ggplot(df, aes(sample, percent, fill=group)) + geom_col() + theme(legend.position="none") + xlab("")+ ylab("% cells per sample") + theme(axis.text.x = element_text(angle=90))
  
  plot7 <- ggplot(df, aes(group, percent, fill=group)) + geom_boxplot() + theme(legend.position="none") + xlab("")+ ylab("% cells per sample")
  
  return(list(plot1,plot2,plot3,plot4,plot5,plot6,plot7))
}

sortCPDB <- function(path) {
  pval_path <- paste0(path,"pvalues.txt")
  sigmean_path <- paste0(path,"significant_means.txt")
  message(paste0("Looking for the following files:\n",pval_path,"\n",sigmean_path))
  pval_full <- read.table(pval_path, sep="\t", header=T)
  sigmean_full <- read.table(sigmean_path, sep="\t", header=T)
  
  #Remove blood and leukocytes
  idx_blood <- colnames(pval_full)[grep("Blood",colnames(pval_full))]
  idx_leukocytes <-colnames(pval_full)[grep("Leukocytes",colnames(pval_full))] 
  
  pval_full <- pval_full[,!colnames(pval_full) %in% idx_blood]
  pval_full <- pval_full[,!colnames(pval_full) %in% idx_leukocytes]
  pval <- pval_full[,12:36]
  rownames(pval) <- pval_full$id_cp_interaction
  pval[pval>0.1] <- 1
  
  sigmean_full <- sigmean_full[,!colnames(sigmean_full) %in% idx_blood]
  sigmean_full <- sigmean_full[,!colnames(sigmean_full) %in% idx_leukocytes]
  sigmean <- sigmean_full[,13:37]
  rownames(sigmean) <- sigmean_full$id_cp_interaction
  
  #Check consistency between data
  message(paste0(length(setdiff(rownames(pval),rownames(sigmean)))," pairs mismatching between matrices."))
  
  #Interaction cell types W/O any significant interactions
  cs_pval <- pval %>% colSums
  idx_cs <- cs_pval[cs_pval==dim(pval_full)[1]] %>% names
  
  pval_full <- pval_full[,!colnames(pval_full) %in% idx_cs]
  sigmean <- sigmean[,!colnames(sigmean) %in% idx_cs]
  sigmean_full <- sigmean_full[,!colnames(sigmean_full) %in% idx_cs]
  
  #Interaction pairs W/O any significant interactions
  rs_pval <- rowSums(pval)
  idx_rs <- names(rs_pval[rs_pval==25])
  
  pval_full <- pval_full[!pval_full$id_cp_interaction %in% idx_rs,]
  sigmean <- sigmean[!rownames(sigmean) %in% idx_rs,]
  sigmean_full <- sigmean_full[!sigmean_full$id_cp_interaction %in% idx_rs,]
  
  message(paste0(sum(colSums(sigmean, na.rm=T)==0)," interactions should be removed columnwise."))
  message(paste0(sum(rowSums(sigmean, na.rm=T)==0)," interactions should be removed rowwise."))
  
  #Save cleaned tables
  message(paste0("Saving tables with ",dim(sigmean)[1]," interaction pairs."))
  write.table(pval_full, paste0(path,"pvalues_clean.txt"), sep="\t", col.names=T, row.names=F)
  sigmean_full[is.na(sigmean_full)] <- ""
  write.table(sigmean_full, paste0(path,"significant_means_clean.txt"), sep="\t", col.names=T, row.names=F)
  message("All done!")
}

renameAnnotation <- function(annotation, old, new) {
  if(!is.factor(annotation)) stop("Annotation must be a factor.")
  
  levels(annotation)[levels(annotation) %in% old] <- new
  
  return(annotation)
}

dotSize <- function(size, alpha=1) {
  ggplot2::guides(colour = guide_legend(override.aes = list(size=size,
                                                            alpha=alpha)))
}

checkDims <- function(cm, con) {
  cat("Dimensions of cm : ",paste((dim(cm)), collapse=" "),"\n")
  cat("Dimensions of con: ",paste((dim(con$embedding)), collapse=" "),"\n")
  
  if(dim(cm)[2]!=dim(con$embedding)[1])
    stop("Dimensions don't match.")
  
  message("All OK!")
}