---
title: "On benchmarking gene set enrichment analysis"
output: html_document
---



## Introduction

Systematic evaluation of gene set enrichment analysis methods 
is a much needed but challenging task. 
Recently Geistlinger and colleagues presented a framework
[GSEABenchmarkeR](https://bioconductor.org/packages/GSEABenchmarkeR)
for extensible and reproducible gene set enrichment benchmarking. 
Here we show that the MalaCards-based phenotype relevance score,
one of the key framework's component,
turns out to be a poor metric for comparing quality of enrichment results.
We demonstrate this by showing that a baseline method that simply orders genes sets by their 
size significantly outperforms all other considered methods, including frequently used
in practice ORA and GSEA methods.


### Loading libraries


```r
library(GSEABenchmarkeR)
library(EnrichmentBrowser)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggpubr)

source("methods/utils.R")
```

### Microarray and TCGA ids


```r
geo2kegg <- loadEData("geo2kegg")
ma.ids <- names(geo2kegg)

rseq.ids <- c("BLCA", "BRCA", "COAD", "HNSC", "KICH", 
              "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", 
              "PRAD", "READ", "STAD", "THCA", "UCEC")
```

### Reading rankings from the article


```r
ea.methods <- sbeaMethods()[1:10]

rseq.kegg.ranks <- readResults("inpdata/TCGA/GSE62944_matched_vst/kegg/", rseq.ids,
                               methods = ea.methods, type="ranking")
rseq.go.ranks <- readResults("inpdata/TCGA/GSE62944_matched_vst/gobp/", rseq.ids,
                               methods = ea.methods, type="ranking")

ma.kegg.ranks <- readResults("inpdata/GEO2KEGG/kegg/perm1k/", ma.ids, 
                             methods = ea.methods, type="ranking")
ma.go.ranks <- readResults("inpdata/GEO2KEGG/go_bp/", ma.ids, 
                             methods = ea.methods, type="ranking")
```

## Benchmarking

Here we introduce a baseline enrichment method `GSS`, which ranks gene sets based only on their size.

### Adding geneset-size based ranking


```r
gss.rseq.kegg <- lapply(rseq.ids, function(rseqId){
    oraRanking <- as.data.table(rseq.kegg.ranks$ora[[rseqId]])
    gssRanking <- oraRanking[, .(GENE.SET, NR.GENES)]
    gssRanking[, PVAL := 1 / NR.GENES]
    return(as.data.frame(gssRanking))
})
names(gss.rseq.kegg) <- rseq.ids

gss.rseq.go <- lapply(rseq.ids, function(rseqId){
    oraRanking <- as.data.table(rseq.go.ranks$ora[[rseqId]])
    gssRanking <- oraRanking[, .(GENE.SET, NR.GENES)]
    gssRanking[, PVAL := 1 / NR.GENES]
    return(as.data.frame(gssRanking))
})
names(gss.rseq.go) <- rseq.ids

gss.ma.kegg <- lapply(ma.ids, function(maId){
    oraRanking <- as.data.table(ma.kegg.ranks$ora[[maId]])
    gssRanking <- oraRanking[, .(GENE.SET, NR.GENES)]
    gssRanking[, PVAL := 1 / NR.GENES]
    return(as.data.frame(gssRanking))
})
names(gss.ma.kegg) <- ma.ids

gss.ma.go <- lapply(ma.ids, function(maId){
    oraRanking <- as.data.table(ma.go.ranks$ora[[maId]])
    gssRanking <- oraRanking[, .(GENE.SET, NR.GENES)]
    gssRanking[, PVAL := 1 / NR.GENES]
    return(as.data.frame(gssRanking))
})
names(gss.ma.go) <- ma.ids

ma.kegg.ranks$gss <- gss.ma.kegg
ma.go.ranks$gss <- gss.ma.go
rseq.kegg.ranks$gss <- gss.rseq.kegg
rseq.go.ranks$gss <- gss.rseq.go
```

### Reading mala cards


```r
data.dir <- system.file("extdata", package="GSEABenchmarkeR")
mala.kegg.file <- file.path(data.dir, "malacards", "KEGG.rds")
mala.go.file <- file.path(data.dir, "malacards", "GO_BP.rds")
mala.kegg <- readRDS(mala.kegg.file)
mala.go <- readRDS(mala.go.file)

d2d.file <- file.path(data.dir, "malacards", "GseId2Disease.txt")
d2d.map <- readDataId2diseaseCodeMap(d2d.file)

d2d.tcga <- rseq.ids
names(d2d.tcga) <- rseq.ids
```

### Evaluating relative scores


```r
ma.kegg.rel.sets <- evalRelevance(ma.kegg.ranks, mala.kegg, d2d.map)
ma.go.rel.sets <- evalRelevance(ma.go.ranks, mala.go, d2d.map)

rseq.kegg.rel.sets <- evalRelevance(rseq.kegg.ranks, mala.kegg, d2d.tcga)
rseq.go.rel.sets <- evalRelevance(rseq.go.ranks, mala.go, d2d.tcga)
```

### Plot results


```r
facetplot(ma.kegg.rel.sets, ma.go.rel.sets, 
          rseq.kegg.rel.sets, rseq.go.rel.sets, 
          ylab="% optimal relevance score", vline = NA)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

## Overall difference between `GSS` and other methods


```r
vapply(setdiff(ea.methods, "gss"), 
        function(m) testMethods("gss", m), 
        numeric(5))
```

```
##                    ora         safe         gsea          gsa        padog
## ma.kegg   7.314432e-07 3.285342e-08 1.442096e-10 4.096359e-10 8.880535e-06
## ma.go     3.138876e-13 6.732820e-15 1.269404e-14 5.842781e-15 4.091227e-13
## rseq.kegg 8.340902e-04 5.777813e-04 8.819120e-06 4.809257e-06 2.701874e-03
## rseq.go   1.408119e-03 1.408119e-03 8.340902e-04 8.340902e-04 3.155375e-03
## overall   2.124573e-23 2.312544e-26 3.607005e-30 2.168215e-30 4.144824e-21
##             globaltest        roast       camera         gsva        samgs
## ma.kegg   8.981595e-13 7.891506e-13 6.944718e-12 1.406870e-13 6.096286e-12
## ma.go     1.073955e-13 5.442263e-15 7.755993e-15 5.842781e-15 1.459795e-14
## rseq.kegg 3.357998e-06 1.289345e-08 1.289345e-08 1.289345e-08 3.105197e-06
## rseq.go   1.289345e-08 1.177172e-05 2.643028e-04 2.643028e-04 1.289345e-08
## overall   7.923105e-32 1.675462e-34 1.144121e-32 1.126819e-34 3.995841e-33
```


## Session info


```r
sessionInfo()
```

```
## R version 4.1.1 (2021-08-10)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Debian GNU/Linux 10 (buster)
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
## LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.3.5.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=C                  LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] ggpubr_0.4.0                ggplot2_3.3.5              
##  [3] reshape2_1.4.4              data.table_1.14.0          
##  [5] EnrichmentBrowser_2.22.2    graph_1.70.0               
##  [7] GSEABenchmarkeR_1.12.1      SummarizedExperiment_1.22.0
##  [9] GenomicRanges_1.44.0        GenomeInfoDb_1.28.1        
## [11] IRanges_2.26.0              S4Vectors_0.30.0           
## [13] MatrixGenerics_1.4.0        matrixStats_0.59.0         
## [15] Biobase_2.52.0              BiocGenerics_0.38.0        
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-7                        bit64_4.0.5                        
##  [3] filelock_1.0.2                      httr_1.4.2                         
##  [5] ggsci_2.9                           Rgraphviz_2.36.0                   
##  [7] tools_4.1.1                         backports_1.2.1                    
##  [9] utf8_1.2.2                          R6_2.5.0                           
## [11] DBI_1.1.1                           colorspace_2.0-2                   
## [13] withr_2.4.2                         tidyselect_1.1.1                   
## [15] bit_4.0.4                           curl_4.3.2                         
## [17] compiler_4.1.1                      DelayedArray_0.18.0                
## [19] labeling_0.4.2                      KEGGgraph_1.52.0                   
## [21] scales_1.1.1                        rappdirs_0.3.3                     
## [23] digest_0.6.27                       stringr_1.4.0                      
## [25] foreign_0.8-81                      rio_0.5.27                         
## [27] XVector_0.32.0                      pkgconfig_2.0.3                    
## [29] highr_0.9                           dbplyr_2.1.1                       
## [31] fastmap_1.1.0                       rlang_0.4.11                       
## [33] readxl_1.3.1                        RSQLite_2.2.7                      
## [35] farver_2.1.0                        generics_0.1.0                     
## [37] dplyr_1.0.7                         zip_2.2.0                          
## [39] car_3.0-11                          RCurl_1.98-1.3                     
## [41] magrittr_2.0.1                      GenomeInfoDbData_1.2.6             
## [43] Matrix_1.3-4                        Rcpp_1.0.7                         
## [45] munsell_0.5.0                       fansi_0.5.0                        
## [47] abind_1.4-5                         lifecycle_1.0.0                    
## [49] stringi_1.6.2                       carData_3.0-4                      
## [51] zlibbioc_1.38.0                     BiocFileCache_2.0.0                
## [53] KEGGdzPathwaysGEO_1.30.0            plyr_1.8.6                         
## [55] grid_4.1.1                          blob_1.2.1                         
## [57] forcats_0.5.1                       crayon_1.4.1                       
## [59] lattice_0.20-44                     Biostrings_2.60.1                  
## [61] haven_2.4.1                         annotate_1.70.0                    
## [63] hms_1.1.0                           KEGGREST_1.32.0                    
## [65] knitr_1.33                          pillar_1.6.2                       
## [67] ggsignif_0.6.2                      KEGGandMetacoreDzPathwaysGEO_1.12.0
## [69] XML_3.99-0.6                        glue_1.4.2                         
## [71] evaluate_0.14                       png_0.1-7                          
## [73] vctrs_0.3.8                         cellranger_1.1.0                   
## [75] gtable_0.3.0                        purrr_0.3.4                        
## [77] tidyr_1.1.3                         assertthat_0.2.1                   
## [79] cachem_1.0.5                        xfun_0.24                          
## [81] openxlsx_4.2.4                      xtable_1.8-4                       
## [83] broom_0.7.9                         rstatix_0.7.0                      
## [85] tibble_3.1.3                        AnnotationDbi_1.54.1               
## [87] memoise_2.0.0                       ellipsis_0.3.2                     
## [89] GSEABase_1.54.0
```
