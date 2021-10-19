# On benchmarking gene set enrichment analysis


Systematic evaluation of gene set enrichment analysis methods 
is a much needed but challenging task. 
Recently Geistlinger and colleagues presented a framework [GSEABenchmarkeR](https://bioconductor.org/packages/GSEABenchmarkeR) for extensible and reproducible gene set enrichment benchmarking. 
Here we show that the MalaCards-based phenotype relevance score,
one of the key framework's component,
turns out to be a poor metric for comparing quality of enrichment results.
We demonstrate this by showing that a baseline method that simply orders genes sets by their 
size significantly outperforms all other considered methods, including frequently used
in practice ORA and GSEA methods.

The built vignette is available at https://ctlab.github.io/gseabenchmarker-counterexample/main.html

The code and data is based on https://github.com/waldronlab/GSEABenchmarking repository.
