library(Seurat)
library(data.table)
library(EWCE)

mat <- fread("https://cells.ucsc.edu/cortex-dev/exprMatrix.tsv.gz")
meta <- data.frame(fread("https://cells.ucsc.edu/cortex-dev/meta.tsv"), row.names=1)
HC_genes = readLines('data_for_Figure6d.txt')

genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)
so <- CreateSeuratObject(counts = mat, project = "Cortex", meta.data=meta)
so <- NormalizeData(so)

so <- AddModuleScore(so, HC_genes)
FeaturePlot(so, features="HC_genes", order=T, max.cutoff="q99", reduction = "tsne") + scale_color_viridis()

# downstream editing was done in the Adobe Illustrator