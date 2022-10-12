# Setup
install.packages("pheatmap")
library(pheatmap)
library(ggplot2)
library(readxl)
library(RColorBrewer)


# Load
asd.genes.all <- read_excel("data/NDD_control_gene_sets_20210512.xlsx")

anno <- read.csv("data/human_ctx_metadata.csv")
anno$cellid <- make.names(anno$cellid)
row.names(anno) <- anno$cellid

expr <- read.csv("data/human_ctx_trimmed_means.csv")
# keep.genes <- which(expr$feature %in% asd.genes.all$gene)
keep.genes <- which(! duplicated(expr$feature))
expr <- expr[keep.genes, ]
row.names(expr) <- expr$feature
expr <- expr[, anno$cellid]
expr <- as.matrix(expr)


#### Cell type gene counts ####
cl.class <- anno$level2class
gene.sets <- c("Control_DNM447", "NDD_SYN884", "Coe253", "ASC102", "DDD285", 
                            "least615", "most138")

exprl <- list()
for (set1 in gene.sets) {
  test.genes <-
    intersect(asd.genes.all$gene[asd.genes.all$set == set1],
              row.names(expr))
  expr.mean.subset <- expr[test.genes, ]
  
  exprl[[set1]] <- 
  t(apply(expr.mean.subset, 1, function(x)
    tapply(x, cl.class, mean)))
}

pw.ks <- function(sets, class1, exprl) {
  nset <- length(sets)
  ks.m <- matrix(NA, nrow = nset, ncol = nset, 
                 dimnames = list(sets, sets))
  for (i in 1:nset) {
    expr1 <- exprl[[sets[i]]][, class1]
    for (j in 1:nset) {
      expr2 <- exprl[[sets[j]]][, class1]
      ks1 <- ks.test(expr1, expr2, alternative = "greater")
      # ks.stat <- ks1$statistic
      ks.stat <- -log10(ks1$p.value)
      ks.m[i, j] <- ks.stat
    }
  }
  return(ks.m)
}

cl.classes <- unique(cl.class)
cl.classes <- setdiff(cl.classes, c("VLMC", "Endothelial", "Pericyte"))  # Remove rare types

ks.hm <- list()
for (class1 in cl.classes) {
  ks.hm[[class1]] <- pw.ks(gene.sets, class1, exprl)
  # pheatmap(ks.hm[[class1]], cluster_rows = FALSE, cluster_cols = FALSE, main = class1)
}


# Viz
pal1 <- colorRampPalette(brewer.pal(8, "Blues"))(100)
# barplot(sort(sapply(ks.hm, function(x) x[1, 4])), horiz = TRUE, las = 2)

ks.v.control <- sapply(ks.hm, function(x) x[1, ])
ks.v.control <- ks.v.control[-1, ]
ks.v.control <- 10^(-ks.v.control) * length(ks.v.control)
ks.v.control <- -log10(ks.v.control)
ks.v.control[ks.v.control < 0] <- 0

pdf(file = "output/kstest_sig_hm.pdf", width = 6, height = 4)
pheatmap(ks.v.control[, order(colnames(ks.v.control))], cluster_rows = FALSE, 
         cluster_cols = TRUE, color = pal1)
dev.off()

class1 <- "LAMP5"
pdf(file = paste0("output/", class1, "_ecdf2.pdf"), width = 6, height = 5)
plot(ecdf(exprl[["Control_DNM447"]][, class1]), xlab = "Log2 expression threshold", 
          ylab = "Prop. of genes with lower expression", main = class1)
lines(ecdf(exprl[["least615"]][, class1]), col = "light blue")
lines(ecdf(exprl[["most138"]][, class1]), col = "blue")
dev.off()




#### Count gene cluster expression ####
analysis.path <- "output/"


cl.class1 <- sub("-", "_", anno$level1class)

cl.cnt <- apply(expr, 1, function(x) {
  tapply(x, cl.class1, function(y) sum(y > 1))
})

marker.scores <- data.frame(gene = row.names(expr), t(cl.cnt))

asd.marker <- merge(marker.scores, asd.genes.all[, c("gene", "set_type", "set")], 
                    by = "gene", all.y = TRUE)
asd.marker <- droplevels(na.omit(asd.marker))

asd.marker.subset <- droplevels(subset(asd.marker, set %in% gene.sets))
write.csv(asd.marker.subset, file = paste0(analysis.path, "asd.gene.expr.csv"), 
          row.names = FALSE)


cell.classes <- unique(cl.class1)

for (cell.class in cell.classes) {
  g.dens.subset <- ggplot(asd.marker.subset, aes_string(x = cell.class,
                                                        color = "set", linetype = "set_type")) +
    stat_ecdf(geom = "step", size = 1) +
    xlab("Number of clusters with expression") +
    ylab("Proportion of genes") +
    ggtitle(cell.class) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  plot(g.dens.subset)
  ggsave(g.dens.subset, filename = paste0(analysis.path, "g.dens.subset_ctx.", cell.class, ".pdf"),
         width = 6, height = 4)
  
  
  # Stats
  print(cell.class)
  
  asd.marker.compare <- subset(asd.marker, set %in% gene.sets)
  formula1 <- as.formula(paste(cell.class, "~ set"))
  kt1 <- kruskal.test(formula1, data = asd.marker.compare)
  pw.diff <- pairwise.wilcox.test(asd.marker.compare[, cell.class],
                                  asd.marker.compare$set)
  print(pw.diff)
  
  set1 <- asd.marker[asd.marker$set == "most138", cell.class]
  set2 <- asd.marker[asd.marker$set == "ASC102", cell.class]
  wilcox1 <- wilcox.test(set1, set2)
  p.cor <- min(wilcox1$p.value * length(cell.classes), 1)  # Bonferroni correction
  print(p.cor)
  
}




#### EWCE ####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EWCE")
library(EWCE)
library(cowplot)
install.packages("limma")
library(limma)
library(readxl)
library(RColorBrewer)
theme_set(theme_cowplot())



# Load
asd.genes.all <- read_excel("data/NDD_control_gene_sets_20210512.xlsx")

anno <- read.csv("data/human_ctx_metadata.csv")
anno$cellid <- make.names(anno$cellid)
row.names(anno) <- anno$cellid

expr <- read.csv("data/human_ctx_trimmed_means.csv")
# keep.genes <- which(expr$feature %in% asd.genes.all$gene)
keep.genes <- which(! duplicated(expr$feature))
expr <- expr[keep.genes, ]
row.names(expr) <- expr$feature
expr <- expr[, anno$cellid]
expr <- as.matrix(expr)



cortex_mrna <- list(exp = expr, annot = anno)

gene="SLC17A7"
cellExpDist = data.frame(e=cortex_mrna$exp[gene,],l1=cortex_mrna$annot[colnames(cortex_mrna$exp),]$level2class)
ggplot(cellExpDist) + geom_boxplot(aes(x=l1,y=e)) + xlab("Cell type") + ylab("Unique Molecule Count") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


exp2 = drop.uninformative.genes(exp = cortex_mrna$exp,
                                level2annot = cortex_mrna$annot$level2class)
annotLevels = list(
  level1class = cortex_mrna$annot$level1class,
  level2class = cortex_mrna$annot$level2class
)
exp3 = generate.celltype.data(exp = exp2, annotLevels = annotLevels, 
                              groupName = "human_ctx", savePath = "data/")

load("data/CellTypeData_human_ctx.rda", verbose = TRUE)


# Test
reps = 1000 # <- Use 100 bootstrap lists so it runs quickly, for publishable analysis use >10000
level = 2 # <- Use level 1 annotations (i.e. Interneurons)

all.genes <- row.names(ctd[[1]]$mean_exp)
# bg <- all.genes
bg <- unique(asd.genes.all$gene[asd.genes.all$set %in% c("most138", "Control_DNM447")])
fg <- asd.genes.all$gene[asd.genes.all$set == "most138"]  # ASC102  most138  NDD_SYN884  Control_DNM447
# fg <- sample(bg, 138)

full_results = bootstrap.enrichment.test(
  sct_data = ctd,
  hits = fg,
  bg = bg,
  genelistSpecies = "human",
  sctSpecies = "human",
  geneSizeControl = TRUE,
  reps = reps,
  annotLevel = level
)


# print(full_results$results[order(full_results$results$p), 3:5][1:6,])

print(ewce.plot(full_results$results, mtc_method = "BH")$plain)

