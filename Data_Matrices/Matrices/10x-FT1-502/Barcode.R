library(Seurat)
library(dplyr)
library(Matrix)
library(DropletUtils)

A <- Read10X(data.dir = "/Critical_Scripts/10x-515/filtered_feature_bc_matrix") 
B <- CreateSeuratObject(counts = A) 
C <- subset(B, subset = PTPRC > 0)
dim(A)
dim(B)
dim(C)
D <-colnames(C)
# 1725 cells

# Now, I will subset to only include cells with common barcodes
E <- as.matrix(A[, D])
E <- A[, D]
dim(E)
save(E,file="/Critical_Scripts/10x-515/E.Robj")
write10xCounts(path ="/Critical_Scripts/10x-515/filtered_feature_bc_matrix", E, barcodes = colnames(E), gene.id = rownames(E), overwrite = TRUE, type = c("sparse"), version = c("3"))
