# monocle.pla_NM.R
# monocle.pla_PE.R
# monocle.deci_NM.R
# monocle.deci_PE.R
# step_monocle.R
# monocle包有很多种读取数据的方式，这里只展示了读取Seurat中的对象的方法
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("monocle")

rm(list=ls())

library(Seurat)
library(monocle)
library(plyr)
library(ggplot2)
pbmc <- readRDS("pbmc20191119.RDS")
# ****** select cluster you want **** #
######### ***************** pla_NM ***************** ###############
exp_data1=subset(x=pbmc,subset = annoClust=="CTB" & loc1_disease=="deci_NM" )
exp_data2=subset(x=pbmc,subset = annoClust=="EVT" & loc1_disease=="deci_NM" )
exp_data3=subset(x=pbmc,subset = annoClust=="STB" & loc1_disease=="deci_NM" )
exp_data12 <- merge(x=exp_data1, y=exp_data2)
exp_data <- m
.erge(x=exp_data12, y=exp_data3)
count_matrix <- exp_data@assays$RNA@counts
cluster <- exp_data@meta.data$annoClust
table(cluster)
names(count_matrix)
######### ***************** pla_PE ***************** ###############


exp_data=as.matrix(count_matrix)


# sample barcode VS cluster
sample_ann <- data.frame(cells=colnames(count_matrix), cellType=cluster)
row.names(sample_ann) = sample_ann$cells
sample_ann <- sample_ann[,-1]

names(sample_ann)=colnames(exp_data)
sample_ann=as.matrix(sample_ann)
sample_ann=as.data.frame(sample_ann)

colnames(sample_ann)='cellType'
head(sample_ann)




# exp_data=as.matrix(pbmc@assays$RNA@scale.data)
#.  use: scale.data or counts ??? .#

all_genes=rownames(exp_data) 

col_data=colnames(exp_data)
names(col_data)=colnames(exp_data) # what ???

col_data=as.data.frame(col_data)

row_data=rownames(exp_data)
names(row_data)=rownames(exp_data)

row_data=as.matrix(row_data)
colnames(row_data)='gene_short_name'
row_data=as.data.frame(row_data)
head(col_data)
head(row_data)
pd <- new("AnnotatedDataFrame", data =sample_ann)
fd <- new("AnnotatedDataFrame", data =row_data)
HSMM=newCellDataSet(exp_data, phenoData = pd, featureData = fd,expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
HSMM_myo=HSMM
remove(HSMM)
# 进行基因过滤
HSMM_myo <- detectGenes(HSMM_myo, min_expr = 0.1)
head(HSMM_myo@featureData@data)
expressed_genes <- row.names(subset(HSMM_myo@featureData@data,
                                    num_cells_expressed >= 5))

length(expressed_genes)

HSMM_myo <- HSMM_myo[expressed_genes,]

# ordering_genes = all_genes
#ordering_genes = c('Olig2','Olig1')

HSMM_myo <- estimateSizeFactors(HSMM_myo)
HSMM_myo = estimateDispersions(HSMM_myo)
# HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
##### 进一步挑基因进行 ###
disp_table <- dispersionTable(HSMM_myo) # 挑有差异的
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1) # 挑表达量不太低的
HSMM_myo <- setOrderingFilter(HSMM_myo, unsup_clustering_genes$gene_id)  # 准备聚类基因名单

head(HSMM_myo@phenoData@data)
pdf('./deci_NM/plot_ordering_genes.pdf')
plot_ordering_genes(HSMM_myo) 
##### 进一步挑基因进行 ###
dev.off()

HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2,  num_dim = 6, reduction_method = 'DDRTree')
HSMM_myo <- orderCells(HSMM_myo) # break the CPU
pdf('./deci_NM/plot_cell_trajectory.State.pdf')
plot_cell_trajectory(HSMM_myo, color_by = "State")
dev.off()
pdf('./deci_NM/plot_cell_trajectory.Pseudotime.pdf')
plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime")
dev.off()

pdf('./deci_NM/plot_cell_trajectory.cellType.pdf')
plot_cell_trajectory(HSMM_myo, color_by = "cellType")
dev.off()

pdf('./deci_NM/fenkaihua.pdf')
plot_cell_trajectory(HSMM_myo, color_by = "cellType" ) + facet_wrap(~cellType,nrow=1)
dev.off()


show_genes <- row.names(subset(fData(HSMM_myo), gene_short_name %in% c("T","APC","NK1")))
pdf('./deci_NM/plot_genes_jitter.pdf')
plot_genes_jitter(HSMM_myo[show_genes,], grouping = "State", min_expr = 0.1, color_by = "State")
dev.off()









# ************************************************************************************  #
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("monocle")
library(monocle)


MARKER=read.table('Marker_GENE.txt',header=T)
THIS=which(MARKER[,7] %in% c('Astro','RG','OPC','IGC') & MARKER[,6]<0.05)
#THIS=which( MARKER[,6]<0.05)
MARKER_GENE=MARKER[THIS,1]
#MARKER_GENE=c('Olig2','Myc')



HSMM_expr_matrix <- read.table('2000_run1716_normalized-FOR QC-after qc_groupasseurat_only135689.txt.tsv',header=T,row.names=1)
HSMM_sample_sheet <- read.delim('2000_run1716_normalized-FOR QC-after qc_groupasseurat_only135689.txt.tsv.pheno',header=F,row.names=1)
colnames(HSMM_sample_sheet)=c('LABEL')
pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)


HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd)
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

print(head(pData(HSMM)))

HSMM_myo <- setOrderingFilter(HSMM,  MARKER_GENE)
plot_ordering_genes(HSMM_myo)

#HSMM_myo <- reduceDimension(HSMM_myo, max_components = 3,method = 'DDRTree')
HSMM_myo <- reduceDimension(HSMM_myo, max_components = 3,method = 'tSNE')
HSMM_myo <- orderCells(HSMM_myo)
plot_cell_trajectory(HSMM_myo)

HSMM_myo <- orderCells(HSMM_myo,root_state=4)
plot_cell_trajectory(HSMM_myo,color_by = "LABEL")
