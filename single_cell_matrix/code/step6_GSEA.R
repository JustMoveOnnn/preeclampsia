library(stringr)
pbmc@meta.data[,c("loc","disease")]] <- str_split_fixed(pbmc@meta.data$loc,"-",2)
tail(pbmc@meta.data)

pbmc@meta.data$loc1 <- apply(pbmc@meta.data, 1, function(x){substr(as.character(x['loc']), 1, nchar(as.character(x['loc']))-1)})



deciNM=which(pbmc@meta.data$loc1=='deci' & pbmc@meta.data$disease=='NM' & pbmc@meta.data$loc !='deci9')

deciNM <- subset(pbmc, subset = loc1 =='deci' & disease == 'NM')
dim(deciNM)
# [1] 31733 13089
deciPE <- subset(pbmc, subset = loc1 =='deci' & disease == 'PE')
dim(deciPE)
# [1] 31733 10976
plaNM <- subset(pbmc, subset = loc1 =='pla' & disease == 'NM')
dim(plaNM)
# [1] 31733  6835
plaPE <- subset(pbmc, subset = loc1 =='pla' & disease == 'PE')
dim(plaPE)
# [1] 31733  7652


TAB=table(pbmc@meta.data$level1, pbmc@meta.data$batch)
.writeTable(TAB, PATH='TAB_hfz.txt')

pbmc@meta.data <- unite(pbmc@meta.data,"loc1_disease",loc1,disease)

############################***************** deci-T  ************############################ 
USE.CELL=which(pbmc@meta.data$level1=='T' & pbmc@meta.data$loc1=='deci' & pbmc@meta.data$loc !='deci9')
TAG=pbmc@meta.data$loc1_disease
EXP=as.matrix(pbmc@assays$RNA@data[,USE.CELL])
VAR=apply(EXP,1,var)
EXP=EXP[which(VAR>0),]
PT=t(as.character(TAG[USE.CELL]))
OUT=cbind(toupper(rownames(EXP)),rep('NO',nrow(EXP)),EXP)
colnames(OUT)[c(1,2)]=c('GENE','DESCRIPTION')
# change OUT colnames
reformat <- function(x){
 p1 = strsplit(as.character(x), split='-')[[1]][1]
 p2 = strsplit(as.character(x), split='-')[[1]][2]

 p11 = substr(p1, 1, nchar(p1)-1)

 merge = paste(p11, p2, sep='')
 return(merge)
}

ori = as.character(colnames(OUT))
new = c()

for (i in ori){
 new = c(new, reformat(i))
}

colnames(OUT) = new
colnames(OUT)[c(1,2)]=c('GENE','DESCRIPTION')
# change PT format
n = dim(PT)[2]
p11 = substr(p1, 1, nchar(p1)-1)
line1 = paste(as.character(n), '2', '1', sep=' ')
line2 = paste(c('#', unique(as.character(PT))), collapse=' ')
f <- file("./GSEA_cell_data/deci_T_PT.cls", "w")
writeLines(line1, f)
writeLines(line2, f)
close(f)

write.table(PT, file="./GSEA_cell_data/deci_T_PT.cls", append=T, quote=F, sep=" ", row.names=F, col.names=F)
####
write.table(OUT,'./GSEA_cell_data/deci_T_EXP.txt',sep='\t',quote=F,row.names=F,col.names=T)
##################################

########################***************** deci-DC_M  ************############################ 
USE.CELL=which(pbmc@meta.data$level1=='DC_M' & pbmc@meta.data$loc1=='deci' & pbmc@meta.data$loc !='deci9')
TAG=pbmc@meta.data$loc1_disease
EXP=as.matrix(pbmc@assays$RNA@data[,USE.CELL])
VAR=apply(EXP,1,var)
EXP=EXP[which(VAR>0),]
PT=t(as.character(TAG[USE.CELL]))
OUT=cbind(toupper(rownames(EXP)),rep('NO',nrow(EXP)),EXP)
colnames(OUT)[c(1,2)]=c('GENE','DESCRIPTION')
# change OUT colnames
reformat <- function(x){
 p1 = strsplit(as.character(x), split='-')[[1]][1]
 p2 = strsplit(as.character(x), split='-')[[1]][2]

 p11 = substr(p1, 1, nchar(p1)-1)

 merge = paste(p11, p2, sep='')
 return(merge)
}

ori = as.character(colnames(OUT))
new = c()

for (i in ori){
 new = c(new, reformat(i))
}

colnames(OUT) = new
colnames(OUT)[c(1,2)]=c('GENE','DESCRIPTION')
# change PT format
n = dim(PT)[2]
line1 = paste(as.character(n), '2', '1', sep=' ')
rename <- function(x){
 p1 = strsplit(as.character(x), split='_')[[1]][1]
 p2 = strsplit(as.character(x), split='_')[[1]][2]

 merge = paste(p1, p2, sep='')
 return(merge)
}
PT = as.character(sapply(as.character(PT), rename))
line3 = paste(PT, collapse=' ')
line2 = paste(c('#', unique(as.character(PT))), collapse=' ')
f <- file("./GSEA_cell_data/deci_DC_M_PT.cls", "w")
writeLines(line1, f)
writeLines(line2, f)
writeLines(line3, f)
close(f)

# write.table(PT, file="./GSEA_cell_data/deci_DC_M_PT.cls", append=T, quote=F, sep=" ", row.names=F, col.names=F)
####
write.table(OUT,'./GSEA_cell_data/deci_DC_M_EXP.txt',sep='\t',quote=F,row.names=F,col.names=T)
# write.table(PT,'./GSEA_cell_data/deci_DC_M_PT.cls',sep=' ',quote=F,row.names=F,col.names=F )
##################################

############################***************** deci-T  ************############################ 

############################***************** deci-Endo  ************############################ 


############################***************** deci-NK  ************############################ 

############################***************** deci-Endo  ************############################ 





