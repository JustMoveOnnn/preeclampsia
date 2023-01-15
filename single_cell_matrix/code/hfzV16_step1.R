hfzV16_step1.R

# step2_202007.R

library(Seurat)
library(cowplot)
source('/home/wzk/tools/BEER-master/BEER.R')
setwd("/home/yjingjing/project/hongfangzi/hongfangzi_MATRIXrawdata")
#.set_python(‘/usr/bin/python’)
.set_python('/usr/bin/python')

print(getwd())
setwd('/home/yjingjing/project/hfz2020adjust/pla_adjust/6th')
# setwd('/home/yjingjing/project/hfz2020adjust/deci_adjust/6th')

#####################
DATA=readRDS('/home/yjingjing/project/hongfangzi/hongfangzi_MATRIXrawdata/DATA.RDS')
BATCH=readRDS('/home/yjingjing/project/hongfangzi/hongfangzi_MATRIXrawdata/BATCH.RDS')
dim(DATA)
# USED=which(BATCH %in% c('decidua0117','decidua0417.2',
#                        'decidua2018c','decidua20190215',
#                        'decidua20190420', 'decidua2019c',
#                        'decidua508','decidua510'
#                        ))

# USED=which(BATCH %in% c('pla1' ,'pla2', 'pla3', 'pla4', 'pla5', 'pla6', 'pla7'))
USED=which(BATCH %in% c('pla1-NM', 'pla2-NM' ,'pla3-NM',  'pla4-PE' ,'pla5-PE' , 'pla6-PE',  'pla7-PE'))
# USED=which(BATCH %in% c('pla1-NM', 'pla2-NM' ,'pla3-NM',  'pla4-PE' ,'pla5-PE' ,  'pla7-PE'))

# USED=which(BATCH %in% c('deci1-NM', 'deci2-NM','deci3-NM','deci4-PE','deci5-NM','deci6-PE','deci7-NM','deci8-PE'))

DATA=DATA[,USED]
BATCH=BATCH[USED]

pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
Idents(pbmc)=BATCH
pbmc@meta.data$batch=BATCH
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


pdf('HFZ_QC.pdf',width=12,height=7)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 25)
pbmc <- subset(pbmc, subset = nFeature_RNA > 700 & nFeature_RNA < 4000 & percent.mt < 25)

dim(pbmc)

#########
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

BATCH=pbmc@meta.data$batch
DATA=as.matrix(pbmc@assays$RNA@counts[,which(colnames(pbmc@assays$RNA@counts) %in% colnames(pbmc@assays$RNA@data))])




# ************************************************************   # 
# **********************      BEER     ***********************   # 
# ************************************************************   # 


source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
BEER <- function(DATA, BATCH,  GNUM=30, PCNUM=150, GN=1000, CPU=4, COMBAT=TRUE, print_step=10, SEED=123, N=2, ROUND=1, RMG=NULL){

    set.seed( SEED)
    RESULT=list()
    library(Seurat)
    #source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
    print('BEER start!')
    print(Sys.time())
    DATA=DATA
    BATCH=BATCH
    RMG=RMG
    COMBAT=COMBAT
    COMBAT.EXP=NULL
    
    
    require(stringi)
    BATCH=stri_replace_all(BATCH, '.',fixed='_')
    CPU=CPU
    GNUM=GNUM
    PCNUM=PCNUM
    UBATCH=unique(BATCH)
    ROUND=ROUND
    
    GN=GN
    N=N
    print_step=print_step
    
    print('Group number (GNUM) is:')
    print(GNUM)
    print('Varible gene number (GN) of each batch is:')
    print(GN)
    print('ROUND is:')
    print(ROUND)
    
    VARG=c()
    i=1
    for(this_batch in UBATCH){
        print(i)
        i=i+1
        print(this_batch)
        this_pbmc=CreateSeuratObject(counts = DATA[,which(BATCH==this_batch)], min.cells = 0, 
                                 min.features = 0, project = this_batch)
        this_pbmc <- NormalizeData(object = this_pbmc, normalization.method = "LogNormalize", 
                           scale.factor = 10000)
        this_pbmc <- FindVariableFeatures(object = this_pbmc, selection.method = "vst", nfeatures = GN)  
        this_varg=VariableFeatures(object = this_pbmc)
        VARG=c(VARG, this_varg)
        }
    VARG=unique(VARG)

    print('Total varible gene number (GN) is:')
    print(length(VARG))

    pbmc=CreateSeuratObject(counts = DATA, min.cells = 0, min.features = 0, project = "ALL") 
    pbmc@meta.data$batch=BATCH
    VariableFeatures(object = pbmc)=VARG

    ########
    if(!is.null(RMG)){
        print('Total removed gene number is:')
        print(length(RMG))
        VariableFeatures(object = pbmc)=VariableFeatures(object = pbmc)[which(! VariableFeatures(object = pbmc) %in% RMG)]
        print('Total used gene number is:')
        print(length(VariableFeatures(object = pbmc)))
        }
    ##########   

    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


    if(COMBAT==FALSE){
        pbmc <- ScaleData(object = pbmc, features = VariableFeatures(object = pbmc),)
    }else{
        ##############
        library(sva)
        library(limma)
        pheno = data.frame(batch=as.matrix(BATCH))
        orig.data=pbmc@assays$RNA@data
        used.gene.index=which(rownames(orig.data) %in% VARG)
        edata = as.matrix(orig.data)[used.gene.index,]
        batch = pheno$batch
        modcombat = model.matrix(~1, data=pheno)
        combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
        rownames(combat_edata)=rownames(edata)
        colnames(combat_edata)=colnames(edata)
        combat_edata=as.matrix(combat_edata)
        combat_edata[which(combat_edata<0)]=0
        combat_edata[which(is.na(combat_edata))]=0
        pbmc@assays$RNA@data=combat_edata
        ######


        data(cc.genes)
        s.genes <- cc.genes$"s.genes"
        g2m.genes <- cc.genes$"g2m.genes"
        # pbmc <- CellCycleScoring(object = pbmc, s.genes = s.genes, g2m.genes = g2m.genes,set.ident = FALSE)
        pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) 
        
        # # scObject <- ScaleData(object = scObject, vars.to.regress = c("nUMI", "percent.mito", "S.Score", "G2M.Score"), do.par = TRUE, num.cores = 8)
        # ######

        # pbmc <- ScaleData(object = pbmc, features = rownames(pbmc), vars.to.regress = c("nCount_RNA","percent.mt",'S.Score','G2M.Score'))
        # # Alternate Workflow
        pbmc$CC.Difference <- pbmc$S.Score - pbmc$G2M.Score
        pbmc <- ScaleData(pbmc, vars.to.regress = "CC.Difference", features = rownames(pbmc))

        # vars.to.regress = c("nUMI", "percent.mito")
        ######
        pbmc@assays$RNA@data=orig.data    
        COMBAT.EXP=combat_edata
        #################
        rm(edata)
        rm(combat_edata)
        rm(orig.data)
        gc()
    }
    print('Calculating PCs ...')
    pbmc <- RunPCA(object = pbmc, seed.use=SEED, npcs=PCNUM, features = VariableFeatures(object = pbmc), ndims.print=1,nfeatures.print=1)
    pbmc <- RunUMAP(pbmc, dims = 1:PCNUM,seed.use = SEED,n.components=N)


    DR=pbmc@reductions$umap@cell.embeddings
    GROUP=rep('NA',length(BATCH))
    for(this_batch in UBATCH){
        this_index=which(BATCH==this_batch)
        this_one=DR[this_index,]
        #CNUM=max(c(5, round(length(this_index)/GNUM) ))
        this_gnum=min(GNUM, (length(this_index)-1))
        this_group=.getGroup(this_one,this_batch,this_gnum)
        GROUP[this_index]=this_group
    }


    pbmc@meta.data$group=GROUP
    
    ##########
    #VP=.getVPnet(pbmc, ROUND)
    VP=.getVPall(pbmc, ROUND)
    ##########
    
    DR=pbmc@reductions$pca@cell.embeddings  

    MAP=rep('NA',length(GROUP))
    MAP[which(GROUP %in% VP[,1])]='V1'
    MAP[which(GROUP %in% VP[,2])]='V2'
    pbmc@meta.data$map=MAP
    
    OUT=.evaluateProBEER(DR, GROUP, VP)

    RESULT=list()
    RESULT$seurat=pbmc
    RESULT$vp=VP
    #RESULT$vpcor=VP_OUT$cor
    RESULT$cor=OUT$cor
    RESULT$pv=OUT$pv
    RESULT$fdr=OUT$fdr
    RESULT$lcor=OUT$lcor
    RESULT$lpv=OUT$lpv
    RESULT$lc1=OUT$lc1
    RESULT$lc2=OUT$lc2
    RESULT$lfdr=OUT$lfdr
    
    ################
    RESULT$ROUND=ROUND
    RESULT$COMBAT=COMBAT
    RESULT$COMBAT.EXP=COMBAT.EXP
    RESULT$RMG=RMG
    RESULT$GNUM=GNUM
    RESULT$GN=GN
    RESULT$PCNUM=PCNUM
    RESULT$SEED=SEED
    RESULT$N=N
    RESULT$APP='BEER'   
    ###############
    
    
    PCUSE=which( (rank(RESULT$cor)>=length(RESULT$cor)/2 | RESULT$cor>0.7 )    & 
                (rank(RESULT$lcor) >=length(RESULT$cor)/2 | RESULT$lcor>0.7)   #&
                #p.adjust(RESULT$lc1,method='fdr') >0.05
               ) 
    
    RESULT$select=PCUSE
    
    print('############################################################################')
    print('BEER cheers !!! All main steps finished.')
    print('############################################################################')
    print(Sys.time())

    return(RESULT)
}

ambient_RNA <- c('PAEP', 'HBG1', 'HBA1', 'HBA2', 'HBM', 'AHSP' , 'HBG2')
RMG <- ambient_RNA
mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=150, ROUND=1, GN=1000, SEED=1, COMBAT=TRUE, RMG=NULL)

# pdf('HFZ1.pdf',width=10,height=10)

PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))
getwd()
# [1] "/home/yjingjing/project/hfz2020adjust/deci_adjust/6th"
saveRDS(mybeer,file='mybeer.RDS')



# Check selected PCs
# PCUSE=mybeer$select
# PCUSE=.selectUSE(mybeer, CUTR=0.9, CUTL=0.9, RR=0.5, RL=0.5)
# PCUSE=.selectUSE(mybeer, CUTR=0.7, CUTL=0.7, RR=0.5, RL=0.5)

# COL=rep('black',length(mybeer$cor))
# COL[PCUSE]='red'
# # pdf('2020_0110hfz.pdf')
# plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
#     xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))
# # dev.off()

# saveRDS(mybeer,'mybeer_decidua.RDS')


