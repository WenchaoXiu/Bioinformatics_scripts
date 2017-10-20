'加载所需要的包'
install.packages('ape')
install.packages('gplots')
library(ape)
library(gplots)
library("RDAVIDWebService")
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
require(edgeR)
require(graphics)


'数据的录入'
setwd('~/Desktop/fpkm')
D0_ctrl_o <- read.table('D0_ctrl_o.fpkm',header=T)
D0_esrrb_o <- read.table('D0_esrrb_o.fpkm',header=T)
D2_ctrl_o <- read.table('D2_ctrl_o.fpkm',header=T)
D2_esrrb_o <- read.table('D2_esrrb_o.fpkm',header=T)
D4_ctrl_o <- read.table('D4_ctrl_o.fpkm',header=T)
D4_esrrb_o <- read.table('D4_esrrb_o.fpkm',header=T)
TSC_1 <- read.table('TSC_1.fpkm',header=T)
TSC_2 <- read.table('TSC_2.fpkm',header=T)
rawdata <- cbind(D0_ctrl_o[,10],D0_esrrb_o[,10],D2_ctrl_o[,10],D2_esrrb_o[,10],D4_ctrl_o[,10],D4_esrrb_o[,10],TSC_1[,10],TSC_2[,10])
index <- !(duplicated(D0_ctrl_o$tracking_id))
delete_index <-rownames(Tet2_fpkmdata)=='__no_feature' | rownames(Tet2_fpkmdata)=='__ambiguous' | rownames(Tet2_fpkmdata)=='__alignment_not_unique' | rownames(Tet2_fpkmdata)=='__too_low_aQual' | rownames(Tet2_fpkmdata)=='__not_aligned'
Tet2_fpkmdata <- Tet2_fpkmdata[!delete_index,]
data <- rawdata[index,]
colnames(data) <- c('D0_ctrl_o','D0_esrrb_o','D2_ctrl_o','D2_esrrb_o','D4_ctrl_o','D4_esrrb_o','TSC_1','TSC_2')
rownames(data) <- D0_ctrl_o[index,1]
tdata <- t(data)
'利用gtf求解基因长度'
library(GenomicRanges)
library(rtracklayer)
GTFfile = "/Users/xiuwenchao/Desktop/mm10.gtf"
GTF <- import.gff(GTFfile, format="gtf", genome="mm10", asRangedData=F, feature.type="exon")
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementLengths(grl))
elementMetadata(reducedGTF)$widths <- width(reducedGTF)
calc_length <- function(x) {
sum(elementMetadata(x)$widths)
}
output <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_length))
output <- t(output)
colnames(output) <- c("Length")



'filter数据,对数据的outlier进行处理,对top10%和bottom10%的数据进行处理'
fpkmdata <- na.omit(t(scale(t(fpkmdata))))
data<-fpkmdata
data<-c(as.matrix(data))
data<-sort(data)
min<-data[1]
max<-data[length(data)]
temp<-data[round(c(0.1,0.5,0.9)*length(data))]
p10<-temp[1]
p50<-temp[2]
p90<-temp[3]
zmin=p10
zmax=p90
fpkmdata[fpkmdata<zmin] <- zmin
fpkmdata[fpkmdata>zmax] <- zmax
filter_data <- fpkmdata




#####################################################
############   1.将有重复ID的一列进行取平均   ###########
#####################################################
'将有重复ID的一列进行取平均'
library(data.table)
dat <- read.table(text='name    value   etc1    etc2
A       9       1       3
A       10      3       2
A       11      5       1
B       2       1       4
C       40      2       5
C       50      2       6',header=TRUE)
dat
X <- as.data.table(dat)
X[,lapply(.SD,max),'name']



'对第一列重复的基因保留最大值，rowname处理'
library(data.table)
dat <- read.table('/Volumes/X/aaaaaaa/lung/merge.txt',header=T)
X <- as.data.table(dat)
keys = colnames(dat)[1]
data <- as.data.frame(X[,lapply(.SD,max),keys])
rownames(data) <- data[,1]
data <- data[,-1]

############################################
############   2.array 数据处理   ###########
###########################################
'array 数据处理'
对于array数据下载好原始数据之后,到GSE看芯片注释标签，然后去brainarray网站找对应的annotation
http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/
然后选择Download Custom CDF Files,拉到最新版本,下载refseq版本,再根据注释label寻找R包，下载c一列的(tar.gz文件)
然后shell下 R CMD install *.tar.gz即可
如果在读数据的时候报错，记得下载brainarry下的affy，然后R CMD install *tar.gz在进行文件读取

source("https://bioconductor.org/biocLite.R")
biocLite("affy")

library(affy)
library(moex10stmmrefseqcdf)
setwd('/Users/xiu/Desktop/GSE34243_RAW')
data <- ReadAffy("GSM845636.CEL","GSM845637.CEL","GSM845638.CEL","GSM845639.CEL","GSM845640.CEL","GSM845641.CEL","GSM845642.CEL","GSM845643.CEL"
,"GSM845644.CEL","GSM845645.CEL","GSM845646.CEL","GSM845647.CEL","GSM845648.CEL","GSM845649.CEL","GSM845650.CEL","GSM845651.CEL","GSM845652.CEL")
data@cdfName <- "moex10stmmrefseqcdf"
exp<-exprs(rma(data))


#####################################################################
############   3.利用基因长度以及count值用edger进行RPKM求解   ###########
####################################################################
'利用基因长度以及count值用edger进行RPKM求解'
rawlength <- output
rawdata <- read.csv('/Users/xiuwenchao/Desktop/GSM1599494_ES_d0_main.csv',header=F,sep=',')

data_d0 <- rawdata
rownames(data_d0) <- data_d0[,1]
data_d0 <- data_d0[,-1]
overlap <- intersect(rownames(data_d0),rownames(rawlength))

data_d0 <- data_d0[overlap,]
gene_length <- data.frame(length=rawlength[overlap,])

require(edgeR)
#	利用count计算rpkm,rpkm/2=fpkm(pair end,single end FPKM=RPKM)
RPKM_calculate <- function(count,length=gene_length){
y <- DGEList(counts=count,genes=length)
y <- calcNormFactors(y)
RPKM <- rpkm(y)
return(RPKM)
}
tmp <- apply(data_d0,2,RPKM_calculate)

rownames(tmp) <- overlap
write.table(tmp,'GSM1599494_ES_d0_main.txt',quote=F,sep='\t',row.names=T,col.names=F)




##########################################
############   4.TCGA ID转换   ###########
#########################################
'转换id'
library(mygene)#change gene id

LIHC_exp <- read.csv('/Volumes/TCGA/春节私活/LIHC肝脏/gdac.broadinstitute.org_LIHC.mRNAseq_Preprocess.Level_3.2016012800.0.0/LIHC.uncv2.mRNAseq_RSEM_all.txt',header=T,sep='\t')
'转换id'
entrez_id <- c()
for (i in 1:nrow(LIHC_exp)){
entrez_id <- c(entrez_id,strsplit(as.character(LIHC_exp[i,1]),'|',fixed = TRUE)[[1]][2])
}
convert <- queryMany(entrez_id, scopes="entrezgene", species="human")
entrez_symbol <- data.frame(entrez = convert[,'query'], symbol = convert[,'symbol'])
index <- !is.na(entrez_symbol[,'symbol'])
LIHC_exp <- LIHC_exp[index,-1]
rownames(LIHC_exp) <- entrez_symbol[index,'symbol']
T_LIHC_exp <- t(LIHC_exp)



'根据entrez ID求解对应的gene symbol'
data <- read.csv('~/Desktop/Mus_musculus_transcription_factors_gene_list.txt',header=T,sep='\t')
library(mygene)#change gene id
convert <- queryMany(data[,2], scopes="entrezgene", species="mouse")
entrez_symbol <- data.frame(entrez = convert[,'query'], symbol = convert[,'symbol'])
TF <- as.character(entrez_symbol[!is.na(entrez_symbol[,'symbol']),'symbol'])
TF <- TF[!duplicated(TF)]
write(TF,file='TF_list.txt')



####################################################
############   5.利用edgeR求解差异表达基因   ###########
#####################################################
'利用edgeR求解差异表达基因,并分别记录下来上下调的refseq号和genesymbol'
'详见edger的使用说明'
edger.fun <- function(data,state,rep_num,label,logfc=1.5,pvalue=0.05){
	require(edgeR)
	y <- DGEList(counts = data,group = rep(state,rep_num),genes = rownames(data))
	keep <- rowSums(cpm(y)>15) >= 2 #sum(keep) filter low expressed genes
	y <- y[keep, , keep.lib.sizes=FALSE]
	y <-calcNormFactors(y) #delete effective library size.
	pdf(paste(label,".pdf",sep=""),width=12,height=4)
	par(mfrow=c(1,3))
	plotMDS(y)#看样本之间的关系
	y <- estimateCommonDisp(y)
	y <- estimateTagwiseDisp(y)
	plotBCV(y)
	et<-exactTest(y, dispersion = "auto")
	write.table(topTags(et,n=40000),file=paste('all',label,".txt",sep="_"),sep="\t",quote=F,row.names=F,col.names=T)#输出所有有表达的基因的logfc,pvalue等
#	write.table(et$table[et$table$PValue < pvalue & et$table$logFC < (-logfc),],file=paste(state[2],'high',state[1],'low',".txt",sep="_"),sep="\t",quote=F,row.names=T,col.names=T)
#	write.table(et$table[et$table$PValue < pvalue & et$table$logFC > (logfc),],file=paste(state[1],'high',state[2],'low',".txt",sep="_"),sep="\t",quote=F,row.names=T,col.names=T)
	de<-decideTestsDGE(et)
	detags<-rownames(y)[as.logical(de)]
	plotSmear(et,de.tags=detags,cex=0.5)#绘制差异表达基因
	abline(h=c(-logfc,logfc),col="blue")
	dev.off()
#	利用count计算rpkm,rpkm/2=fpkm(pair end,single end FPKM=RPKM)
#	y <- DGEList(counts=counts,genes=data.frame(Length=GeneLength))
#	y <- calcNormFactors(y)
#	RPKM <- rpkm(y)	
}

edger.fun(Tet2_countdata[,1:6],c('WT','Tet2_homo'),c(2,4),'WT_Tet2_homo',logfc=1)

##########################################
############   6.绘制hclust图   ###########
##########################################
'绘制hclust'
hcdup <- hclust(dist(tdata), "ave")
pdf('hcluster-dup.pdf')
plot(as.phylo(hcdup),type="unrooted",cex=0.5 )
dev.off()


#######################################
############   7.绘制PCA图   ###########
#######################################
'绘制PCA'
pdf("PCA-dup.pdf")
par(mar=c(5, 5, 3, 9))
paint <- rainbow(14)
pca <- prcomp(tdata)
colo <- c(paint[1],rep(paint[2:12],c(rep(2,11))),paint[13],rep(paint[14],4))
plot(pca$x[,1],pca$x[,2],pch="",main = "PCA of the MiNPC-dup",xlab = "PC1",ylab = "PC2")
text(pca$x[,1],pca$x[,2],labels=rownames(pca$x),col = colo, cex = 1)
usr <- par("usr")
x <- usr[2]*1.02
y <- usr[4]*1 
legend(x, y, rownames(pca$x), pch=c(16,16,16,16), col=colo, box.lty=0, cex=1, xpd=TRUE)
dev.off()

'top/bottom gene in PCA component score'
pca<-princomp(datanodup)
score_matrix<-data.frame(pca[6])
pca_comp1 <- score_matrix[order(score_matrix[,1]),]
pca_comp1_genem <- rownames(pca_comp1)[1:100]
pca_comp1_genep <- c(rownames(pca_comp1)[(length(pca_comp1[,1])-99):length(pca_comp1[,1])])
write(pca_comp1_genem,file='Top_100_genes_in_pca_component1_minus.txt')
write(pca_comp1_genep,file='Top_100_genes_in_pca_component1_plus.txt')

'3D PCA图的绘制'
pdf("PCA-dup-3D.pdf")
opar <- par(no.readonly=TRUE)
par(las=2)
paint <- rainbow(14)
pca <- prcomp(tdata)
colo <- c(paint[1],rep(paint[2:12],c(rep(2,11))),paint[13],rep(paint[14],4))
scatterplot3d(pca$x[,1],pca$x[,2],pca$x[,3],color=colo,pch=20,mar=c(5, 5, 5, 7),main = "PCA of the MiNPC-dup",xlab = "PC1",ylab = "PC2",zlab="PC3",angle=20)
usr <- par("usr")
x <- usr[2]*1.05
y <- usr[4]*1 
legend(x, y, rownames(pca$x), pch=c(16,16,16,16), col=colo, box.lty=0, cex=0.95, xpd=TRUE)
par(opar)
dev.off()


#################################################
############   8.利用heatmap2绘制热图   ###########
#################################################
'绘制各个样本间的correlation的heatmap图,并聚类'
library(gplots)
ColorRamp<-colorRampPalette(c("blue","white","red"), bias=1)(100)
datacor<-cor(datanodup, method="pearson", use="pairwise.complete.obs")
data.lab<-matrix(as.character(round(datacor, 2)), ncol=dim(datacor)[2])
pdf("Heatmap-com.pdf")
heatmap.2(datacor, col=ColorRamp,cellnote= data.lab, trace="none", cexRow=0.5, cexCol=0.5)
dev.off()


'对基因表达数据进行heatmap绘制,利用kmeans进行聚类'
'转化refseq到基因名'
#refseq2genenm <- function(vect,fpkm){
#	index <- c()
#	for (i in vect){index <- c(index,which(fpkm[,'tracking_id']==i))}
#	return(fpkm[index,'gene_id'])
#}

library(gplots)
set.seed(0928)
km = kmeans(filter_data,7)
kmdata =  cbind(filter_data,km$cluster)
kmdata1 = kmdata[kmdata[,length(kmdata[1,])]==1,]
kmdata2 = kmdata[kmdata[,length(kmdata[1,])]==2,]
kmdata3 = kmdata[kmdata[,length(kmdata[1,])]==3,]
kmdata4 = kmdata[kmdata[,length(kmdata[1,])]==4,]
kmdata5 = kmdata[kmdata[,length(kmdata[1,])]==5,]
kmdata6 = kmdata[kmdata[,length(kmdata[1,])]==6,]
kmdata7 = kmdata[kmdata[,length(kmdata[1,])]==7,]
write(rownames(kmdata1),file='filter_genes_km1.txt')
write(rownames(kmdata2),file='filter_genes_km2.txt')
write(rownames(kmdata3),file='filter_genes_km3.txt')
write(rownames(kmdata4),file='filter_genes_km4.txt')
write(rownames(kmdata5),file='filter_genes_km5.txt')
write(rownames(kmdata6),file='filter_genes_km6.txt')
write(rownames(kmdata7),file='filter_genes_km7.txt')
o = order(kmdata[,length(kmdata[1,])])
kmdata = kmdata[o,1:(length(kmdata[1,])-1)]
ColorRamp<-colorRampPalette(c("darkblue", "white", "red"), bias=1)(1000)
pdf('Kmeans_heatmap_Esrrb.pdf')
hm_file<-heatmap.2(kmdata, Colv=F, Rowv=F, col=ColorRamp, labRow=NA, labCol=colnames(kmdata),cexCol = 0.8, trace="none", scale="none")
title('Esrrb kmeans heatmap')
dev.off()


##############################################
############   9.利用image绘制热图   ###########
##############################################
'利用image进行heatmap图的绘制'
image_ranked_gene <- function(matrix, M, T, geneID)
{
	pdf(T)
	exprdata <- c()
	for (i in geneID){exprdata <- rbind(exprdata,matrix[which(matrix[1]==i),])}
	expr <- exprdata[,-1]
	temp <- c()
	for (i in 1:nrow(expr)){
		Normcol <- (expr[i,]-min(expr[i,]))/(max(expr[i,])-min(expr[i,]))
		temp <- rbind(temp,Normcol)
	}
	expr <- temp
	pdf(T,height=8,width=3.5)
#	par(oma=c(0.5,0.5,0.5,0.5),mar=c(2,5,3,3))
#	layout(matrix(seq(2),nrow=2),height=c(4,0.75))
	ColorRamp <- colorRampPalette(c("white","red"), bias=1)(10000)   #color list
	image(1:ncol(expr), 1:nrow(expr), t(expr), axes=F, col=ColorRamp, xlab="", ylab="",main=M)
	axis(side=1,1:ncol(expr),labels=colnames(expr),cex.axis=0.8,las=2)
	axis(side=2,1:nrow(expr),labels=exprdata[,1],cex.axis=0.8,las=2);box()
	dev.off()
}
image_ranked_gene(datanodup,"Endoderm-related-genes","Endoderm-related-genes.pdf",c("Sox7","Sox17","Afp"))

'利用image绘制heatmap复杂版'
pdf(file="iNPC_heatmap.pdf")
par(mar=c(4,3,3,8))
zmax <- quantile(newdata,0.95)
zmin <- quantile(newdata,0.05)
ColorRamp <- colorRampPalette(c("white","red"), bias=1)(10000)  
ColorLevels <- seq(to=zmax,from=zmin, length=10000) 
newdata[newdata<zmin] <- zmin
newdata[newdata>zmax] <- zmax
ColorRamp_ex <- ColorRamp[round( (min(newdata)-zmin)*10000/(zmax-zmin) ) : round( (max(newdata)-zmin)*10000/(zmax-zmin) )]
image(1:ncol(newdata), 1:nrow(newdata), t(newdata), axes=F, col=ColorRamp_ex, xlab="", ylab="",useRaster=T,main="heatmap of iNPC process")

axis(side=2)
axis(side=1,at=c(1,2,3,4,5),labels=c("MEF.P","Day1","Day14","P6","Brain-NPC"),las=2)
axis(side=4,at=c(
length(upmidmid)/2+0.5,
length(upmidmid)+length(midupup)/2+0.5,
length(upmidmid)+length(midupup)+length(midupmid)/2+0.5,
length(upmidmid)+length(midupup)+length(midupmid)+length(midmidup)/2+0.5,
length(upmidmid)+length(midupup)+length(midupmid)+length(midmidup)+length(downmidmid)/2+0.5,
length(upmidmid)+length(midupup)+length(midupmid)+length(midmidup)+length(downmidmid)+length(middowndown)/2+0.5,
length(upmidmid)+length(midupup)+length(midupmid)+length(midmidup)+length(downmidmid)+length(middowndown)+length(middownmid)/2+0.5,
length(upmidmid)+length(midupup)+length(midupmid)+length(midmidup)+length(downmidmid)+length(middowndown)+length(middownmid)+length(midmiddown)/2+0.5),las=2,
labels=c("upmidmid","midupup","midupmid","midmidup","downmidmid","middowndown","middownmid","midmiddown"))

abline(v=1+0.5)
abline(v=2+0.5)
abline(v=3+0.5)
abline(v=4+0.5)
abline(v=5+0.5)

abline(h=length(upmidmid))
abline(h=length(upmidmid)+length(midupup))
abline(h=length(upmidmid)+length(midupup)+length(midupmid))
abline(h=length(upmidmid)+length(midupup)+length(midupmid)+length(midmidup))
abline(h=length(upmidmid)+length(midupup)+length(midupmid)+length(midmidup)+length(downmidmid))
abline(h=length(upmidmid)+length(midupup)+length(midupmid)+length(midmidup)+length(downmidmid)+length(middowndown))
abline(h=length(upmidmid)+length(midupup)+length(midupmid)+length(midmidup)+length(downmidmid)+length(middowndown)+length(middownmid))
box()
dev.off()



################################################
############   10.利用smooth绘制散点图   ###########
################################################

'利用smooth scatter进行散点图的绘制'
corrplot <- function(dataframe,colnm1,colnm2){
data <- dataframe
logic <- is.na(data$col1) | is.na(data$col2)
x <- data$col1[!logic]
y <- data$col2[!logic]
pdf(paste(colnm1,colnm2,'correlation.pdf',sep="_"))
smoothScatter(x,y,xlab=colnm1,ylab=colnm2)
title(main=paste(colnm1,colnm2,'correlation plot',sep=" "))
legend("topleft",inset=0.05,paste("cor = ",round(cor(scale(data1[4]),scale(data2[4])),4),sep=""),bty = "n",cex=0.5)
dev.off()
}
corrplot(datanodup,'P6','Brain-NPC')



'绘制散点图并且绘制拟合的直线，并在拟合直线上下各加1条平行线'
pdf('WT_vs_Tet2_homo_KO.pdf')
plot(log(Tet2_countdata[!delete_index,'WT']+1,2),log(Tet2_countdata[!delete_index,'Tet2_homo']+1,2),
xlab='WT expression (log2(count))',ylab='Tet2 homo KO expression (log2(count))',col='grey',pch=16,main='WT vs Tet2_homo_KO')
model=lm(log(Tet2_countdata[!delete_index,'Tet2_homo']+1,2)~log(Tet2_countdata[!delete_index,'WT']+1,2))
abline(model)
abline(1.02846,0.92426,lty=2)
abline(-0.92846,0.92426,lty=2)
legend('topleft',paste('R2',round(cor(log(Tet2_countdata[!delete_index,'WT']+1,2),log(Tet2_countdata[!delete_index,'Tet2_homo']+1,2)),2)
,sep='='),bty='n')
dev.off()





#################################################
############   11.绘制火山图表示差异基因   ##########
#################################################

'绘制火山图展示差异基因'
setwd('/Volumes/X/2.Tet12_hmc_mc/2016.12.09(results)/1.edgeR/genesymbol')
all_WT_Tet2_hetero <- read.table('WT_Tet2_hetero.txt',header=T,sep='\t')
all_WT_Tet2_homo <- read.table('WT_Tet2_homo.txt',header=T,sep='\t')
Tet2_hetero_low <- WT_Tet2_hetero[WT_Tet2_hetero[,'logFC']>1 & WT_Tet2_hetero[,'PValue']<0.05,]
Tet2_hetero_up <- WT_Tet2_hetero[WT_Tet2_hetero[,'logFC']<(-1) & WT_Tet2_hetero[,'PValue']<0.05,]
Tet2_homo_low <- WT_Tet2_homo[WT_Tet2_homo[,'logFC']>1 & WT_Tet2_homo[,'PValue']<0.05,]
Tet2_homo_up <- WT_Tet2_homo[WT_Tet2_homo[,'logFC']<(-1) & WT_Tet2_homo[,'PValue']<0.05,]
dup_Tet2_hetero_low <- !duplicated(Tet2_hetero_low[,'genes'])
dup_Tet2_hetero_up <- !duplicated(Tet2_hetero_up[,'genes'])
dup_Tet2_homo_low <- !duplicated(Tet2_homo_low[,'genes'])
dup_Tet2_homo_up <- !duplicated(Tet2_homo_up[,'genes'])
Tet2_hetero_low <- Tet2_hetero_low[dup_Tet2_hetero_low,]
Tet2_hetero_up <- Tet2_hetero_up[dup_Tet2_hetero_up,]
Tet2_homo_low <- Tet2_homo_low[dup_Tet2_homo_low,]
Tet2_homo_up <- Tet2_homo_up[dup_Tet2_homo_up,]

Tet2_hetero_up <- nrow(Tet2_hetero_up)
Tet2_hetero_down <- nrow(Tet2_hetero_low)
all_WT_Tet2_hetero[all_WT_Tet2_hetero[,'PValue']<0.000001,'PValue'] <- 0.000001
pdf('WT_vs_Tet2_hetero_KO.pdf')
plot(all_WT_Tet2_hetero[,'logFC'],(-log(all_WT_Tet2_hetero[,'PValue'],10)),
xlab='Log2 fold change : WT/Tet2_hetero',ylab='-Log10(P-value)',col='lightblue',pch=16,main='WT vs Tet2_hetero_KO',xlim=c(-10,10))
legend(-10.5,5,paste('Tet2_hetero_up : ',Tet2_hetero_up,sep=''),bty='n')
legend(2.5,5,paste('Tet2_hetero_down : ',Tet2_hetero_down,sep=''),bty='n')
abline(v=1,lty=2)
abline(v=-1,lty=2)
abline(h=1.30103,lty=2)
dev.off()




#################################################
############   12.绘制venn图表示overlap   #########
#################################################
'利用两个数据集vectror,进行绘制韦恩图'
venn_plot <- function(file1,file2,pdf,label1,label2,overlapf1,overlapf2,overlapf3){
	library('VennDiagram')
	file1_dif <- as.character(file1)
	file2_dif <- as.character(file2)
	overlap <- intersect(file1_dif,file2_dif)
	pdf(pdf)
	venn.plot <- draw.pairwise.venn(length(file1_dif), length(file2_dif), length(overlap), c(label1, label2),col=c(rainbow(15)[9],'purple'),fill=c(rainbow(15)[9],'purple'),ext.text=T,cat.pos = c(350, 180),alpha = rep(0.5, 2))
	write.table(setdiff(file1_dif,overlap),overlapf1,quote=F,row.names=F,col.names=F)
	write.table(overlap,overlapf2,quote=F,row.names=F,col.names=F)
	write.table(setdiff(file2_dif,overlap),overlapf3,quote=F,row.names=F,col.names=F)
	dev.off()
}

venn_plot(Tet2_hetero_low,Tet2_homo_low,'Homo_Hetero_down_overlap.pdf','Tet2_hetero','Tet2_homo','Tet2_hetero.txt','Tet2_hetero_homo.txt','Tet2_homo.txt')




#####################################################
############   13.绘制折线图，修改坐标轴上的值   #########
#####################################################
'绘制折线图，修改坐标轴上的值'
plot(c(1,2,3,4),c(d0_matrix, d1_EPO_matrix, d3_EPO_matrix, d6_EPO_matrix),xlab='stage',ylab='dispersity',main='EPO_process',pch=16,type='b',col='red',xaxt='n')
axis(side=1,at=c(1,2,3,4),labels=c("d0","d1","d3","d6"))

'绘制折线图'
plot(c(1,4),c(-1.5,1.5),pch="",ylim=c(min(data_mean),max(data_mean)),xlim=c(1,4),xlab='',ylab='',xaxt='n')
for (i in 1:nrow(data_mean)){lines(seq(4),data_mean[i,],col="red",pch=16,type='l')}
axis(1,c(1,2,3,4),c('E14','E16','E18','AT2'))


##########################################
############   14.绘制堆叠条形图   #########
##########################################
'绘制堆叠条形图'
pdf('probability.pdf')
#par(mfrow=c(1,3))
#### TS
barplot(c(1,1),col='gold',main='TS',ylab='Probability',names.arg=c('A','B'))
barplot(c(0.46,0.51),col='darkgreen',add=TRUE)
abline(h=0.5,lty=3)


##########################################################################
############   15.计算correlation以及对应的p-value,并绘制线性拟合线   #########
##########################################################################
'计算correlation以及对应的p-value,并绘制线性拟合线'
cor_plot <- function(data_exp,genex,geney,pdfnm,color){
	pdf(pdfnm)
	plot(log(data_exp[,genex]+1,2),log(data_exp[,geney]+1,2),xlab=genex,ylab=geney,col=color,pch=16,main=paste('correlation between',genex,'and',geney,sep=' '))
	legend('topleft',paste('r',round(cor(log(data_exp[,genex]+1,2),log(data_exp[,geney]+1,2)),3),sep=' = '),bty='n')
	cortest <- cor.test(log(data_exp[,genex]+1,2),log(data_exp[,geney]+1,2), method = "pearson", alternative = "two.sided")
	if(cortest$p.value==0){
	legend('topright',paste('p-value',2* pt(cortest$statistic,  df = cortest$parameter, lower.tail=FALSE),sep=' = '),bty='n')}else {
	legend('topright',paste('p-value',cortest$p.value,sep=' = '),bty='n')}
	fit <- lm(log(data_exp[,geney]+1,2)~log(data_exp[,genex]+1,2))
	lines(log(data_exp[,genex]+1,2),fitted(fit),col='red') 	
	dev.off()
}

cor_plot(T_LIHC_exp,'FANCD2','PCNA','FANCD2_PCNA_correlation.pdf','gold')


##############################################
############   16.绘制带点的boxplot   #########
#############################################
'boxplot带点的那种'
boxplot(lis,ann=F,las=1,cex.axis=0.8,outpch=NA,frame.plot=F,ylab="log_rpkm",xaxt="n",main=gene,ylim=c(0,12))
  stripchart(lis,vertical = TRUE, method = "jitter", 
             cex=0.8,pch = 19, col = rep(color,7),add = TRUE) 

  axis(1,1:length(DAY),cex.axis=1,labels = DAY)

##############################################
############   17.删除inf，nan，na值   #########
##############################################
'删除inf，nan，na值'
data_scale <- as.data.frame(t(scale(t(data_mean))))
data_scale[data_scale==Inf|data_scale==-Inf] <- NA
na_index <- is.na(data_scale[,1])
data_scale <- data_scale[!na_index,]



#####################################
############   18.熵的计算   #########
#####################################
'计算entroy'
http://stackoverflow.com/questions/27254550/calculating-entropy
https://cran.r-project.org/web/packages/entropy/entropy.pdf

'分段计算熵的值,根据数据分成8段进行计算'
library("entropy")
data <- read.table('merge_fpkm.txt',header=T)
data <- data.frame(data$gene,(data$ESC_1.1+data$ESC_1.2)/2.0,(data$ESC_2.1+data$ESC_2.2)/2.0,(data$ESC_3.1+data$ESC_3.2)/2.0
,(data$ESC_4.1+data$ESC_4.2)/2.0,(data$ESC_5.1+data$ESC_5.2)/2.0,(data$ESC_6.1+data$ESC_6.2)/2.0,(data$ESC_7.1+data$ESC_7.2)/2.0
,(data$ESC_8.1+data$ESC_8.2)/2.0,(data$ESC_9.1+data$ESC_9.2)/2.0,(data$ESC_10.1+data$ESC_10.2)/2.0,(data$ESC_11.1+data$ESC_11.2)/2.0
,(data$ESC_12.1+data$ESC_12.2)/2.0,(data$ESC_13.1+data$ESC_13.2)/2.0,(data$ESC_14.1+data$ESC_14.2)/2.0,(data$ESC_15.1+data$ESC_15.2)/2.0
,(data$ESC_16.1+data$ESC_16.2)/2.0,(data$ESC_17.1+data$ESC_17.2)/2.0,(data$ESC_18.1+data$ESC_18.2)/2.0,(data$ESC_19.1+data$ESC_19.2)/2.0
,(data$ESC_20.1+data$ESC_20.2)/2.0,(data$ESC_21.1+data$ESC_21.2)/2.0,(data$ESC_22.1+data$ESC_22.2)/2.0,(data$ESC_23.1+data$ESC_23.2)/2.0
,(data$ESC_24.1+data$ESC_24.2)/2.0,(data$ESC_25.1+data$ESC_25.2)/2.0,(data$ESC_26.1+data$ESC_26.2)/2.0,(data$ESC_27.1+data$ESC_27.2)/2.0
,(data$ESC_28.1+data$ESC_28.2)/2.0,(data$ESC_29.1+data$ESC_29.2)/2.0,(data$ESC_30.1+data$ESC_30.2)/2.0,(data$ESC_31.1+data$ESC_31.2)/2.0
,(data$ESC_32.1+data$ESC_32.2)/2.0,(data$ESC_33.1+data$ESC_33.2)/2.0,(data$ESC_34.1+data$ESC_34.2)/2.0,(data$ESC_35.1+data$ESC_35.2)/2.0
,(data$ESC_36.1+data$ESC_36.2)/2.0,(data$ESC_37.1+data$ESC_37.2)/2.0,(data$ESC_38.1+data$ESC_38.2)/2.0,(data$ESC_39.1+data$ESC_39.2)/2.0
,(data$ESC_40.1+data$ESC_40.2)/2.0,(data$ESC_41.1+data$ESC_41.2)/2.0,(data$ESC_42.1+data$ESC_42.2)/2.0,(data$ESC_43.1+data$ESC_43.2)/2.0
,(data$ESC_44.1+data$ESC_44.2)/2.0,(data$ESC_45.1+data$ESC_45.2)/2.0,(data$ESC_46.1+data$ESC_46.2)/2.0,(data$ESC_47.1+data$ESC_47.2)/2.0
,(data$ESC_48.1+data$ESC_48.2)/2.0)
colnames(data) <- c('gene','ESC_1','ESC_2','ESC_3','ESC_4','ESC_5','ESC_6','ESC_7','ESC_8','ESC_9','ESC_10','ESC_11','ESC_12','ESC_13','ESC_14'
,'ESC_15','ESC_16','ESC_17','ESC_18','ESC_19','ESC_20','ESC_21','ESC_22','ESC_23','ESC_24','ESC_25','ESC_26','ESC_27','ESC_28','ESC_29','ESC_30','ESC_31'
,'ESC_32','ESC_33','ESC_34','ESC_35','ESC_36','ESC_37','ESC_38','ESC_39','ESC_40','ESC_41','ESC_42','ESC_43','ESC_44','ESC_45','ESC_46','ESC_47','ESC_48')

index <- c()
for (i in 1:nrow(data)){
if (sum(as.numeric(as.character(data[i,-1]))>=0.0001)>=24){
	index <- c(index,i)
	}
}
data <- data[index,]
entroy_vect <- c()
for (i in 1:nrow(data)){
x1 = as.numeric(as.character(data[i,-1]))
y1 = discretize(x1, numBins=8)
entroy_vect <- c(entroy_vect,entropy(y1, unit="log"))
}
data <- cbind(data,'entroy'=entroy_vect)
data <- data[order(data[,ncol(data)],decreasing=T),]
write.table(data,'merge_sc_fpkm_entroy.txt',quote=F,sep='\t',row.name=F,col.names=T)


##################################################################################
############   19.利用秩来计算两两基因之间的关系，欧氏距离越大，两个基因差的越远   #########
##################################################################################
'利用秩来计算两两基因之间的关系，欧氏距离越大，两个基因差的越远，一个高一个低'
data <- read.table('/Volumes/X/3.cell_decision/dataset/3.ESC_RNA-seq_data/parse_reslut/single_cell_ESC/merge_sc_fpkm_entroy.txt',header=T,row.names=1)
data <- data[,-49]
for (i in 1:nrow(data)){
	data[i,]=rank(data[i,])
}
dis <- dist(data,p=2)#europian distance
dis_matrix <- data.frame(as.matrix(dis))
gene_dis_rank_N <- colnames(dis_matrix)[order(dis_matrix['Nanog',],decreasing=T)];which(gene_dis_rank_N=='Gata6')#distance max rank
gene_dis_rank_S <- colnames(dis_matrix)[order(dis_matrix['Sox2',],decreasing=T)];which(gene_dis_rank_S=='Gata6')#distance max rank
gene_dis_rank_P <- colnames(dis_matrix)[order(dis_matrix['Pou5f1',],decreasing=T)];which(gene_dis_rank_P=='Gata6')#distance max rank

plot(seq(48),data['Nanog',],col='red',pch=16,type='b')
points(seq(48),data['Sdc4',],col='blue',pch=16,type='b')





#############################################
############   19.直方图模拟正态分布   #########
#############################################
'直方图模拟正态分布'
probability <- function(data,name){
	pdf(paste(name,'.pdf',sep=''))
	hist(data,freq=F,main=paste('Histgram of ',name,sep=''))
	curve(dnorm(x, mean = mean(data), sd = sd(data)), min(data), max(data), add = TRUE, col = "red", lwd = 2)
	pop_sd <- sd(data)*sqrt((length(data)-1)/(length(data)))
	pop_mean <- mean(data)
	legend('topleft',paste(paste(round(pop_mean,2)),' ± ',paste(round(pop_sd,2)),sep=''),bty='n')
	dev.off()
#	print(c(pop_mean,pop_sd))
}






#######################################################################################
############   20.利用kmeans进行绘制不同类别的曲线,以中位数为中线,上下是百分之40和60   #########
#######################################################################################
'利用kmeans进行绘制不同类别的曲线,以中位数为中线,上下是百分之40和60'
library(amap)

set.seed(1)
km <- Kmeans(diffexp_0_12H_median,k,iter.max=100,method = "correlation")
pdf(file = paste("Cor_kmeans_cluster_",k,".pdf",sep=""), width = 5, height = 3.5);
clusterCol <- rainbow(k)
for (each in seq(k)){
     genes = names(which(km$cluster==each))
     v1 = apply(diffexp_0_12H_median[genes,],2,median)
     v2 = apply(diffexp_0_12H_median[genes,],2,quantile,probs=c(0.4,0.6))[1,]
     v3 = apply(diffexp_0_12H_median[genes,],2,quantile,probs=c(0.4,0.6))[2,]
     plot(v1,lwd=3,type="l",col=clusterCol[each],main=paste("cluster_",each,"--",length(genes),"genes ",sep=""),ylim=c(min(c(v1,v2,v3)),max(c(v1,v2,v3))),xlab=NA,ylab="Averaged expression level(Log TPM)",xaxt="n")
     axis(side=1,1:length(time_num),names(time_num),las=2);axis(side=2);box()
     polygon(c(1,1:length(time_num),length(time_num):2),c(v2[1],v3,v2[length(time_num):2]),col=adjustcolor("grey", alpha.f = 0.4),border=NA)
#     plotMatrix <- rbind(plotMatrix,nData[names(which(km$cluster==each)),])
#     tmp_boundary_sum <- tmp_boundary_sum+length(which(km$cluster==each))
#     clusterBoundary <- c(clusterBoundary,tmp_boundary_sum)
#     fileName = paste(k,"cluster_",each,"_gene.txt",sep="")
#     write.table(cbind(modGenes),file=fileName,row.names=F,col.names=F,quote=F,sep="\t")
}
dev.off()





#################################################
##########    21.计算transition index   ##########
#################################################
'计算transition index'
#transition index mean(abs(cor(genes)))/mean(cor(cells))

#除去那些在n%个细胞都不表达的基因
rem <- function(x,n=1){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
  remove_logis <- !(r >= dim(x)[2]*n)
  return(x[remove_logis,])
}

#这里要改一下,去除batch effect
#http://blog.sina.com.cn/s/blog_72512a1d0102wo6r.html
#https://github.com/johnstantongeddes/sim-transcriptome/blob/master/TPM.R
norm <- function(x){
	return(x*1000000/x[length(x)])
}
TPM <- function(x){
  x <- as.matrix(x)
  sumx <- apply(x,2,sum)
  x <- rbind(x,sumx)
  x <- apply(x,2,norm)
  x <- x[-dim(x)[1],]
  return(x)
}

#计算transition index
transition_index <- function(data){
data_cell_cor <- ((sum(cor(data))-ncol(data))/2.0)/((ncol(data)^2-ncol(data))/2)
data_gene_cor <- ((sum(abs(cor(t(data))))-nrow(data))/2.0)/((nrow(data)^2-nrow(data))/2)
return(data_gene_cor/data_cell_cor)
}


rawdata <- read.csv('~/Desktop/hESC_differentiation/GSE75748_sc_time_course_ec.csv',header=T,row.names=1,sep=',')
data <- TPM(rawdata)
log2data <- log(data+1)
#不同时间点对应的样本数
time_num <- table(substr(colnames(data),4,6))
log2data_00h <- rem(log2data[,c(seq(92))])
log2data_12h <- rem(log2data[,c(seq(93,194))])
log2data_24h <- rem(log2data[,c(seq(195,260))])
log2data_36h <- rem(log2data[,c(seq(261,432))])
log2data_72h <- rem(log2data[,c(seq(433,570))])
log2data_96h <- rem(log2data[,c(seq(571,758))])

transition_index(log2data_00h)




###########################################
#######   22.dtwclust进行时序分类   #########
###########################################


















