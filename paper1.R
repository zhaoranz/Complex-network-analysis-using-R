options(stringsAsFactor=F)
library(WGCNA)
library(multtest)
library(tnet)
library(glmnet)
library(GSEABase)
library(GOstats)

#pre processing 

load("cpmALL.RData")
#1. remove masked genes
ga <- read.table("gene.annotator.txt", header=T)
gID <- ga$ID
gIDt <- transferidT(gID)
keep.genes <- gIDt[ga$masked==FALSE]
expr <- cpmALL[keep.genes,]

#2. extract samples from middle and moulting stages

sampo  <- colnames(expr)
mid1 <- grep("middle",sampo)
mid2 <- grep("pad1.m.",sampo)
cpm.mid<- expr[,c(mid1,mid2)]
old1 <- grep("old", sampo)
old2 <- grep("pad1.o", sampo)
mo1 <- grep("moult",sampo)
cpm.omo <- expr[,c(old1,old2, mo1)]
datE.momo <- cbind(cpm.mid, cpm.omo)
#have a look 
dim(datE.momo)
#[1] 9493   45

#3. set cutoff to filter genes with low expression: cpm should larger than 3 in at least 3 samples
f1 <- rowSums(datE.momo>3)
f2 <- f1>=3
datEM <- datE.momo[f2,]# 7108  45
middle <- t(datEM[,1:18])
omo <- t(datEM[,19:45])

# 4 check whether there are any outlying samples
A <- adjacency(t(omo), type="distance")#moulting
k <- as.numeric(apply(A,2,sum))-1
zk <- scale(k)
zk[zk<(-2.5)]#threshold is often set as -2.5
#no outlier for moulting samples,
## analoguous process for middle samples
A <- adjacency(t(middle), type="distance")#middle
k <- as.numeric(apply(A,2,sum))-1
zk <- scale(k)
zk[zk<(-2.5)]
#no outlier for middle samples,


#construct network
enableWGCNAThreads()
s.mid <- abs(bicor(middle))
s.olm <- abs(bicor(omo))
powers = c(c(1:10), seq(from = 12, to=20, by=2))
powers=seq(1:20)
sft <- pickSoftThreshold(middle, powerVector=powers, verbose=3)#estimated=7
sfto <- pickSoftThreshold(omo, powerVector=powers, verbose=3)#old+moulting, estimated=7
sftmo <- pickSoftThreshold(momo, powerVector=powers, verbose=3)#momo
sizeGrWindow(16, 9)
par(mfrow = c(2,2));

png("scale-free-parameter-estimate.png", width=2000, height=1125)
par(mfrow = c(3,2));

cex1 = 0.9; cex1=2
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence \n(Middle network)"), cex.lab=1.5, cex.main=cex1);
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity \n(Middle network)"), cex.lab=1.5, cex.main=cex1)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# Scale-free topology fit index as a function of the soft-thresholding power
plot(sfto$fitIndices[,1], -sign(sfto$fitIndices[,3])*sfto$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence \n(Moulting network)"), cex.lab=1.5, cex.main=cex1);
text(sfto$fitIndices[,1], -sign(sfto$fitIndices[,3])*sfto$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sfto$fitIndices[,1], sfto$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity \n(Moulting network)"), cex.lab=1.5, cex.main=cex1)
text(sfto$fitIndices[,1], sfto$fitIndices[,5], labels=powers, cex=cex1,col="red")



#global network
plot(sftmo$fitIndices[,1], -sign(sftmo$fitIndices[,3])*sftmo$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence \n(Global network)"), cex.lab=1.5, cex.main=cex1);
text(sftmo$fitIndices[,1], -sign(sftmo$fitIndices[,3])*sftmo$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sftmo$fitIndices[,1], sftmo$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity \n(Global network)"), cex.lab=1.5, cex.main=cex1)
text(sftmo$fitIndices[,1], sftmo$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

beta <- 7
a.mid <- s.mid^beta
a.olm <- s.olm^beta


#identify modules for network based on middle-stage samples, and draw dendrogram
w <- 1-a.mid
geneTree <- hclust(as.dist(w), method="average")
modLabel<- cutreeDynamic(dendro=geneTree, distM=w, deepSplit=4, pamRespectsDendro=FALSE,minClusterSize=30)#84, grey=203
mC.mid <- labels2colors(modLabel)#84
MEs.mid <- moduleEigengenes(middle, mC.mid)$eigengenes
sizeGrWindow(16,9)
png("dendro-middle.png", width=2000, height=1125)
plotDendroAndColors(geneTree, mC.mid, "Dynamic Tree Cut",cex.main=1.5,
dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors \n(middle network)")
dev.off()

#identify modules for network based on moult-stage samples, and draw dendrogram
w <- 1-a.olm
geneTree <- hclust(as.dist(w), method="average")
modLabel<- cutreeDynamic(dendro=geneTree, distM=w, deepSplit=4, pamRespectsDendro=FALSE,minClusterSize=30)#61, grey=444
mC.omo <- labels2colors(modLabel)#old+moult mod colors
MEs.olm <- moduleEigengenes(omo, mC.omo)$eigengenes

sizeGrWindow(16,9)
png("dendro-moulting.png", width=2000, height=1125)
plotDendroAndColors(geneTree, mC.omo, "Dynamic Tree Cut",cex.main=1.5,
dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors \n(moulting network)")
dev.off()

#global network dendrogram drawing
w <- 1- a.momo
geneTree <- hclust(as.dist(w), method="average")



png("dendro-global.png", width=2000, height=1125)
plotDendroAndColors(geneTree, mc.momo, "Dynamic Tree Cut",cex.main=1.5,
dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors \n(global network)")
dev.off()




# construct a GCN based on all samples from both middle stages and moulting stages
momo <- rbind(middle, omo)
cl <- rep(c(0,1), times=c(18,27))
s.momo <- abs(bicor(momo))
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(momo, powerVector=powers, verbose=3)#es



plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence \n(global)"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity \n(global)"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.print(pdf, file="scale-connectivity.pdf")
a.momo <- s.momo^7
w <- 1-a.momo
mesmomo <- moduleEigengenes(momo, mc.momo)$eigengene

#############################################module preservation analysis##############################################################################################
#examine whether modules found in middle stages can be also found in moultinig stages (preserved?)
etiketti <- c("middle", "moulting")
multiE <- vector("list",2)
multiE[[1]] <- list(data=as.data.frame(middle))
multiE[[2]] = list(data = as.data.frame(omo))
names(multiE) <- etiketti
multiColor <- list(middle=mC.mid)
mp.midomo  <- modulePreservation(multiE, multiColor,referenceNetworks = 1, nPermutations = 200,networkType = "unsigned", corFnc = "bicor",
                              randomSeed = 1, quickCor = 0, verbose = 4, checkData=FALSE)#may take some time


save(mp.midomo, file="New-data/mp-midomo.RData")

#same process for modules in moulting network
etiketti2 <- c("moulting", "middle")
multiE <- vector("list",2)
multiE[[1]] <- list(data=as.data.frame(omo))
multiE[[2]] = list(data = as.data.frame(middle))
names(multiE) <- etiketti2
#check: lapply(multiE, lapply, dim)
multiColor <- list(moulting=mC.omo)
mp.omomid  <- modulePreservation(multiE, multiColor,referenceNetworks = 1, nPermutations = 200,networkType = "unsigned", corFnc = "bicor",
                              randomSeed = 1, quickCor = 0, verbose = 4, checkData=FALSE)

save(mp.omomid, file="New-data/mp-omomid.RData")


##############################################DE analysis##############################################################################################
#DE
cl <- rep(c(0,1), times=c(18,27))
#check whether the expression profile of genes follow normal distribution:
shtp <- function(x){
shapiro.test(x)$p}
shapiro.p <- apply(momo, 2, shtp)
table(shapiro.p < 0.05)

mt.momo <- mt.maxT(t(momo),classlabel=cl, test="wilcoxon", side="abs",B=100000)
ord.momo <- order(mt.momo$index)
adjp.momo <- mt.momo$adjp[ord.momo]
names(adjp.momo) <- colnames(middle)#looks OK

s <- adjp.momo[adjp.momo<.05]#significantly DE
deg <- names(adjp.momo[adjp.momo<.05])#1364
tfomo <- intersect(tfgt, colnames(momo))
deg1 <- setdiff(deg,tfomo)#1343#?

DEexp.mid <- middle[,deg1]
DEexp.omo <- omo[,deg1]
emid <- apply(DEexp.mid,2,mean)
eomo <- apply(DEexp.omo,2,mean)
#why must use the relactant tfgt? try tfomo
r2tf.mid <- (s.mid[getIndex(tfomo, rownames(s.mid)), getIndex(deg1, rownames(s.mid))])^2
r2tf.omo <- (s.olm[getIndex(tfomo, rownames(s.olm)), getIndex(deg1, rownames(s.olm))])^2
> dim(r2tf.omo)
[1]  135 1343
e2m <- emid^2
e2o <- eomo^2 
RIF1 <- rep(NA,135)
for ( i in 1:135){
rmid <- r2tf.mid[i,]
romo <- r2tf.omo[i,]
RIF1[i] <- sum(sapply(deg1, function(x) romo[x]*e2o[x]-rmid[x]*e2m[x]))
}
names(RIF1) <- rownames(r2tf.mid)


RIF2 <- rep(NA,135)
ndeg1 <- length(deg1)
r.mid <- s.mid[getIndex(tfomo, rownames(s.mid)), getIndex(deg1, rownames(s.mid))]
r.omo <- s.olm[getIndex(tfomo, rownames(s.olm)), getIndex(deg1, rownames(s.olm))]
DW <- r.omo-r.mid#differentially wired
averageE <- apply(rbind(DEexp.mid,DEexp.omo),2,mean)
stat.deg1 <- stat.momo[deg1]
RIF2 <- apply(DW,1, function(x) mean(averageE*stat.deg1*(x)^2))
#############################################have a look at the distribution of DE genes in the all-sample network################################################
cv1 <- cv.glmnet(as.matrix(mesmomo), cl, family="binomial", type.measure="class",alpha=.5)
n0mods <- colnames(mesmomo)[which(coef(cv1, s="lambda.min")!=0)-1]
moultmods <- substring(colnames(mesmomo)[which(coef(cv1, s="lambda.min")>0)-1],3)
del <- vector("list", length=length(moultmods))
for (i in 1: length(del)){
del[[i]] <- modGenes(moultmods[i], mc.momo, momo)
}
names(del) <- moultmods
sapply(del, function(x) length(intersect(x,deg))/length(x))
n0m <- substring(n0mods,3)
n0l <- vector("list", length=length(n0m))
for (i in 1: length(n0l)){
n0l[[i]] <- modGenes(n0m[i], mc.momo, momo)
}
names(n0l) <- n0m
sapply(n0l, function(x) length(intersect(x,deg))/length(x))
length(intersect(deg, unlist(n0l)))/length(deg)

#blast against FlyBase and query in genomeRNAi
cn<-c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
r <- read.csv("GenomeRNAi_v16_FrequentHitters.txt",sep="\t",skip=11,header=F)#USE READ.CSV!!!
d <- read.table("licefly.txt",sep="\t")#blast lice protein sequence against fly sequences
colnames(d) <- cn

droso <- read.csv("GenomeRNAi_v16_Drosophila_melanogaster.txt", sep="\t", skip=19,header=T)
g2p <- read.csv("fbgn_fbtr_fbpp_fb_2017_05.tsv",header=T,sep="\t",skip=4)#maybe fbgn and fbpp

maxflyID <- function(gene,d){
gp <- transferid(gene)
hits <- d[grep(gp, d$qseqid),c(1,2,11,12)]
hb <- hits$bitscore
max.idx <- which(hb==max(hb))
fly <- hits[max.idx,2]
return(fly)
}#only use the max of max idx, 



gpList <- vector("list", length=length(genes))
names(gpList) <- genes
for (i in 1:length(genes)){
gene <- genes[i]
gpList[[i]] <- maxflyID(gene,d)}# got a list of homolog protein ID from FlyBase, for each sea louse protein sequence
#are there any protein sequences having no homologs?
idx0 <- sapply(gpList, length)==0# if a sequence has no homologs, it will be annotated as "nohits" in the final phenotypes list.
#check whether there are any NA)
h1 <- gpList[!idx0]
h1.na <- sapply(h1, function(x) any(is.na(x)))
h0 <- gpList[idx0]
#no NA
 


uniFB <- function(pp,d,g2p){
if (length(pp) ==0){ return("noid")}
else {
wc <- rep(NA, length=length(pp))
for ( i in 1: length(pp)){
gnid <- g2p[grep(pp[i], g2p$FlyBase_FBpp),1]
if(length(gnid)>0) wc[i] <- gnid
}
names(wc) <- pp
return(unique(wc))
}
}
homoG <- lapply(gpList, uniFB, d=d, g2p=g2p)#a list of homolog Gene ID from FlyBase, for each sea louse gene ID, is ltiu
idx1 <- sapply(homoG, length)==1
#check the list
homo1 <- homoG[idx1]
noid <- sapply(homoG, function(x) "noid" %in% x)
all.equal(names(noidl), names(h0))
#OK
idx.na.gene <- sapply(homoG, function(x)all(is.na(x))) #no gene id were found

allGn <- unique(unlist(homoG))
allGid <- allGn[!is.na(allGn)]
table(allGid %in% droso$Gene.ID) #237 cannot be found

feni <- function(px,droso){

if (length(px)==1) {
if(is.na(px)){feno<-NA}
else{
feno <- droso[grep(px, droso$Gene.ID),7]
feno <- feno[feno!="none"]
return(feno)}}

else{
wc <- rep(NA, 1)
for ( i in 1: length(px)){
if(is.na(px[i])){feno<-NA}
else{
feno <- droso[grep(px[i], droso$Gene.ID),7]
wct <- feno[feno!="none"]
wc <- c(wc,wct)
}}
return (wc[-1])
}
}

lfeno <- lapply(ltiu,feni,droso=droso)
phenoListori <- lapply(homoG, pheno, droso=droso)
phenoList <- phenoListori
phenoList[noid] <- "nohits"
phenoList[idx.na.gene] <- "noGeneID"
num.pheno <- sapply(phenoListori, length)
num.pheno[noid] <- "nohits"
num.pheno[idx.na.gene] <- "noGeneID"


#lice-lice blastp
ll <- read.table("licelice.txt")
names(ll) <- cn
homolice <- function(x){
homo1 <- unique(ll[grep(x, ll$qseqid),2])
homo <- homo1[!(homo1 %in% x)]
return(homo)}


homoLice <- lapply(gp, homolice)
names(homoLice) <- genes #colnames(omo)
save(gpList, homoG, phenoListori, num.pheno, homoLice, file="2blast.RData")
############whether phenotype or lethal phentoype is enriched in a module?#############################
mods <- unique(mC.omo)
lfenlen <- sapply(phenoListori, length)
fgenes <- names(lfenlen[lfenlen>0])
N <- length(phenoListori)
M <- sum(lfenlen>0)


enrichPheno <- function(mod){
modgenes <- modGenes(mod, mC.omo, omo)
m <- length(modgenes)
k <- sum(modgenes %in% fgenes)
p <- phyper(k,M, N-M, m, lower.tail=F)
return(p)}
mod.p <- sapply(mods, enrichPheno) 
letha <- function(x){
y <- grep("lethal", x)
if (length(y) >0) return(TRUE)
else return(FALSE)}
lethal <- sapply(phenoListori, letha)
lethali <- lethal[lethal==TRUE]
M2 <- length(lethali)#778#, note, M changed
lethalg<- names(lethali)


#or change M to M2 here
enrichPhenoL <- function(mod){
modgenes <- modGenes(mod, mC.omo, omo)
m <- length(modgenes)
k <- sum(modgenes %in% lethalg)
p <- phyper(k,M2, N-M2, m, lower.tail=F)
return(p)
}

mod.pl <- sapply(mods, enrichPhenoL)#14 lethal enriched
################for middle network###########################################################
mmods <- unique(mC.mid)

enrichPheno <- function(mod){
modgenes <- modGenes(mod, mC.mid, middle)
m <- length(modgenes)
k <- sum(modgenes %in% fgenes)
p <- phyper(k,M, N-M, m, lower.tail=F)
return(p)}

enrichPhenoL <- function(mod){
modgenes <- modGenes(mod, mC.mid, middle)
m <- length(modgenes)
k <- sum(modgenes %in% lethalg)
p <- phyper(k,M2, N-M2, m, lower.tail=F)
return(p)
}

mid.p <- sapply(mmods, enrichPheno)
mid.pl <- sapply(mmods, enrichPhenoL)
save(mod.p, mod.pl, mid.p, mid.pl, file="Enrich-feno-p.RData")


########################GO enrichment analysis###################################################

genes <- colnames(omo);genep <- transferid(genes)
universeP <- genep[genep%in% slimFrame$gene_id]#4657
omo.alldfList <- vector("list",length=length(unique(mC.omo)))#61 mods
names(omo.alldfList)<- unique(mC.omo)
for ( i in 1:length(omo.alldfList)){
  modName <- names(omo.alldfList)[i]
  omo.alldfList[[i]]<- IdentifyGOP(modName, universeP,mC.omo, omo)
}

mid.alldfList <- vector("list",length=length(unique(mC.mid)))#84 mods
names(mid.alldfList)<- unique(mC.mid)
for ( i in 1:length(mid.alldfList)){
  modName <- names(mid.alldfList)[i]
  mid.alldfList[[i]]<- IdentifyGOP(modName, universeP,mC.mid, middle)
}

momo.alldfList <- vector("list",length=length(unique(mc.momo)))
names(momo.alldfList)<- unique(mc.momo)
for ( i in 1:length(momo.alldfList)){
  modName <- names(momo.alldfList)[i]
  momo.alldfList[[i]]<- IdentifyGOP(modName, universeP,mc.momo, momo)
}

##########################Final candidate selection table#########################################
dr <- BGSKC("darkred", mC.omo, omo, a.olm, a.mid, MEs.olm)
TF <- names(RIF1)
homolice <- sapply(homoLice, length)
others <- function(moddf){
g <- rownames(moddf)
moddf$phenoN <- num.pheno[g]
moddf$DE <- g %in% deg
moddf$TF <- g %in% TF
moddf$Lethal <- lethal[g]
moddf$homoLice <- homolice[g]
return(moddf)
}

dr <- others(dr)
write.table(dr, file="darkred.txt", quote=F, sep="\t")

interest.mods <- c("violet","darkred", "yellowgreen", "lavenderblush3")
for (i in 1:length(interest.mods)){
filn <- paste(interest.mods[i], "txt", sep=".")
mdf <- BGSKC(interest.mods[i], mC.omo, omo, a.olm, a.mid, MEs.olm)
mdf.new <- others(mdf)
write.table(mdf.new, file=filn, quote=F, sep="\t")

steeldf <- BGSKC("steelblue", mC.mid,middle, a.mid, a.olm, MEs.mid)
steeldf.new <- others(steeldf)



















fillist <- function(x){if(length(x)==0) {x="noid"}
else (x=x)}
lfifeno <- lapply(lfeno, fillist)#this is useful only when dcal.
########################################################final candidate selection results##############################################################################
candi <- c(4777,12631,5000,0591,1767,8887,6670,1059,7048,0043,5129,1007,5299,2021,3179,3849,8045,1144,2499,1436)
candidate <- sprintf("EMLSAT000000%05d",candi)
