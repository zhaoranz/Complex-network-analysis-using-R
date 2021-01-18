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


save(mp.midomo, file="mp-midomo.RData")

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

save(mp.omomid, file="mp-omomid.RData")

