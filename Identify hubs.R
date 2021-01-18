#datE is the data frame of gene expression profiles across different rows(samples), each column represent a gene. 
#modLabel is a vector with the same length of datE, corresponding to the module assignments for each of gene in datE.
modList <- function(modLabel,datE){ 
mods <- unique(modLabel)
mlist <- vector("list", length=length(mods))
names(mlist) <- mods
for (i in 1: length(mods)){
mlist[[i]] <- modGenes(mods[i], modLabel, datE)
}
return (mlist)
}



mod_list <- modList(modLabel,datE)
#connectivity 

k <- function(x,A){
ax <- A[x,x]
conn <- sort(apply(ax,1,sum)-1, decreasing =T)
return(conn)
}


hub_k <- lapply(mod_list, k, A=A) # A is the adjacency matrix, same dimension with the number of genes

#Module membership:

MEs <- moduleEigengenes(datE, modLabel)$eigengenes
colnames(MEs) <- substring(colnames(MEs), 3)
mods <- unique(modLabel) #all the modules we have
mod_mem <- vector("list", length=length(mods)) 
names(mod_mem) <- mods

idxme <- match(names(mod_mem),colnames(MEs))
mes <- MEs[,idxme]#OK

for (i in 1: length(mod_mem)) {
exprm <- datE[,mod_list[[i]]]
mod_mem[[i]] <- ssort(bicor(exprm, mes[,i])[,1])
}
#calculation betweenness centraliy in weighted network
library(tnet)
bcw <- function(genes,adjMatrix,alphai){
a1 <- adjMatrix[genes,genes]
diag(a1) <- 0
net <- as.tnet(a1, type="weighted one-mode tnet")

# Get node "names" from the matrix
nodeLabels <- colnames(a1)

# Compute metrics
b1 <- data.frame(betweenness_w(net,alpha=alphai))
bw1 <- b1[,2]
names(bw1) <- nodeLabels#same with genes
return(bw1)
}

#weighted betweenness centrality
wbc <- lapply(mod_list, function(x) bcw(x, A, 1))
