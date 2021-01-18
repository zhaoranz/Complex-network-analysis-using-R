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
