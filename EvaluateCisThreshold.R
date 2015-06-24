rm(list=ls())
options(stringsAsFactors = F)
library(igraph)

#setting work directory
setwd("~/Dropbox/Network Methods Evaluation/")
#loading Param 
load("SimulatedDataParameters.RData")
str(Param)

##### FUNCTIONS: ######
getParentNodes <- function(network){
  # given the network file, returns vector of nodes w/ 0, 1 if gene is a parent
  # parent is defined as only source, cant be a sink. Activation and Inactivation not taken into accout. 
  
  nodes = unique(c(network$node1,network$node2))
  parent = c()
  for(node in nodes){
    is_source = is.element(node,network$node1) 
    is_sink = is.element(node,network$node2)
    if(is_source & !is_sink){
      parent = c(parent,node)
    }
  }
  is_parent = rep(0, length(nodes))
  names(is_parent) = nodes
  is_parent[parent] = 1
  return(is_parent)
  #g <-graph.empty() + vertices(unique(c(network$node1,network$node2)))
  #for(i in 1:nrow(network)){g <- g + edge(network$node1[i],network$node2[i])}
  #for(node in )
}
fdr_cutoffs = c(0.01,0.05,0.10,0.15)
all_res = c()
for(dataName in Param$NETID){
  cat("data --",dataName,"\n")
  cat("run ",which(Param$NETID == dataName ), "out of",length(Param$NETID),"\n")
  
  network = read.delim(paste("SimulatedData/",dataName,"/network.txt",sep=""),header=F,sep=",")
  cis = read.delim(paste("SimulatedData/",dataName,"/cis.txt",sep=""))
  genePos = read.delim(paste("SimulatedData/",dataName,"/genepos.txt",sep=""),header=T,sep=",")
  
  colnames(network)=c("source","sink","act_inhib")
  network$node1 = genePos$symbol[network$source]
  network$node2 = genePos$symbol[network$sink]
  
  is_parent_node = getParentNodes(network)
  
  fdr_parents =  sapply(fdr_cutoffs, function(x){
      num_parents_inFDR = sum(is.element(unique(cis$gene[cis$FDR<x]),names(is_parent_node[is_parent_node==1])))
      num_FDR_genes  = length(unique(cis$gene[cis$FDR<x]))
      c(x,num_parents_inFDR,num_FDR_genes,num_parents_inFDR/num_FDR_genes,num_parents_inFDR/sum(is_parent_node==1))
    })
  fdr_parents = as.data.frame(t(fdr_parents))
  names(fdr_parents) = c("fdr","parents_in_FDR","total_FDR_genes","percentage_of_FDR_genes","percentage_of_parents")
  fdr_parents$NETID = dataName
  all_res = rbind(all_res,fdr_parents)
  
}

Param2 = merge(Param,all_res)

save(Param2,file="Parameters_cisQTL_FDR_added.RData")
summary  = sapply(fdr_cutoffs,function(x){print(x); print(summary(Param2$percentage_of_FDR_genes[Param2$fdr == x ]))})
colnames(summary) = fdr_cutoffs
summary

summary  = sapply(fdr_cutoffs,function(x){print(x); print(summary(Param2$percentage_of_parents[Param2$fdr == x ]))})
colnames(summary) = fdr_cutoffs
summary
