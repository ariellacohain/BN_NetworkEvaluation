rm(list=ls())
options(stringsAsFactors = F)

args <- commandArgs(trailingOnly=TRUE)

workingDir = "/sc/orga/scratch/cohaia01/NetworkMethodsEvaluation" #"~/Dropbox/Network Methods Evaluation/"

dataName = args[[1]]

#setting work directory
setwd(workingDir)
#loading Param 
load("SimulatedDataParameters.RData")
# str(Param)
# dataName = "Simdata_0000000000000A" ### EDIT THIS OUT // JUST FOR TESTING

genePos = read.delim(paste("SimulatedData/",dataName,"/genepos.txt",sep=""),header=T,sep=",")

network = read.delim(paste("SimulatedData/",dataName,"/network.txt",sep=""),sep=',',header=F)
network$node1 = genePos$symbol[network[,1]]
network$node2 = genePos$symbol[network[,2]]
network$act_inhib = "act"
network$act_inhib[network[,3]== -1] = "inhib"

cit  = read.delim(file = paste("SimulatedData/",dataName,"/cit_results_filtered.txt",sep=""),sep=" ",header=T)

inCITresults <- function(cit, node1, node2){
	p_act = cit$p_cit[cit$cis_gene == node1 & cit$trans_gene == node2]
	fdr_act = cit$fdr[cit$cis_gene == node1 & cit$trans_gene == node2]
	
	p_inh = cit$p_cit[cit$cis_gene == node2 & cit$trans_gene == node1]
	fdr_inh = cit$fdr[cit$cis_gene == node2 & cit$trans_gene == node1]

	if(sum(cit$cis_gene == node1 & cit$trans_gene == node2)>0){
		return(cbind(p_act,fdr_act,"act"))
	}else if(sum(cit$cis_gene == node2 & cit$trans_gene == node1)>0){
		return(cbind(p_inh,fdr_inh,"inh"))
	}else{
		return(NA)
	}

}

pval = c()
fdr = c()
type = c()
for(i in 1:nrow(network)){

	inCit = inCITresults(cit,network$node1[i],network$node2[i])
	if(!is.na(a)){
		pval= c(pval, inCit[1])
		fdr = c(fdr, inCit[2])
		type = c(type, inCit[3])
	}else{
		pval =c(pval,"NA")
		fdr = c(fdr, "NA")
		type = c(type, "NA")
	}
	
}
network$pval = pval
network$fdr = fdr
network$type = type

write.table(network, file = paste("SimulatedData/",dataName,"/network_cit_added.txt",sep=""),sep="\t",row.names=F,col.names=T)