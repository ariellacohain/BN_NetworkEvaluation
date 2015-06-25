# cd /sc/orga/scratch/cohaia01/NetworkMethodsEvaluation; for i in `ls SimulatedData/*/cit_results_filtered.txt |cut -f2 -d'/'`; do Rscript ~/thirdparty/BN_NetworkEvaluation/CIT_accuracy.R $i; done

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
	p_act_cit = min(cit$p_cit[cit$cis_gene == node1 & cit$trans_gene == node2])
	fdr_act = min(cit$fdr[cit$cis_gene == node1 & cit$trans_gene == node2])
	p_act_inh_cit = min(cit$p_cit_reactive[cit$cis_gene == node1 & cit$trans_gene == node2])
	bonfPval_act = min(cit$bonfPval[cit$cis_gene == node1 & cit$trans_gene == node2])

	p_inh_cit = min(cit$p_cit[cit$cis_gene == node2 & cit$trans_gene == node1])
	fdr_inh = min(cit$fdr[cit$cis_gene == node2 & cit$trans_gene == node1])
	p_inh_act_cit = min(cit$p_cit_reactive[cit$cis_gene == node2 & cit$trans_gene == node1])
	bonfPval_inh = min(cit$bonfPval[cit$cis_gene == node2 & cit$trans_gene == node1])
	return(cbind(p_act_cit,fdr_act,p_act_inh_cit,bonfPval_act,p_inh_cit,fdr_inh,p_inh_act_cit,bonfPval_inh))

	# if(sum(cit$cis_gene == node1 & cit$trans_gene == node2)>0){
	# 	return(cbind(p_act,fdr_act,p_act_inh,"act"))
	# }else if(sum(cit$cis_gene == node2 & cit$trans_gene == node1)>0){
	# 	return(cbind(p_inh,fdr_inh,p_inh_act,"inh"))
	# }else{
	# 	return(NA)
	# }

}

cit$bonfPval = p.adjust(cit$p_cit)
write.table(cit,file = paste("SimulatedData/",dataName,"/cit_results_filtered.txt",sep=""),sep=" ")
# pval = c()
# fdr = c()
# p_rev = c()
# type = c()
res=c()
for(i in 1:nrow(network)){

	inCit = inCITresults(cit,network$node1[i],network$node2[i])
	res = rbind(res,inCit)
	# if(!is.na(inCit)){
	# 	pval= c(pval, inCit[1])
	# 	fdr = c(fdr, inCit[2])
	# 	p_rev = c(p_rev, inCit[3])
	# 	type = c(type, inCit[4])
	# }else{
	# 	pval =c(pval,"NA")
	# 	fdr = c(fdr, "NA")
	# 	type = c(type, "NA")
	# }
	
}

network = cbind(network, res)
# network$pval = pval
# network$fdr = fdr
# network$type = type

write.table(network, file = paste("SimulatedData/",dataName,"/network_cit_added.txt",sep=""),sep="\t",row.names=F,col.names=T)