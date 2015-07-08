# cd /sc/orga/scratch/cohaia01/NetworkMethodsEvaluation; for i in `ls SimulatedData/*/cit_results.RData |cut -f2 -d'/'`; do echo $i; Rscript ~/thirdparty/BN_NetworkEvaluation/plot_ROC.R $i 0.05; done
rm(list=ls())
library(ROCR)
library(reshape)

options(stringsAsFactors=F)
setwd("/sc/orga/scratch/cohaia01/NetworkMethodsEvaluation")

files = list.files("SimulatedData/",pattern ="cit_results.RData",recursive=T)
Names = sapply(files,function(X){strsplit(X,"/",fixed =T)[[1]][1]})
# dataName = args[[1]]
# fdr = as.numeric(as.character(args[[2]]))
fdr = 0.05

#####--------FUNCTIONS START-------------
fill_matrix <- function( positions, matrix){
	if(ncol(positions) == 3){
		for(i in 1:nrow(positions)){
			matrix[positions[i,2],positions[i,1]]  = positions[i,3]
		}
	}else{
		for(i in 1:nrow(positions)){
			matrix[positions[i,2],positions[i,1]]  = 1 
		}
	}
	
	return(matrix)

}

plot_ROC <- function(scores, values, color="black",title="ROCR"){
	pred <- prediction(scores,values)
	perf <- performance(pred,"tpr","fpr")
	plot(perf,col=color,lty=3, lwd=3,main = title)

	abline(0,1,col="red")
	auc <- performance(pred,"auc") 

	legend("bottomright",paste("AUC: ",round(auc@y.values[[1]],3)))
	
	return(auc@y.values)
}
#####---------FUNCTIONS END------------
all_auc = c()
pass_auc = c()
top_auc = c()
no_fdr = c()

for(dataName in Names){
	print(dataName)
	cat("--on run", which(Names == dataName),"out of",length(Names),"\n")
	load(paste("SimulatedData/",dataName,"/cit_results.RData",sep=""))
	if(sum(colnames(output)=="fdr")==0){
		cat(dataName,"had no fdr calc \n")
		no_fdr = c(no_fdr,dataName)
	}else{

		genePos = read.delim(paste("SimulatedData/",dataName,"/genepos.txt",sep=""),header=T,sep=",")

		network = read.delim(paste("SimulatedData/",dataName,"/network.txt",sep=""),sep=',',header=F)
		network = read.delim(paste("SimulatedData/",dataName,"/network.txt",sep=""),sep=',',header=F)
		network$node1 = genePos$symbol[network[,1]]
		network$node2 = genePos$symbol[network[,2]]
		network$act_inhib = "act"
		network$act_inhib[network[,3]== -1] = "inhib"

		output$pos1 = sapply(output$cis_gene, function(X){which(genePos$symbol == X)})
		output$pos2 = sapply(output$trans_gene, function(X){which(genePos$symbol == X)})



		real_network_matrix = matrix(0,nrow= nrow(genePos),ncol= nrow(genePos))
		colnames(real_network_matrix) = genePos$symbol
		rownames(real_network_matrix) = genePos$symbol

		real_network_matrix = fill_matrix(network[,c(1,2)],real_network_matrix)

		top_data = output[output$passSnpThresh == TRUE & output$fdr <= fdr & output$p_cit_reactive >0.05,]

		cit_network = matrix(0,nrow= nrow(genePos),ncol= nrow(genePos))
		colnames(cit_network) = genePos$symbol
		rownames(cit_network) = genePos$symbol


		cit_all = fill_matrix(output[output$p_TassocG_pre <0.05,c('pos1','pos2','fdr')],cit_network)
		cit_passSNP  = fill_matrix(output[output$passSnpThresh == TRUE,c('pos1','pos2','fdr')],cit_network)
		cit_topData = fill_matrix(top_data[,c('pos1','pos2')],cit_network)

		m_cit_all = melt(cit_all)
		m_cit_passSNP = melt(cit_passSNP)
		m_cit_topData = melt(cit_topData)
		m_real = melt(real_network_matrix)

		pdf(paste("SimulatedData/",dataName,"/ROC.pdf",sep=""))
		par(mfrow= c(1,3))
		all = plot_ROC(m_cit_all$value,m_real$value,title="All CIT Results")
		pass = plot_ROC(m_cit_passSNP$value,m_real$value,title = "Pass SNP threshold Results")
		top = plot_ROC(m_cit_topData$value,m_real$value,title = paste("FDR <",fdr*100,"%"))
		dev.off()

		all_auc= c(all_auc,all)
		pass_auc= c(pass_auc,pass)
		top_auc= c(top_auc,top)
	}	
}

# Simdata_1000100010000A 
# which(top_auc == max(unlist(top_auc)))



# dataName =  "Simdata_1000110000100C"
# fdr = 0.05



# # calculating the values for ROC curve
# pred <- prediction(m_cit_all$value,m_real$value)
# perf <- performance(pred,"tpr","fpr")
# # changing params for the ROC plot - width, etc

# par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex.axis=1.3,cex.lab=1.4)
# # plotting the ROC curve
# plot(perf,col="black",lty=3, lwd=3)
# points()
# # calculating AUC
# auc <- performance(pred,"auc") 
# auc
# abline(0,1,col="red")
# legend("bottomright",paste("AUC: ",round(auc@y.values[[1]],3)))
# dev.off()


# # save(top_data,auc,fdr,file = paste("SimulatedData/",dataName,"/ROC.RData",sep="")