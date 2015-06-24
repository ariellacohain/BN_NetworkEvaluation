rm(list=ls())
options(stringsAsFactors = F)

#setting work directory
setwd("~/Dropbox/Network Methods Evaluation/")
#loading Param 
load("SimulatedDataParameters.RData")
str(Param)

# dataName = Param$NETID[1]

for(dataName in Param$NETID){

  cat("data --",dataName,"\n")
  cat("run ",which(Param$NETID == dataName ), "out of",length(Param$NETID),"\n")
  
  expr = read.delim(paste("SimulatedData/",dataName,"/exp.txt",sep=""),header=F,sep=",")
  genePos = read.delim(paste("SimulatedData/",dataName,"/genepos.txt",sep=""),header=T,sep=",")
  network = read.delim(paste("SimulatedData/",dataName,"/network.txt",sep=""),header=F,sep=",")
  snp  = read.delim(paste("SimulatedData/",dataName,"/snp.txt",sep=""),header=F,sep=",")
  snpPos = read.delim(paste("SimulatedData/",dataName,"/snppos.txt",sep=""),header=T,sep=",")
  
  ## call eQTLs for simultation data: 
  library(MatrixEQTL)
  
  ## Setting up Matrix eQTL for running all linear models:
  # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  errorCovariance = numeric();
  cisDist = 1e6;
  pvOutputThreshold_cis = 0.05;
  pvOutputThreshold_tra = 0.05; 
  
  snp = as.matrix(snp)
  rownames(snp) = snpPos$name
  snpsME = SlicedData$new()  
  snpsME$CreateFromMatrix(snp)
  
  expr = as.matrix(expr)
  rownames(expr) = genePos$symbol
  genesME = SlicedData$new()  
  genesME$CreateFromMatrix(expr)
  
  ## setting up the snp positions and gene positions
  names(snpPos)=c("marker","pos","chr")
  snpPos$chr= as.numeric(snpPos$chr)
  snpPos$pos= as.numeric(snpPos$pos)
  
  genePos = genePos[,-c(2)]
  names(genePos) = c("gene","start_coord","end_coord","chr")
  # Running Matrix EQTL: 
  snpPos= snpPos[,c(1,3,2)]
  genePos = genePos[,c(1,4,2,3)]
  
  me = Matrix_eQTL_main(
    snps = snpsME, #matrix eqtl formmatted sliced data of snp of interest
    gene = genesME, #all genes in matrix eqtl format
    useModel = useModel, #linear model
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name     = NULL,#output_file_name_cis, #output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    output_file_name.cis = NULL,#output_file_name_tra,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpPos,
    genepos = genePos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
  );
  
  write.table(me$cis$eqtls, file = paste("SimulatedData/",dataName,"/cis.txt",sep=""), sep="\t",quote=F,row.names=F)
  write.table(me$trans$eqtls, file = paste("SimulatedData/",dataName,"/trans.txt",sep=""), sep="\t",quote=F,row.names=F)
  
}