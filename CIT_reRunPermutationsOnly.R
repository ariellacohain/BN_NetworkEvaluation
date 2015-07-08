# cd /sc/orga/scratch/cohaia01/NetworkMethodsEvaluation; for data in `ls SimulatedData/ `; do echo $data;Rscript ~/thirdparty/BN_NetworkEvaluation/CIT_reRunPermutationsOnly.R $data; done
# cd /sc/orga/scratch/cohaia01/NetworkMethodsEvaluation; for data in `ls SimulatedData/ `; do echo $data;submitjob 4 -c 10 -m 2 -q low -P acc_STARNET Rscript ~/thirdparty/BN_NetworkEvaluation/CIT_reRunPermutationsOnly.R $data; done

rm(list=ls())
options(stringsAsFactors = F)

args <- commandArgs(trailingOnly=TRUE)

workingDir = "/sc/orga/scratch/cohaia01/NetworkMethodsEvaluation" #"~/Dropbox/Network Methods Evaluation/"

dataName = args[[1]]

threads = 5
threads_cit = 10
perm_threads = 10
numPermTests = 1000

cis_fdr_cutoff = 0.20
trans_fdr_cutoff = 0.20
cor_pval_cutoff = 0.05

library(citpp)
library(foreach)
library(multicore)
library(doMC)
options(cores = multicore:::detectCores())
registerDoMC(threads) 


if_ran = try(load(paste("SimulatedData/",dataName,"/cit_results.RData",sep="")))

if_ran_pre = try(load(paste("SimulatedData/",dataName,"/Pre_cit_results.RData",sep="")))

if(!inherits(if_ran,"try-error")){
  cat("RAN ALREADY!")

  #setting work directory
  setwd(workingDir)
  #loading Param 
  load("SimulatedDataParameters.RData")
  # str(Param)
  # dataName = Param$NETID[1] ### EDIT THIS OUT // JUST FOR TESTING

  cat("data --",dataName,"\n")
  cat("run ",which(Param$NETID == dataName ), "out of",length(Param$NETID),"\n")

  expr = read.delim(paste("SimulatedData/",dataName,"/exp.txt",sep=""),header=F,sep=",")
  genePos = read.delim(paste("SimulatedData/",dataName,"/genepos.txt",sep=""),header=T,sep=",")
  expr = t(expr)
  colnames(expr) = genePos$symbol

  snp  = read.delim(paste("SimulatedData/",dataName,"/snp.txt",sep=""),header=F,sep=",")
  snpPos = read.delim(paste("SimulatedData/",dataName,"/snppos.txt",sep=""),header=T,sep=",")
  snp = t(snp)
  colnames(snp) = snpPos$name

  cis = read.delim(paste("SimulatedData/",dataName,"/cis.txt",sep=""))
  trans = read.delim(paste("SimulatedData/",dataName,"/trans.txt",sep=""))

  sig_cis = cis
  sig_trans = trans


  ### NEXT STEP IS TO THEN RUN PERMUTATIONS 
  permute <- function(snp,gene1,gene2,ExprData, time = Sys.time(), maxit=50000,  n=10000,threads){
    set.seed(as.numeric(time))

    # rand_snps=sapply(1:n,function(X){snp[sample(seq(1:length(snp)))]})
    pos_gene2 = which(colnames(ExprData)==gene2)
    trans = ExprData[,pos_gene2]
    L = sapply(snp, round) #round or ceiling

    shuffTrans = sapply(1:n, function(x){
      res = rep(Inf, length(trans))
      for(j in unique(L)){
        res[L==j] = sample(trans[L==j])
      }
      return(res)
      })

    # trans_gene = sapply(1:n,function(X){ExprData[sample(seq(1:nrow(geneExpr))),pos_gene2]})

    l_index = 1
    c_index = rep(which(colnames(ExprData)== gene1),n) 
    t_index = 1:n

    trios = cbind(l_index,c_index,t_index)
    cit = cit(as.matrix(L),ExprData,shuffTrans,trios,maxit=maxit,threads=threads)

    return(cit)

  }

  perm_pvals = c()
  time = as.numeric(Sys.time())
  for(i in 1:nrow(output)){
    # if(i %% 10 ==0 ){

      cat("on run",i, "out of",nrow(output),"\n")
      print(Sys.time())

    # }
    snp_i = snp[,c(output$snp_id[i])]
    gene1 = output$cis_gene[i]
    gene2 = output$trans_gene[i]

    #run permutation test on snp, gene1, gene2 - returns matrix from cit with N rows corresponding to the number permutations
    permuted_res <- permute(snp_i,gene1,gene2,expr,time = time, maxit =50000, n = numPermTests,threads=perm_threads)

    # calculating how many permuted tests are more significant than observed (b+1)/(N+1) where b = # permuted < observed
    permuted_pval = (sum(permuted_res$p_cit < output$p_cit[i])+1)/ (numPermTests  +1)
    perm_pvals = c(perm_pvals,permuted_pval)
  }

  output$perm_pval = perm_pvals
  output$fdr = p.adjust(perm_pvals,"fdr")

  output$passSnpThresh = output$cond_pval <=0.05 | (output$cond_pval>0.05 & output$cond_pval_MaxGvnGwas >0.05) | output$snp_id == output$max_snp

  save(L,Cis,Trans,trios,output,time, file=paste("SimulatedData/",dataName,"/cit_results.RData",sep=""))

  sig_output = output[output$p_cit_reactive >0.05 & output$passSnpThresh, ]
  write.table(sig_output,file = paste("SimulatedData/",dataName,"/cit_results_filtered.txt",sep=""))



}else if(!inherits(if_ran_pre,"try-error")){
   cat("RAN ALREADY!\t too many trios -- keeping cis <10e-4\n")

  #setting work directory
  setwd(workingDir)
  #loading Param 
  load("SimulatedDataParameters.RData")
  # str(Param)
  # dataName = Param$NETID[1] ### EDIT THIS OUT // JUST FOR TESTING

  cat("data --",dataName,"\n")
  cat("run ",which(Param$NETID == dataName ), "out of",length(Param$NETID),"\n")

  expr = read.delim(paste("SimulatedData/",dataName,"/exp.txt",sep=""),header=F,sep=",")
  genePos = read.delim(paste("SimulatedData/",dataName,"/genepos.txt",sep=""),header=T,sep=",")
  expr = t(expr)
  colnames(expr) = genePos$symbol

  snp  = read.delim(paste("SimulatedData/",dataName,"/snp.txt",sep=""),header=F,sep=",")
  snpPos = read.delim(paste("SimulatedData/",dataName,"/snppos.txt",sep=""),header=T,sep=",")
  snp = t(snp)
  colnames(snp) = snpPos$name

  cis = read.delim(paste("SimulatedData/",dataName,"/cis.txt",sep=""))
  trans = read.delim(paste("SimulatedData/",dataName,"/trans.txt",sep=""))

  sig_cis = cis
  sig_trans = trans


  ### NEXT STEP IS TO THEN RUN PERMUTATIONS 
  permute <- function(snp,gene1,gene2,ExprData, time = Sys.time(), maxit=50000,  n=10000,threads){
    set.seed(as.numeric(time))

    # rand_snps=sapply(1:n,function(X){snp[sample(seq(1:length(snp)))]})
    pos_gene2 = which(colnames(ExprData)==gene2)
    trans = ExprData[,pos_gene2]
    L = sapply(snp, round) #round or ceiling

    shuffTrans = sapply(1:n, function(x){
      res = rep(Inf, length(trans))
      for(j in unique(L)){
        res[L==j] = sample(trans[L==j])
      }
      return(res)
      })

    # trans_gene = sapply(1:n,function(X){ExprData[sample(seq(1:nrow(geneExpr))),pos_gene2]})

    l_index = 1
    c_index = rep(which(colnames(ExprData)== gene1),n) 
    t_index = 1:n

    trios = cbind(l_index,c_index,t_index)
    cit = cit(as.matrix(L),ExprData,shuffTrans,trios,maxit=maxit,threads=threads)

    return(cit)

  }
  cat("before filtering number of trios:",nrow(output),"\n")
  
  output =   output[output$cis_pvalue <=10e-4,]
  
  cat("after filtering number of trios:",nrow(output),"\n")

  perm_pvals = c()
  time = as.numeric(Sys.time())
  for(i in 1:nrow(output)){
    # if(i %% 10 ==0 ){

      cat("on run",i, "out of",nrow(output),"\n")
      print(Sys.time())

    # }
    snp_i = snp[,c(output$snp_id[i])]
    gene1 = output$cis_gene[i]
    gene2 = output$trans_gene[i]

    #run permutation test on snp, gene1, gene2 - returns matrix from cit with N rows corresponding to the number permutations
    permuted_res <- permute(snp_i,gene1,gene2,expr,time = time, maxit =50000, n = numPermTests,threads=perm_threads)

    # calculating how many permuted tests are more significant than observed (b+1)/(N+1) where b = # permuted < observed
    permuted_pval = (sum(permuted_res$p_cit < output$p_cit[i])+1)/ (numPermTests  +1)
    perm_pvals = c(perm_pvals,permuted_pval)
  }

  output$perm_pval = perm_pvals
  output$fdr = p.adjust(perm_pvals,"fdr")

  output$passSnpThresh = output$cond_pval <=0.05 | (output$cond_pval>0.05 & output$cond_pval_MaxGvnGwas >0.05) | output$snp_id == output$max_snp

  save(L,Cis,Trans,trios,output,time, file=paste("SimulatedData/",dataName,"/cit_results.RData",sep=""))

  sig_output = output[output$p_cit_reactive >0.05 & output$passSnpThresh, ]
  write.table(sig_output,file = paste("SimulatedData/",dataName,"/cit_results_filtered.txt",sep=""))

  
}else{
  cat("data NOT RUN BEFORE --",dataName,"\n")
  q()
}
