# cd /sc/orga/scratch/cohaia01/NetworkMethodsEvaluation; for data in `ls SimulatedData/ `; do echo $data;Rscript ~/thirdparty/BN_NetworkEvaluation/CIT_synthetic_data.R $data; done
# cd /sc/orga/scratch/cohaia01/NetworkMethodsEvaluation; for data in `ls SimulatedData/ `; do echo $data;submitjob 2 -c 10 -m 2 -q low -P acc_STARNET Rscript ~/thirdparty/BN_NetworkEvaluation/CIT_synthetic_data.R $data; done
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

if_ran = try(load(paste("SimulatedData/",dataName,"/cit_results.RData",sep="")))
if(!inherits(if_ran,"try-error")){
  cat("RAN ALREADY!")
  quit()
}

library(citpp)
library(foreach)
library(multicore)
library(doMC)
options(cores = multicore:::detectCores())
registerDoMC(threads) 

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
#sig_cis = cis[cis$FDR <= cis_fdr_cutoff, ]
#sig_cis = cis[cis$FDR <= cis_fdr_cutoff, ]


### --- FUNCTIONS --- #####

getConditionalResults<-function(cis_eqtl_data,snpNameInterest,geneME = geneME,snpME = snpInterest){
  # cis_eqtl_data = me_cis[me_cis$gene==cis.gene,]
  # snpInterest = paste(snp.of.interest[2],snp.of.interest[3],sep=":")
  
  if(sum(cis_eqtl_data$snps == snpNameInterest)==0){
    cat("GWAS SNP NOT PRESENT\n")
    return(NA)
  }else{
    cis_gene  = unique(cis_eqtl_data$gene)
    gwas_pval = cis_eqtl_data$pvalue[cis_eqtl_data$snp==snpNameInterest]
    most_sigSnp = cis_eqtl_data$snps[cis_eqtl_data$pvalue==min(cis_eqtl_data$pvalue,na.rm=T)]
    gene_expr = unlist(geneME[,colnames(geneME)== cis_gene]) #unlist(geneME$FindRow(cis_gene)$row)[1,]
    gwas_snp_expr = unlist(snpME[,colnames(snpME)==snpNameInterest]) #unlist(snpME$FindRow(snpNameInterest)$row)[1,]
    most_sig_snp_expr = unlist(snpME[,colnames(snpME)==most_sigSnp[1]]) #unlist(snpME$FindRow(most_sigSnp[1])$row)[1,]
    # if(is.element(most_sigSnp,snpNameInterest)){
    #     res = cbind(cis_gene,gwas_pval, snpNameInterest, most_sigSnp[1],min(cis_eqtl_data$pvalue,na.rm=T), 0,1)
    #     colnames(res) = c("cis_gene","gwas_cis_pvalue","gwas_snp","most_sig_snp","sig_snp_pval","cond_pval","cond_rsqd")
    #     return(res)
    #   }else{
    lm_conditional = lm(gene_expr~most_sig_snp_expr,na.action = na.exclude)
    resid.vec = naresid(lm_conditional$na.action, lm_conditional$residuals)
    a = try(lm(resid.vec[!is.na(resid.vec) & !is.na(gwas_snp_expr) & resid.vec != "X" & gwas_snp_expr != "X"]~as.numeric(gwas_snp_expr[!is.na(resid.vec) & !is.na(gwas_snp_expr) & resid.vec != "X" & gwas_snp_expr != "X"])), TRUE)  
    
    if (!inherits(a, "try-error"))
    {
      a = summary(a)
      # cat(counter, " ", esnps.mod[counter,"gene"], "\t", esnps.mod[counter,"dbsnp_id"], "\t", a$r.squared, "\t", 1-pf(a$fstatistic, a$df[1]-1, a$df[2])[1], "\n", sep="")
      cond_rsqd = a$r.squared
      cond_pval = 1-pf(a$fstatistic, a$df[1]-1, a$df[2])[1]
    }else{
      cond_pval = NA
      cond_rsqd = NA
    }
    
    lm_conditional_MaxGvnGwas = lm(gene_expr~gwas_snp_expr,na.action = na.exclude)
    resid.vec_MaxGvnGwas = naresid(lm_conditional_MaxGvnGwas$na.action, lm_conditional_MaxGvnGwas$residuals)
    a_MaxGvnGwas = try(lm(resid.vec_MaxGvnGwas[!is.na(resid.vec_MaxGvnGwas) & !is.na(most_sig_snp_expr) & resid.vec_MaxGvnGwas != "X" & most_sig_snp_expr != "X"]~as.numeric(most_sig_snp_expr[!is.na(resid.vec_MaxGvnGwas) & !is.na(most_sig_snp_expr) & resid.vec_MaxGvnGwas != "X" & most_sig_snp_expr != "X"])), TRUE)  
    
    if (!inherits(a_MaxGvnGwas, "try-error"))
    {
      a_MaxGvnGwas = summary(a_MaxGvnGwas)
      # cat(counter, " ", esnps.mod[counter,"gene"], "\t", esnps.mod[counter,"dbsnp_id"], "\t", a$r.squared, "\t", 1-pf(a$fstatistic, a$df[1]-1, a$df[2])[1], "\n", sep="")
      cond_rsqd_MaxGvnGwas = a_MaxGvnGwas$r.squared
      cond_pval_MaxGvnGwas = 1-pf(a_MaxGvnGwas$fstatistic, a_MaxGvnGwas$df[1]-1, a_MaxGvnGwas$df[2])[1]
    }else{
      cond_pval_MaxGvnGwas = NA
      cond_rsqd_MaxGvnGwas = NA
    }
    
    
    res = cbind(cis_gene, gwas_pval, snpNameInterest, most_sigSnp[1],min(cis_eqtl_data$pvalue,na.rm=T), cond_pval,cond_rsqd, cond_pval_MaxGvnGwas,cond_rsqd_MaxGvnGwas)
    colnames(res) = c("cis_gene","gwas_cis_pvalue","snp","most_sig_snp","sig_snp_pval","cond_pval","cond_rsqd","cond_pval_MaxGvnGwas","cond_rsqd_MaxGvnGwas")
    return(res)
  }
  
}

### ---------------- #####

cond_results_all = c()
correlations = c()
all_trios = c()
for(snp_id in unique(sig_cis$snps)){
  #checking for every snp that has a cis eQTL, how many genes are called cis, and how many trans
  cis_genes = sig_cis$gene[sig_cis$snps == snp_id]
  trans_genes = sig_trans$gene[sig_trans$snps == snp_id]
  if(which( unique(sig_cis$snps) == snp_id)%%10 ==0){
    cat("--RUN",which(snp_id == unique(sig_cis$snps)),"out of",length(unique(sig_cis$snps)),"--\n")
    cat("---",snp_id,"number cis:",length(cis_genes),"number trans:",length(trans_genes),"---\n" )
      
  }

  
  #if SNP has at least 1 cis and at least 1 trans genes continue
  if(length(cis_genes)>0 & length(trans_genes)>0){
    #running conditional analysis for every cis gene that is present:
      cond_res = foreach(i  = 1:length(cis_genes),combine=rbind)%dopar%{
        cis.gene = cis_genes[i]
        # cat("on run ",i," / ",length(cis_genes),"\n")
        cisQTL_data = sig_cis[sig_cis$gene == cis.gene,]
        getConditionalResults(cisQTL_data,snp_id,expr,snp)  
      }
    
    cond_results = do.call(rbind,cond_res)
    if(nrow(cond_results) == 1){
      cond_results = data.frame(t(cond_results[!is.na(cond_results[,1]),]))
      
    }else{
      cond_results = data.frame(cond_results[!is.na(cond_results[,1]),])
    }
    cond_results$gwas_cis_pvalue = as.numeric(cond_results$gwas_cis_pvalue)
    cond_results$sig_snp_pval = as.numeric(cond_results$sig_snp_pval)
    cond_results$cond_rsqd = as.numeric(cond_results$cond_rsqd)
    cond_results$cond_pval = as.numeric(cond_results$cond_pval)
    cond_results$cond_rsqd_MaxGvnGwas = as.numeric(cond_results$cond_rsqd_MaxGvnGwas)
    cond_results$cond_pval_MaxGvnGwas = as.numeric(cond_results$cond_pval_MaxGvnGwas)
  
    # calculate trans cis correlation:
    cis_data = as.matrix(expr[,is.element(colnames(expr),cis_genes)])
    if(length(cis_genes)==1){
      colnames(cis_data) = cis_genes
    }
    trans_data = as.matrix(expr[,is.element(colnames(expr),trans_genes)])
    if(length(trans_genes)==1){
      colnames(trans_data) = trans_genes
    }
    for(i in 1:ncol(cis_data)){
      g.vec = unlist(cis_data[,i])
      cor = apply(trans_data, 2, FUN=function(e,g)
      {
        a = summary(lm(as.numeric(e[!is.na(e) & !is.na(g) & e != "X" & g != "X"])~as.numeric(g[!is.na(e) & !is.na(g) & e != "X" & g != "X"])))
        return(as.list(c(a$r.squared,1-pf(a$fstatistic, a$df[1]-1, a$df[2])[1])))
      }, g.vec)
      cor.rsqr = unlist(lapply(cor, FUN=function(x) { return(x[[1]]) }))
      cor.pv = unlist(lapply(cor, FUN=function(x) { return(x[[2]]) }))
      correlations=rbind(correlations,cbind(colnames(cis_data)[i],colnames(trans_data),cor.rsqr,cor.pv))
    } # end of correlation for loop
    
    trio = list(cis = cis_genes, trans = trans_genes)
    all_trios[[snp_id]] = trio
    cond_results_all = rbind(cond_results_all,cond_results)
    
    } # end of cis and trans at least 1 hit if statement 
  

} # end of loop over all snps with cis QTL


correlations = as.data.frame(unique(correlations))
names(correlations)=c("Cis_gene","Trans_gene","cor.rsqr","cor.pv")
cat("correlations finished.\n")
correlations$cor.pv = as.numeric(correlations$cor.pv)


## setting CITPP indexes and which tests to run:
snp_position = seq(1:ncol(snp))
names(snp_position) = colnames(snp)

expr_position = seq(1:ncol(expr))
names(expr_position) = colnames(expr)

L = snp
Cis = expr
Trans = expr

trios_tmp = sapply(names(all_trios),function(x){
  snp_col = snp_position[x]
  # cat(names(x))
  cis_position = expr_position[all_trios[[x]]$cis]
  trans_position = expr_position[all_trios[[x]]$trans]
  cbind(snp_col,cis_position,trans_position)
  })

trios = data.frame(do.call(rbind,trios_tmp))

trios$cor_pval = unlist(lapply(seq(1:nrow(trios)),function(x){
  correlations$cor.pv[correlations$Cis_gene == colnames(expr)[trios$cis_position[x]] & correlations$Trans_gene==colnames(expr)[trios$trans_position[x]]]
  }))
trios$cor_rsqr = unlist(lapply(seq(1:nrow(trios)),function(x){
  correlations$cor.rsqr[correlations$Cis_gene == colnames(expr)[trios$cis_position[x]] & correlations$Trans_gene==colnames(expr)[trios$trans_position[x]]]
  }))

cat("number of trios (pre filtering):",nrow(trios),"\n")
if(sum(trios$cor_pval <= cor_pval_cutoff)>0){
  trios = trios[trios$cor_pval <= cor_pval_cutoff,]
  cat("removing correlation high trios")
}else{
  cat("all trios removed\n")
  write.table("NO_TRIOS_LEFT",file=paste("SimulatedData/",dataName,"/NOcit_results.RData",sep=""))
  q()
}

cat("number of trios (post filtering):",nrow(trios),"\n")
output = cit(L,Cis,Trans, trios[,c(1,2,3)], threads = threads_cit)

cat("cit finished.\n")

cit_LTC = cit(L,Trans,Cis,trios[,c(1,3,2)],threads=threads_cit)
colnames(cit_LTC)= paste(colnames(cit_LTC),"reactive",sep="_")
cat("running cit L->T->C finished.\n")
output= cbind(output,cit_LTC)

## CHECK THAT THIS IS CORRECT || SPOT CHECK!!
#output = cit
output$snp_id = colnames(L)[output$L_index] #rep(snp.of.interest[1],nrow(cit))
output$cis_gene = colnames(Cis)[output$G_index]
output$trans_gene = colnames(Trans)[output$T_index]

output$max_snp = sapply(1:nrow(output),function(x){
  cond_results_all$most_sig_snp[
    cond_results_all$snp == output$snp_id[x] &
    cond_results_all$cis_gene == output$cis_gene[x]
  ] })

output$max_snp_pval = sapply(1:nrow(output),function(x){
  cond_results_all$sig_snp_pval[
    cond_results_all$snp == output$snp_id[x] &
    cond_results_all$cis_gene == output$cis_gene[x]
  ] })

output$cis_pvalue = sapply(1:nrow(output),function(x){
  cond_results_all$gwas_cis_pvalue[
    cond_results_all$snp == output$snp_id[x] &
    cond_results_all$cis_gene == output$cis_gene[x]
  ] })


output$cond_pval = sapply(1:nrow(output),function(x){
  cond_results_all$cond_pval[
    cond_results_all$snp == output$snp_id[x] &
    cond_results_all$cis_gene == output$cis_gene[x]
  ] })


output$cond_rsqd = sapply(1:nrow(output),function(x){
  cond_results_all$cond_rsqd[
    cond_results_all$snp == output$snp_id[x] &
    cond_results_all$cis_gene == output$cis_gene[x]
  ] })

output$cond_rsqd_MaxGvnGwas = sapply(1:nrow(output),function(x){
  cond_results_all$cond_rsqd_MaxGvnGwas[
    cond_results_all$snp == output$snp_id[x] &
    cond_results_all$cis_gene == output$cis_gene[x]
  ] })

output$cond_pval_MaxGvnGwas = sapply(1:nrow(output),function(x){
  cond_results_all$cond_pval_MaxGvnGwas[
    cond_results_all$snp == output$snp_id[x] &
    cond_results_all$cis_gene == output$cis_gene[x]
  ] })

output$rsqr_TassocG = trios$cor_rsqr
output$p_TassocG_pre = trios$cor_pval

output$cisQTL_fdr =sapply(1:nrow(output),function(x){
  sig_cis$FDR[
    sig_cis$snps == output$snp_id[x] &
    sig_cis$gene == output$cis_gene[x]
  ] })
output$transQTL_fdr = sapply(1:nrow(output),function(x){
  sig_trans$FDR[
    sig_trans$snps == output$snp_id[x] &
    sig_trans$gene == output$trans_gene[x]
  ] })

save(L,Cis,Trans,trios,output, file=paste("SimulatedData/",dataName,"/Pre_cit_results.RData",sep=""))

# cat("running REACTIVE cit (L->T->C)...\n")
# cit_LTC = cit(L,Trans,Cis,trios[,c(1,3,2)],threads=threads)
# colnames(cit_LTC)= paste(colnames(cit_LTC),"reactive",sep="_")
# cat("running cit L->T->C finished.\n")
# output= cbind(output,cit_LTC)

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
  if(i %% 10 ==0 ){

    cat("on run",i, "out of",nrow(output),"\n")
    as.numeric(print(Sys.time()))

  }
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

