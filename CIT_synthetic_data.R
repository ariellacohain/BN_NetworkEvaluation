rm(list=ls())
options(stringsAsFactors = F)

threads = 1
threads_cit = 10
dataName = c() #Param$NETID[1]
cis_fdr_cutoff = 0.20
trans_fdr_cutoff = 0.20
cor_pval_cutoff = 0.05

library(citpp)
library(foreach)
library(multicore)
library(doMC)
options(cores = multicore:::detectCores())
registerDoMC(threads) 

#setting work directory
setwd("~/Dropbox/Network Methods Evaluation/")
#loading Param 
load("SimulatedDataParameters.RData")
str(Param)
dataName = Param$NETID[1] ### EDIT THIS OUT // JUST FOR TESTING

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
for(snp_id in unique(sig_cis$snps)){
  #checking for every snp that has a cis eQTL, how many genes are called cis, and how many trans
  cis_genes = sig_cis$gene[sig_cis$snps == snp_id]
  trans_genes = sig_trans$gene[sig_trans$snps == snp_id]
  cat("---",snp_id,"number cis:",length(cis_genes),"number trans:",length(trans_genes),"---\n" )
  #if SNP has at least 1 cis and at least 1 trans genes continue
  if(length(cis_genes)>0 & length(trans_genes)>0){
    #running conditional analysis for every cis gene that is present:
      cond_res = foreach(i  = 1:length(cis_genes),combine=rbind)%dopar%{
        cis.gene = cis_genes[i]
        cat("on run ",i," / ",length(cis_genes),"\n")
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
    
    } # end of cis and trans at least 1 hit if statement 
  cond_results_all = rbind(cond_results_all,cond_results)

} # end of loop over all snps with cis QTL


correlations = as.data.frame(unique(correlations))
names(correlations)=c("Cis_gene","Trans_gene","cor.rsqr","cor.pv")
cat("correlations finished.\n")
correlations$cor.pv = as.numeric(correlations$cor.pv)


## setting CITPP indexes and which tests to run:

L = 
Cis = 
Trans = 

trios = 

output = cit(L,Cis,Trans, trios, threads = threads_cit)

cat("cit finished.\n")

cit_LTC = cit(L,Trans,Cis,trios[,c(1,3,2)],threads=threads)
colnames(cit_LTC)= paste(colnames(cit_LTC),"reactive",sep="_")
cat("running cit L->T->C finished.\n")
output= cbind(output,cit_LTC)

## CHECK THAT THIS IS CORRECT || SPOT CHECK!!
#output = cit
output$snp_id = colnames(L)[output$L_index] #rep(snp.of.interest[1],nrow(cit))
output$cis_gene = colnames(Cis)[output$G_index]
output$trans_gene = colnames(Trans)[output$T_index]

output$max_snp_chr_pos = cond_results$most_sig_snp[output$G_index]
output$max_snp_pval = cond_results$sig_snp_pval[output$G_index]
output$cis_pvalue = cond_results$gwas_cis_pvalue[output$G_index]
output$cond_pval = cond_results$cond_pval[output$G_index]
output$cond_rsqd = cond_results$cond_rsqd[output$G_index]

output$cond_rsqd_MaxGvnGwas = cond_results$cond_rsqd_MaxGvnGwas[output$G_index]
output$cond_pval_MaxGvnGwas = cond_results$cond_pval_MaxGvnGwas[output$G_index]

output$rsqr_TassocG = unlist(lapply(seq(1:nrow(output)),function(x){correlations$cor.rsqr[correlations$Cis_gene == output$cis_gene[x] & correlations$Trans_gene==output$trans_gene[x]]}))
output$p_TassocG_pre = unlist(lapply(seq(1:nrow(output)),function(x){correlations$cor.pv[correlations$Cis_gene == output$cis_gene[x] & correlations$Trans_gene==output$trans_gene[x]]}))

output$cisQTL_fdr = 
output$transQTL_fdr = 

cat("running REACTIVE cit (L->T->C)...\n")
cit_LTC = cit(L,Trans,Cis,trios[,c(1,3,2)],threads=threads)
colnames(cit_LTC)= paste(colnames(cit_LTC),"reactive",sep="_")
cat("running cit L->T->C finished.\n")
output= cbind(output,cit_LTC)

save(L,Cis,Trans,trios,output, file=paste("SimulatedData/",dataName,"/cit_results.RData",sep=""))

#testing commit
### NEXT STEP IS TO THEN RUN PERMUTATIONS 
