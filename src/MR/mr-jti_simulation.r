
#ml GCC/8.2.0  OpenMPI/3.1.4 Intel/2019.1.144  IntelMPI/2018.4.274 R/3.6.0

#MR-JTI simulation

args<-as.numeric(commandArgs(TRUE))

library(data.table)
library(foreach)
library(doParallel)
library(dplyr)
library(MendelianRandomization)
library(glmnet)
library(HDCI)

source('/home/z*/script/gtex/github/MR-JTI.function')

alpha_list = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
alpha = alpha_list[args]

#load snp in dosage

#load LD-pruned snp list
snp_pruned<-read.table(paste0('/data/c*/z*/data/gtex/geno/v8/liftover/pruning/prune0.1.prune.in'),stringsAsFactors = F)

#load empirical eqtl
emp_eqtl = data.frame(fread('/data/c*/z*/data/gtex/eqtl_v8/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz'))
emp_eqtl$gene_id = sapply(emp_eqtl$gene_id, function(x) strsplit(x,"[.]")[[1]][1])

#get the gene list (100 genes)
gene_list_geno = sub('.rds','',dir('/data/c*/z*/data/gtex/geno/v8/gene/dosage_1m/'))
gene_list_exp = sub('.rds','',dir('/data/c*/z*/data/gtex/exp/v8/weighted/Whole_Blood/'))
gene_list_all = intersect(intersect(gene_list_exp,gene_list_geno),
                          unique(emp_eqtl$gene_id))
set.seed(2021)
gene_list = sample(gene_list_all,100,replace = F)


#load ldsc
ldsc<-read.table('/data/c*/z*/data/gtex/geno/v8/ldsc/gtex_pruned.score.ld',header = T,stringsAsFactors = F)
ldsc<-ldsc[,c(1,8)]

#register for parallele running
registerDoParallel(cores = detectCores())

#for test only
i = 1

####parallel running####
output <- foreach(i = 1:length(gene_list), .combine=rbind) %dopar% {
    
  gene_id = gene_list[i]
  
  #load dosage (GTEx)
  dosage = readRDS(paste0('/data/c*/z*/data/gtex/geno/v8/gene/dosage_1m/',gene_id,'.rds'))
  dosage = dosage[,c(1,which(colnames(dosage) %in% snp_pruned[,1]))]
  dosage[,-1] = apply(dosage[,-1], MARGIN = 2, function(x) scale(x))
  
  #load expression (GTEx)
  exp<-readRDS(paste0('/data/c*/z*/data/gtex/exp/v8/weighted/Whole_Blood/',gene_id,'.rds')) #covariates adjusted expression level
  exp<-exp[exp$tissue=='Whole_Blood',] #use WB
  exp$sampleid = sub('GTEX.','',exp$sampleid)
  
  eqtl = merge(exp[,c('sampleid','exp')], dosage, by = 1)
  
  #estimate empirical eqtl effect size
  res = data.frame(beta_eqtl =NA, se_eqtl = NA, z_eqtl = NA, p_eqtl = NA)
  for(j in 1:(ncol(eqtl)-2)){
    rsid = colnames(eqtl)[j+2]
    
    ans = summary(lm(eqtl[,'exp']~eqtl[,rsid]))
    res[j,] = ans$coefficients[2,]
  }
  
  res$rsid = colnames(eqtl)[-c(1,2)]
  
  dosage = dosage[,-1] #rm the sample id column
  
  
  #generate the genotype for a large number of samples (as GWAS dataset)
  n_gwas = 50000
  dosage_gwas = dosage
  for(k in 1:ncol(dosage_gwas)){
    dosage_gwas[1:n_gwas,k] = sample(dosage[,k], n_gwas, replace = T)
  }
  
  #generate expression levels 
  sigma_e = sqrt(1-sum(res$beta_eqtl^2))   #matching the empirical data. var(expression) = 1
  set.seed(i)
  exp_gwas = as.matrix(dosage_gwas[,colnames(dosage)]) %*% res$beta_eqtl + rnorm(nrow(dosage_gwas), mean = 0, sd = sigma_e)
  
  #estimate beta_eqtl
  set.seed(i+100)
  sample_for_est_eqtl = sample(n_gwas, nrow(exp), replace = F)
  df_eqtl = data.frame(beta_eqtl =NA, se_eqtl = NA, z_eqtl = NA, p_eqtl = NA)
  for(j in 1:ncol(dosage_gwas)){
    rsid = colnames(dosage_gwas)[j]
    
    ans = summary(lm(exp_gwas[sample_for_est_eqtl] ~ dosage_gwas[sample_for_est_eqtl,rsid]))
    df_eqtl[j,] = ans$coefficients[2,]
  }
  
  
  #assume that 20% of the eQTLs (IVs) have pleiotropy effect
  n_pleiotropy = round(length(res$rsid)/5) #no of IVs have pleiotropy effect
  
  #generate pleiotropy effect
  if(n_pleiotropy == 0){
    pleiotropy_effect = rep(0,nrow(dosage_gwas))
  }else{
    pos_pleiotropy = sample(ncol(dosage),
                            n_pleiotropy,
                            replace = T)
    dosage_pleiotropy = dosage_gwas[,pos_pleiotropy]
    
    #assuming that the 'extra-effect is 2 times than its original size'
    pleiotropy_effect = as.matrix(dosage_pleiotropy) %*% abs(res$beta_eqtl[pos_pleiotropy]) * alpha*2
  }
  
  #generate the outcome
  set.seed(i*100 + args)
  outcome = alpha * exp_gwas + pleiotropy_effect + rnorm(nrow(dosage_gwas),0,1)
  
  #estimate beta_gwas 
  df_outcome = data.frame(beta_outcome = NA, se_outcome = NA, z_outcome = NA, p_outcome = NA)
  
  for(j in 1:ncol(dosage_gwas)){
    rsid = colnames(dosage_gwas)[j]
    
    ans = summary(lm(outcome ~ dosage_gwas[,rsid]))
    df_outcome[j,] = ans$coefficients[2,]
  }
  
  #combine beta_eqtl beta_gwas
  df = as.data.frame(cbind(df_eqtl,df_outcome))
  df$SNP = colnames(dosage_gwas)
  
  
  #--MR-IVW--
  
  #output
  ans = list()
  
  #orientation
  df$beta_outcome = ifelse(df$beta_eqtl>0,df$beta_outcome,df$beta_outcome*-1)
  df$beta_eqtl = abs(df$beta_eqtl)
  
  #regular ivw-mr
  res = mr_ivw(mr_input(bx = df$beta_eqtl, bxse = df$se_eqtl, by = df$beta_outcome, byse = df$se_outcome))
  ans$beta_ivw = res@Estimate
  ans$p_ivw = res@Pvalue
  
  #--MR-JTI--
  
  #merge with ldsc
  df<-merge(df,ldsc,by = 'SNP')
  
  #no of sig gwas loci (only consider these for pleiotropy control)
  n_gwas_sig<-length(which(df$p_outcome<0.05))
  
  if(n_gwas_sig ==0){
    #in case there is no IV significantly associated with the outcome, reduce to ivw
    ans$beta_median_mrjti = res@Estimate
    ans$beta_mean_mrjti = res@Estimate
    ans$p_mrjti = res@Pvalue
  
  #skip genes with less than 20 obs (no way for cv)
  }else if(nrow(df)>20){
  
  #weighted by the precision
  #weights=1/(df$se_outcome*df$se_eqtl)^2
  
  #penalty factor, no penalty for the first two features (beta_eqtl and ldsc)
  penalty.factor<-c(0,0,rep(1,n_gwas_sig))
  
  #include the identity matrix for IVs with potential pleiotropy effects
  df<-df[order(df$p_outcome),]
  df<-df[,c('SNP','beta_outcome','beta_eqtl','ldscore')]
  if(n_gwas_sig>0){
    df[,(ncol(df)+1):(ncol(df)+n_gwas_sig)]<-diag(nrow(df))[,1:n_gwas_sig]
  }
  
  #---assign x and y---
  y=df[,'beta_outcome'];x=as.matrix(df[,3:ncol(df)])
  #scaling
  #y=scale(y)
  #x=apply(x, 2, function(x) scale(x))
  
  #---real data bootstrap lasso unweighted---
  set.seed(i)
  fit<-try(TRB_LASSO(x=x,y=y,intercept = T,standardize = T,nfolds = 5,alpha=0.05/length(gene_list),weights=rep(1,length(y)),penalty.factor=penalty.factor,B=500)) ######
  
  ans$beta_median_mrjti = fit$beta_median[1]
  ans$beta_mean_mrjti = fit$beta_mean[1]
  ans$p_mrjti = fit$pvalue
  
  }else{
    
    #not enough IVs, reduce to ivw
    ans$beta_median_mrjti = res@Estimate
    ans$beta_mean_mrjti = res@Estimate
    ans$p_mrjti = res@Pvalue
  }
  
  print(i)
  
  unlist(ans)

}

output = as.data.frame(output)
output$gene_id = gene_list
output$alpha = alpha

stopImplicitCluster()

write.table(output,paste0('/data/g*/z*/covid19_neanderthals/mrjti/simu/raw/alpha_',alpha,'.txt'),quote = F,sep='\t',row.names = F)
















