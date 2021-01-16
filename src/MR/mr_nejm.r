
args<-as.numeric(commandArgs(TRUE))

model='xt' #based on JTI model
options(digits = 4)
nf=5 #number of fold for cross-validation
r2_cutoff=0.1 #ld clumping

.libPaths("~/R/rlib-3.4_new")

library("glmnet")
#library("mvtnorm")
library('HDCI')
library('MRPRESSO')

library('MendelianRandomization')
#library('selectiveInference')
#library('hdi')

source('/home/zhoud2/script/gtex/github/MR-JTI.function')

#non_palindromic
#bim<-read.table('/data/c***/z***/data/gtex/geno/v8/liftover/pruning/gtex_rsid_prune0.9.bim',header = F,stringsAsFactors = F)
#bim<-bim[-which(paste0(bim$V5,bim$V6) %in% c('AT','TA','GC','CG')),]
#write.table(bim,'/data/c***/z***/data/gtex/geno/v8/liftover/pruning/gtex_rsid_prune0.9_non_palindromic.bim',quote = F,col.names = F,row.names = F,sep='\t')

#clumping
#cmd<-paste0('plink --bfile /data/c***/z***/data/gtex/geno/v8/liftover/pruning/gtex_rsid_prune0.9 --extract /data/c***/z***/data/gtex/geno/v8/liftover/pruning/gtex_rsid_prune0.9_non_palindromic.bim --clump /data/g***/z***/covid19_neanderthals/covid19_gwas/raw/GCST90000256.rsid.txt --clump-field p_value --clump-snp-field rsid --clump-p1 1 --clump-r2 0.1 --out /data/g***/z***/covid19_neanderthals/mrjti/clumping/nejm_0.1')
#system(cmd,wait = T)

clumped<-read.table(paste0('/data/g***/z***/covid19_neanderthals/mrjti/clumping/nejm_',r2_cutoff,'.clumped'),header = T,stringsAsFactors = F)

#load gwas result
gwas<-readRDS('/data/g***/z***/covid19_neanderthals/covid19_gwas/raw/GCST90000256.rsid.rds')
gwas<-gwas[,c('rsid','effect_allele','other_allele','beta','standard_error','p_value')]
colnames(gwas)<-c('rsid','effect_allele_gwas','ref_allele_gwas','beta_gwas','se_gwas','p_gwas')
gwas<-gwas[which(gwas$rsid %in% clumped$SNP),]  #all or pruned

#remove palindromic snps
gwas<-gwas[-which(paste0(gwas$effect_allele_gwas,gwas$ref_allele_gwas) %in% c('AT','TA','CG','GC')),]

#load TWAS results (only run for genes with significant TWAS results)
twas<-read.csv(paste0('/data/g***/z***/covid19_neanderthals/asso/covid19/xt_covid19.csv'),header = T,stringsAsFactors = F)
twas<-twas[p.adjust(twas$pvalue,method = 'BH')<0.05,] #only run MR for TWAS-significant genes
gene_list<-twas$gene

#load ldsc
ldsc<-read.table('/data/c***/z***/data/gtex/geno/v8/ldsc/gtex_pruned.score.ld',header = T,stringsAsFactors = F)
ldsc<-ldsc[,c(1,8)]

#built output df
output<-data.frame('geneid'=NA)

#subjob start and end id
i_start=(args-1)*10+1
i_end=args*10
if (i_end>length(gene_list)){
  i_end=length(gene_list)
}

i=1

if (i_start<length(gene_list)){
  for (i in i_start:i_end){
    print(i)
    output[i,'geneid']<-gene<-gene_list[i]
    output[i,'tissue']<-tissue<-twas[i,'tissue']
    output[i,'gene_name']<-gene_name<-twas[i,'gene_name']

    #load genotype in dosage and allele info
    gtex_geno_path<-'/data/c***/z***/data/gtex/geno/v8/gene/dosage_1m/'
    gtex_geno_info_path<-'/data/c***/z***/data/gtex/geno/v8/gene/info_1m/'
    gtex_geno<-readRDS(paste0(gtex_geno_path,gene,'.rds'))
    gtex_info<-readRDS(paste0(gtex_geno_info_path,gene,'.rds'))
    
    #overlap snps with gwas results
    if(length(which(colnames(gtex_geno) %in% gwas$rsid))==0){next}
    gtex_geno<-gtex_geno[,c(1,which(colnames(gtex_geno) %in% gwas$rsid))]
    gtex_info<-gtex_info[which(gtex_info$rsid %in% gwas$rsid),c(1,4)]
    
    #load expression
    gtex_exp_path<-'/data/c***/z***/data/gtex/exp/v8/weighted/'
    gtex_exp<-readRDS(paste0(gtex_exp_path,tissue,'/',gene,'.rds'))
    
    #run eqtl, get beta_eqtl
    gtex_geno[,1]<-paste0('GTEX.',gtex_geno[,1])
    gtex_exp<-gtex_exp[gtex_exp$tissue==tissue,c(1,2)]
    d<-merge(gtex_exp,gtex_geno,by=1)
    eqtl<-data.frame('rsid'=gtex_info$rsid,'eqtl_beta'=NA,'eqtl_se'=NA,'eqtl_p'=NA,stringsAsFactors = F)
    for (k in 1:nrow(eqtl)){
      eqtl_fit<-summary(lm(d$exp~d[,k+2]))
      eqtl[k,2]<-eqtl_fit$coefficients[2,1]
      eqtl[k,3]<-eqtl_fit$coefficients[2,2]
      eqtl[k,4]<-eqtl_fit$coefficients[2,4]
    }
    eqtl<-merge(eqtl,gtex_info,by='rsid')
    
    #merge eqtl results with gwas results
    df<-merge(eqtl,gwas,by='rsid')
    
    #match eqtl and gwas effect alleles
    df$gwas_beta<-ifelse(df$counted_allele==df$effect_allele_gwas,df$beta_gwas,(df$beta_gwas*-1))
    
    df_mr<-df[df$eqtl_p<0.05,]
    
    #SNP-SNP cor matrix, to remove remaining correlations
    df_cor<-d[,sapply(df_mr$rsid, function(x) which(colnames(d)==x))]
    cor_matrix<-cor(df_cor)
    cor_matrix<-matrix(data=as.numeric(cor_matrix),ncol = ncol(df_cor),nrow=ncol(df_cor))
    
    output[i,'n_eqtl']<-nrow(df_mr)
    
    if (nrow(df_mr)>=3){
      #----MR-EGGER---
      set.seed(as.numeric(sub('ENSG','',gene))+100) 
      ans_egger<-mr_egger(mr_input(bx = df_mr$eqtl_beta, bxse = df_mr$eqtl_se, by = df_mr$gwas_beta, byse = df_mr$se_gwas,corr=cor_matrix))
      output[i,'egger_beta']<-ans_egger@Estimate
      output[i,'egger_se']<-ans_egger@StdError.Est
      output[i,'egger_p']<-ans_egger@Pvalue.Est

      #----MR-Weighted_median----
      set.seed(as.numeric(sub('ENSG','',gene))+100) 
      ans_median_weighted<-mr_median(mr_input(bx = df_mr$eqtl_beta, bxse = df_mr$eqtl_se, by = df_mr$gwas_beta, byse = df_mr$se_gwas,corr=cor_matrix),weighting = "weighted", iterations = 100)
      output[i,'median_weighted_beta']<-ans_median_weighted@Estimate
      output[i,'median_weighted_se']<-ans_median_weighted@StdError
      output[i,'median_weighted_p']<-ans_median_weighted@Pvalue
    }

    
    #----MR-JTI----
    
    # #orientation (not necessary)
    # df$gwas_beta<-ifelse(df$eqtl_beta>0,df$gwas_beta,df$gwas_beta*-1)
    # df$eqtl_beta<-abs(df$eqtl_beta)

    #merge with ldsc
    df<-merge(df,ldsc,by=1)
    
    #skip genes with less than 20 obs (no way for cv)
    output[i,'n_ivs']<-nrow(df)
    if (nrow(df)<20){next}
    
    #no of sig gwas loci (only consider these for pleiotropy control)
    n_gwas<-length(which(df$p_gwas<0.05))
    
    #weighted by the precision
    #weights=1/(df$se_gwas)^2
    
    #penalty factor, no penalty for the first two features (eqtl_beta and ldsc)
    penalty.factor<-c(0,0,rep(1,n_gwas))
    
    #identity matrix
    df<-df[order(df$p_gwas),]
    df<-df[,c('rsid','gwas_beta','eqtl_beta','ldscore')] 
    if(n_gwas>0){
      df[,(ncol(df)+1):(ncol(df)+n_gwas)]<-diag(nrow(df))[,1:n_gwas]
    }
    
    #---assign x and y---
    y=df[,'gwas_beta'];x=as.matrix(df[,3:ncol(df)])
    #scale
    y=scale(y)
    x=apply(x, 2, function(x) scale(x))
    

    #---real data bootstrap lasso unweighted---
    set.seed(as.numeric(sub('ENSG','',gene))+100) 
    fit<-try(TRB_LASSO(x=x,y=y,intercept = T,standardize = T,nfolds = nf,alpha=0.05/length(gene_list),weights=rep(1,nrow(df)),penalty.factor=penalty.factor,B=1000)) ######
    if(!('try-error' %in% class(fit))){
      output[i,'mrjti_beta_median']<-fit$beta_median[1]
      #output[i,'beta_mean']<-fit$beta_mean[1]
      output[i,'mrjti_lower']<-fit$interval[1,1]
      output[i,'mrjti_upper']<-fit$interval[2,1]
      output[i,'mrjti_pvalue']<-fit$pvalue
    }
      
  }
  
  #sig or not according to 95%CI bonf adjusted
  output$mrjti_sig<-ifelse(output$mrjti_lower*output$mrjti_upper>0,'sig','nonsig')
  
  #rm NA lines
  output<-output[!is.na(output[,1]),]
  
  #output
  write.csv(output,paste0('/data/g***/z***/covid19_neanderthals/mr/result/covid_nejm_r2_',r2_cutoff,'_',args,'.csv'),quote = F,row.names = F)
  
}






















