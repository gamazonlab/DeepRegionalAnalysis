
#local h2 estimation

library('MASS')
library('Matrix')

#segment position
start_pos=45859651
stop_pos=45909024

#sample size (ascertained)
n_cases=7885
n_controls=961804

#load MR results
dataset='covid19hg_112020'
mr<-read.csv(paste0('/data/g***/z***/covid19_neanderthals/mr/result/',dataset,'_r2_0.1.csv'),header = T,stringsAsFactors = F)
mr<-mr[mr$mrjti_sig=='sig',]
mr<-mr[!is.na(mr$mrjti_beta_mean),]

#gene annotation
anno<-read.table('/data/c***/z***/anno/gencode/37/gencode.v32.GRCh37.txt',header = T,stringsAsFactors = F)
mr<-merge(mr,anno,by='geneid')

#load variants
variants<-read.table(paste0('/data/g***/z***/covid19_neanderthals/expression_mediated_h2/pooled_gwas.clumped'),header = T,stringsAsFactors = F)

#EA samples
ea_sample<-read.table('/data/g***/z***/gtex/info/v8_sample_ea.fam',header = F,stringsAsFactors = F)
ea_sample$sampleid=sub('GTEX-','',ea_sample$V1)

#sample size
ascertained_prevalence=n_cases/(n_cases+n_controls)
N_effective = 4*ascertained_prevalence*(1-ascertained_prevalence)*(n_cases+n_controls)

i=1 #for test

for (i in 1:nrow(mr)){
  
  print(i)
  
  geneid=mr[i,'geneid']
  tissue=mr[i,'tissue']
  alpha=mr[i,'mrjti_beta_mean']
  alpha_se=mr[i,'mrjti_beta_se']
  
  #load snp info
  snp_info<-readRDS(paste0('/data/c***/z***/data/gtex/geno/v8/gene/info_1m/',geneid,'.rds'))
  snp_info$pos=sapply(snp_info$chr_bp,function(x) strsplit(x,"[_]")[[1]][2])
  snp_info$in_segment<-ifelse((snp_info$pos>start_pos & snp_info$pos<stop_pos),1,0)
  snp_info<-snp_info[which(snp_info$rsid %in% variants$SNP),]
  in_segment_pos<-which(snp_info$in_segment==1)
  
  #load genotype in dosage
  dosage<-readRDS(paste0('/data/c***/z***/data/gtex/geno/v8/gene/dosage_1m/',geneid,'.rds'))
  dosage<-dosage[which(dosage[,1] %in% ea_sample$sampleid),] #exclude non EA samples
  dosage[,-1]<-sapply(dosage[,-1], function(x) scale(x)) #scale genotype in dosage
  dosage<-dosage[,c('sampleid',snp_info$rsid)]
  
  #load expression file (residual)
  exp<-readRDS(paste0('/data/c***/z***/data/gtex/exp/v8/weighted/',tissue,'/',geneid,'.rds'))
  exp<-exp[exp$tissue==tissue,]
  exp$sampleid=sub('GTEX.','',exp$sampleid)
  
  #estimate beta_eQTL
  eqtl_df<-merge(exp[,c('sampleid','exp')],dosage,by='sampleid')
  for (j in 1:nrow(snp_info)){
    rsid=snp_info[j,'rsid']
    ans<-summary(lm(eqtl_df[,'exp']~eqtl_df[,rsid]))
    snp_info[j,'beta_eqtl']=ans$coefficients['eqtl_df[, rsid]','Estimate']
    snp_info[j,'se_eqtl']=ans$coefficients['eqtl_df[, rsid]','Std. Error']
  }
  
  
  #---full---
  
  full<-snp_info
  
  #estimate snp-snp var-cov matrix
  full_c<-cov(as.matrix(dosage[,-1]))
  
  #Moore-Penrose pseudoinverse
  inverse_full_c<-ginv(full_c)
  
  #q=rank(C)  the maximum number of linearly independent vectors in rows or columns
  q=rankMatrix(inverse_full_c)
  q=as.numeric(q)
  
  #pseudo gwas beta
  alpha=mr[i,'mrjti_beta_mean']
  full$pseudo_gwas_beta = alpha * full$beta_eqtl
  
  #pseudo gwas se
  alpha_se=mr[i,'mrjti_beta_se']
  full$pseudo_gwas_se = ((alpha_se)^2 * (full$se_eqtl)^2 +
                               (alpha_se)^2 * (full$beta_eqtl)^2 +
                               (alpha)^2 * (full$se_eqtl)^2)^0.5
  
  #n
  mr[i,'full_n'] = n = N_effective
  #p
  mr[i,'full_p'] = p = nrow(full)
  #h2
  mr[i,'full_h2_local_unadj'] = (t(full$pseudo_gwas_beta) %*% inverse_full_c %*% full$pseudo_gwas_beta - q/n) / (n-q) *n


  #---segment---
  
  segment<-full[full$in_segment==1,]
  
  #estimate snp-snp var-cov matrix
  segment_c<-full_c[which(full$in_segment==1),which(full$in_segment==1)]
  
  #Moore-Penrose pseudoinverse
  inverse_segment_c<-ginv(segment_c)
  
  #q=rank(C)  the maximum number of linearly independent vectors in rows or columns
  q=rankMatrix(inverse_segment_c)
  q=as.numeric(q)
  
  #pseudo gwas beta
  alpha=mr[i,'mrjti_beta_mean']
  segment$pseudo_gwas_beta = alpha * segment$beta_eqtl
  
  #pseudo gwas se
  alpha_se=mr[i,'mrjti_beta_se']
  segment$pseudo_gwas_se = ((alpha_se)^2 * (segment$se_eqtl)^2 +
                               (alpha_se)^2 * (segment$beta_eqtl)^2 +
                               (alpha)^2 * (segment$se_eqtl)^2)^0.5
  
  #n
  mr[i,'segment_n']=n=N_effective
  #p
  mr[i,'segment_p']=p=nrow(segment)
  #h2
  mr[i,'segment_h2_local_unadj']=(t(segment$pseudo_gwas_beta) %*% inverse_segment_c %*% segment$pseudo_gwas_beta - q/n) / (n-q) *n

  
  
}

mr$pi_c=mr$segment_h2_local_unadj/mr$full_h2_local_unadj

mr$concentration_ratio=mr$segment_h2_local_unadj*mr$full_p/mr$segment_p/mr$full_h2_local_unadj

write.table(mr,paste0('/data/g***/z***/covid19_neanderthals/hess/result/h2_segment_',dataset,'.txt'),quote = F,sep='\t',row.names = F)




