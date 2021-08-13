
#local h2 estimation
library('RiskPortfolios')
library('MASS')
library('Matrix')
.libPaths("~/R/rlib-3.4_new")
library('pdist')
library('corpcor')


#segment position
start_pos=45859651
stop_pos=45909024

#sample size (ascertained)
n_cases=1610
n_controls=2205

#load esimated effect size from MR results
dataset='covid_nejm'
df<-read.csv(paste0('/data/g*/z*/covid19_neanderthals/mr/result/',dataset,'_r2_0.1.csv'),header = T,stringsAsFactors = F)
df<-df[df$mrjti_sig=='sig',]
df<-df[!is.na(df$mrjti_beta_mean),]

#gene annotation
anno<-read.table('/data/coxvgi/z*/anno/gencode/37/gencode.v32.GRCh37.txt',header = T,stringsAsFactors = F)
df<-merge(df,anno,by='geneid')

#sample size
ascertained_prevalence=n_cases/(n_cases+n_controls)
N_effective = 4*ascertained_prevalence*(1-ascertained_prevalence)*(n_cases+n_controls)

i=1 #for test

for (i in 1:nrow(df)){
  
  print(i)
  
  geneid=df[i,'geneid']
  tissue=df[i,'tissue']
  alpha=df[i,'mrjti_beta_mean']
  alpha_se=df[i,'mrjti_beta_se']
  
  #get chr pos
  chr=as.numeric(sub('chr','',anno[which(anno$geneid == geneid),'chr']))
  left=max(anno[which(anno$geneid == geneid),'left']-1e6,1)
  right=anno[which(anno$geneid == geneid),'right']+1e6
  
  #extract snp using plink
  cmd<-paste0('plink --bfile /data/g*/z*/covid19_neanderthals/expression_mediated_h2/geno/geno --chr ',chr,' --from-bp ',left,' --to-bp ',right,' --recode A --out /data/g*/z*/covid19_neanderthals/expression_mediated_h2/geno/tmp/',geneid,'_1m; plink --bfile /data/g*/z*/covid19_neanderthals/expression_mediated_h2/geno/geno --chr ',chr,' --from-bp ',left,' --to-bp ',right,' --make-bed --out /data/g*/z*/covid19_neanderthals/expression_mediated_h2/geno/tmp/',geneid,'_1m')
  system(cmd,wait = T)
  
  #load dosage for
  dosage_raw<-try(read.table(paste0('/data/g*/z*/covid19_neanderthals/expression_mediated_h2/geno/tmp/',geneid,'_1m.raw'),header = T,stringsAsFactors = F))
  if('try-error' %in% class(dosage_raw)){next}
  dosage<-dosage_raw[,-c(1,3:6)]  #rm useless cols
  
  #load info
  info_raw<-read.table(paste0('/data/g*/z*/covid19_neanderthals/expression_mediated_h2/geno/tmp/',geneid,'_1m.bim'),stringsAsFactors = F)
  
  #rm raw file
  cmd<-paste0('rm /data/g*/z*/covid19_neanderthals/expression_mediated_h2/geno/tmp/',geneid,'_1m*'); system(cmd,wait = T)
  
  #info
  info_raw$counted_allele<-sapply(colnames(dosage)[-1],function(x) strsplit(x,"[_]")[[1]][2])
  info_raw$ref_allele<-ifelse(info_raw$counted_allele==info_raw$V5,info_raw$V6,info_raw$V5)
  info_raw$chr_bp<-paste0(info_raw$V1,'_',info_raw$V4)
  colnames(info_raw)[2]<-'rsid'
  snp_info<-info_raw[,c('rsid','chr_bp','ref_allele','counted_allele')]
  snp_info$pos=sapply(snp_info$chr_bp,function(x) strsplit(x,"[_]")[[1]][2])
  snp_info$in_segment<-ifelse((snp_info$pos>start_pos & snp_info$pos<stop_pos),1,0)
  in_segment_pos<-which(snp_info$in_segment==1)
  
  #dosage
  dosage$IID<-sub('^.....','',dosage$IID)
  colnames(dosage)<-c('sampleid',sapply(colnames(dosage)[-1],function(x) strsplit(x,"[_]")[[1]][1]))
  
  #post imputation imputation replace by mean
  if(ncol(dosage)>2){
    dosage[,-1]<-round(apply(dosage[,-1], 2, function(x) ifelse(is.na(x),mean(x,na.rm=T),x)),3)
  }
  
  #scale genotype in dosage
  dosage[,-1]<-sapply(dosage[,-1], function(x) scale(x)) 
  dosage<-dosage[,c('sampleid',snp_info$rsid)]
  
  #load expression file (residual)
  exp<-readRDS(paste0('/data/coxvgi/z*/data/gtex/exp/v8/weighted/',tissue,'/',geneid,'.rds'))
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
  full_c <- covEstimation(as.matrix(dosage[,-1]), control = list(type = 'cor'))

  #Moore-Penrose pseudoinverse
  inverse_full_c<-ginv(full_c)
  
  #q=rank(C)  the maximum number of linearly independent vectors in rows or columns
  q=rankMatrix(inverse_full_c)
  q=as.numeric(q)
  
  #pseudo gwas beta
  alpha=df[i,'mrjti_beta_mean']
  full$pseudo_gwas_beta = alpha * full$beta_eqtl
  
  #pseudo gwas se
  alpha_se=df[i,'mrjti_beta_se']
  full$pseudo_gwas_se = ((alpha_se)^2 * (full$se_eqtl)^2 +
                           (alpha_se)^2 * (full$beta_eqtl)^2 +
                           (alpha)^2 * (full$se_eqtl)^2)^0.5
  
  #n
  df[i,'full_n'] = n = N_effective
  #p
  df[i,'full_p'] = p = nrow(full)
  #h2
  df[i,'full_h2_local'] = (t(full$pseudo_gwas_beta) %*% inverse_full_c %*% full$pseudo_gwas_beta - q/n) / (n-q) *n
  
  
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
  alpha=df[i,'mrjti_beta_mean']
  segment$pseudo_gwas_beta = alpha * segment$beta_eqtl
  
  #pseudo gwas se
  alpha_se=df[i,'mrjti_beta_se']
  segment$pseudo_gwas_se = ((alpha_se)^2 * (segment$se_eqtl)^2 +
                              (alpha_se)^2 * (segment$beta_eqtl)^2 +
                              (alpha)^2 * (segment$se_eqtl)^2)^0.5
  
  #n
  df[i,'segment_n']=n=N_effective
  #p
  df[i,'segment_p']=p=nrow(segment)
  #h2
  df[i,'segment_h2_local']=(t(segment$pseudo_gwas_beta) %*% inverse_segment_c %*% segment$pseudo_gwas_beta - q/n) / (n-q) *n
  
  
}

df$pi_c=df$segment_h2_local/df$full_h2_local

df$concentration_ratio=df$segment_h2_local*df$full_p/df$segment_p/df$full_h2_local

write.table(df,paste0('/data/g*/z*/covid19_neanderthals/expression_mediated_h2/result/h2_segment_',dataset,'.txt'),quote = F,sep='\t',row.names = F)







