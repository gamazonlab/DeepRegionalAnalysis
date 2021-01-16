#ml GCC/5.4.0-2.26  OpenMPI/1.10.3 pandas/0.18.1-Python-2.7.12 numpy/1.11.1-Python-2.7.12 scipy/0.17.0-Python-2.7.12 R
#Biovu complications

args=as.numeric(commandArgs(TRUE))  #421

#extract segment region
#cmd='plink --bfile /data/g***/z***/biovu/90k/geno/gtex/ea --chr 3 --from-bp 42859651 --to-bp 48909024 --make-bed --out /data/g***/z***/covid19_neanderthals/biovu/geno/seg'
#system(cmd,wait = T)

#info
biovu<-read.table('/data/c***/z***/data/biovu/pheno/k90/from_tyne/sample_size_info.txt',header = T,stringsAsFactors = F)  #1766 jobs
biovu<-biovu[biovu$N_cases>=100,]
biovu$phecode_num<-as.numeric(sub('X','',biovu$phecode))

pheinfo<-read.csv('/data/c***/z***/anno/phecode/phecode_definitions1.2.csv',header = T,stringsAsFactors = F)

pheinfo<-pheinfo[which(pheinfo$category %in% c("neurological", "respiratory", "circulatory system","endocrine/metabolic")),]
colnames(pheinfo)[1]<-'phecode_num'
pheinfo<-merge(biovu,pheinfo,by='phecode_num')


#gwas
phecode=pheinfo[args,'phecode']
cmd=paste0('plink2 --bfile /data/g***/z***/covid19_neanderthals/biovu/geno/seg  --allow-no-sex  --logistic hide-covar --pheno /data/c***/z***/data/biovu/pheno/k90/from_tyne/Phetable_MEGA_Acute_CEU_Covs.txt --pheno-name ',phecode,' --covar /data/c***/z***/data/biovu/pheno/k90/from_tyne/Phetable_MEGA_Acute_CEU_Covs.txt --covar-name is.male,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,Genotype.Batch --covar-variance-standardize --out /data/g***/z***/covid19_neanderthals/biovu/gwas/raw/seg')
#system(cmd,wait = T)
#remove #
#gwas<-read.table(paste0('/data/g***/z***/covid19_neanderthals/biovu/gwas/raw/seg.',phecode,'.glm.logistic'),header = T,stringsAsFactors = F,comment.char = '&')
#write.table(gwas,paste0('/data/g***/z***/covid19_neanderthals/biovu/gwas/raw/seg.',phecode,'.glm.logistic'),quote = F,sep='\t',row.names = F)

#twas
tissue_list<-dir('/data/c***/z***/projects/gtex/weights/st/')
model='xt'
for (i in 1:length(tissue_list)){
  tissue=tissue_list[i]
  
  cmd=paste0('python2.7 /data/c***/z***/projects/cross_tissue/metaxcan/MetaXcan/software/MetaXcan.py   --model_db_path /data/c***/z***/projects/gtex/db/combined/',model,'_',tissue,'.db   --covariance /data/c***/z***/projects/gtex/cov/combined/',model,'_',tissue,'.txt.gz   --gwas_file /data/g***/z***/covid19_neanderthals/biovu/gwas/raw/seg.',phecode,'.glm.logistic   --snp_column ID   --effect_allele_column ALT   --non_effect_allele_column REF   --or_column OR   --pvalue_column P   --output_file /data/g***/z***/covid19_neanderthals/biovu/asso/raw/xt_',phecode,'_',tissue,'.csv')
  system(cmd,wait = T)
  
  tmp<-read.csv(paste0('/data/g***/z***/covid19_neanderthals/biovu/asso/raw/xt_',phecode,'_',tissue,'.csv'),header = T,stringsAsFactors = F)
  tmp$tissue=tissue
  tmp$phecode=phecode
  
  if(i==1){
    result=tmp
  }else{
    result=rbind(result,tmp)
  }

}

#summarize twas
write.csv(result,paste0('/data/g***/z***/covid19_neanderthals/biovu/asso/trait/xt_',phecode,'.csv'),quote = F,row.names = F)
















