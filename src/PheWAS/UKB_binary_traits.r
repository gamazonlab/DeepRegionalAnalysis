#ml GCC/5.4.0-2.26  OpenMPI/1.10.3 pandas/0.18.1-Python-2.7.12 numpy/1.11.1-Python-2.7.12 scipy/0.17.0-Python-2.7.12 R


args=as.numeric(commandArgs(TRUE))  #254

phecode=dir('/data/c*/z*/data/ukbb/Lee/gwas_pruned/')
phecode=sub('.txt','',phecode)

tissue_list<-dir('/data/c*/z*/projects/gtex/weights/st/')


pheinfo<-read.csv('/data/c*/z*/anno/phecode/phecode_definitions1.2.csv',header = T,stringsAsFactors = F)

pheinfo<-pheinfo[which(pheinfo$category %in% c("neurological", "respiratory", "circulatory system","endocrine/metabolic")),]

phecode = intersect(phecode,pheinfo$phecode)

phecode=phecode[args]

#twas

model='xt'

for (i in 1:length(tissue_list)){
  tissue=tissue_list[i]
  
  cmd=paste0('python2.7 /data/c*/z*/projects/cross_tissue/metaxcan/MetaXcan/software/MetaXcan.py   --model_db_path /data/c*/z*/projects/gtex/db/combined/',model,'_',tissue,'.db   --covariance /data/c*/z*/projects/gtex/cov/combined/',model,'_',tissue,'.txt.gz   --gwas_file /data/c*/z*/data/ukbb/Lee/gwas_pruned/',phecode,'.txt   --snp_column ID   --effect_allele_column ALT   --non_effect_allele_column REF   --beta_column beta   --pvalue_column pval   --output_file /data/g_gamazon_lab/z*/covid19_neanderthals/asso/ukbb_phewas/raw/X',phecode,'_',tissue,'_',model,'.csv')
  system(cmd,wait = T)
  
  setwd('/data/c*/z*/data/ukbb/Lee/twas/')
  cmd=paste0('gzip /data/c*/z*/data/ukbb/Lee/twas/X',phecode,'_',tissue,'_',model,'.csv')
  system(cmd,wait = T)
}


