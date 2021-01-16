
#twas

#ml GCC/5.4.0-2.26  OpenMPI/1.10.3 pandas/0.18.1-Python-2.7.12 numpy/1.11.1-Python-2.7.12 scipy/0.17.0-Python-2.7.12 R

args=as.numeric(commandArgs(TRUE))

info<-read.table('/data/g***/z***/data/ukbb/bn_blood_cell/info.txt',header = T,stringsAsFactors = F,sep='\t')

pheno=info[args,'Phenotype_Code']

tissue_list<-dir('/data/c***/z***/projects/gtex/weights/st/')

model='xt'

for (i in 1:length(tissue_list)){
  tissue=tissue_list[i]
  
  cmd=paste0('python2.7 /data/c***/z***/projects/cross_tissue/metaxcan/MetaXcan/software/MetaXcan.py   --model_db_path /data/c***/z***/projects/gtex/db/combined/',model,'_',tissue,'.db   --covariance /data/c***/z***/projects/gtex/cov/combined/',model,'_',tissue,'.txt.gz   --gwas_file /data/g***/z***/data/ukbb/bn_blood_cell/processed/',pheno,'_chr3.txt   --snp_column rsid   --effect_allele_column eff_allele   --non_effect_allele_column ref_allele   --beta_column beta   --pvalue_column pval   --output_file /data/g***/z***/covid19_neanderthals/asso/ukbb/raw/',pheno,'_',tissue,'_',model,'.csv')
  system(cmd,wait = T)
  
  tmp<-read.csv(paste0('/data/g***/z***/covid19_neanderthals/asso/ukbb/raw/',pheno,'_',tissue,'_',model,'.csv'),header = T,stringsAsFactors = F)
  tmp$tissue=tissue
  tmp$phenoid=pheno
  if(i==1){
    df<-tmp
  }else{
    df<-rbind(df,tmp)
  }
  
}


write.table(df,paste0('/data/g***/z***/covid19_neanderthals/asso/ukbb/trait/xt_',pheno,'.csv'),quote =F,sep=',',row.names = F)

