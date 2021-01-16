
#twas

#ml GCC/5.4.0-2.26  OpenMPI/1.10.3 pandas/0.18.1-Python-2.7.12 numpy/1.11.1-Python-2.7.12 scipy/0.17.0-Python-2.7.12 R

args=as.numeric(commandArgs(TRUE))

#info<-read.table('/data/g***/z***/data/ukbb/bn_blood_cell/info.txt',header = T,stringsAsFactors = F,sep='\t')

pheno='covid19'

tissue_list<-dir('/data/c***/z***/projects/gtex/weights/st/')

for (model in c('st','xt')){ #st=PrediXcan xt=JTI

tissue=tissue_list[args]

cmd=paste0('python2.7 /data/c***/z***/projects/cross_tissue/metaxcan/MetaXcan/software/MetaXcan.py   --model_db_path /data/c***/z***/projects/gtex/db/combined/',model,'_',tissue,'.db   --covariance /data/c***/z***/projects/gtex/cov/combined/',model,'_',tissue,'.txt.gz   --gwas_file /data/g***/z***/covid19_neanderthals/covid19_gwas/raw/GCST90000256.rsid.txt   --snp_column rsid   --effect_allele_column effect_allele   --non_effect_allele_column other_allele   --beta_column beta   --pvalue_column p_value  --output_file /data/g***/z***/covid19_neanderthals/asso/ukbb/raw/nejm/',pheno,'_',tissue,'_',model,'.csv')
system(cmd,wait = T)

}
