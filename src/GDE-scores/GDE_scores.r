
pheno='covid19hg_112020'

df<-read.csv(paste0('/data/g***/z***/covid19_neanderthals/mr/result/',pheno,'_r2_0.1.csv'),header = T,stringsAsFactors = F)

df$egger_sig<-ifelse(p.adjust(df$egger_p)<0.05,1,0)
df$median_weighted_sig<-ifelse(p.adjust(df$median_weighted_p)<0.05,1,0)
df$mrjti_sig<-ifelse(p.adjust(df$mrjti_pvalue)<0.05,1,0)

df$sig=df$egger_sig+df$median_weighted_sig+df$mrjti_sig

#distribution of GReX across populations

library(data.table)
library(ggplot2)

model='xt'

gene_tissue=df[df$sig>=1,]

#load weights
for (i in 1:nrow(gene_tissue)){
  tissue=gene_tissue[i,'tissue']
  geneid=gene_tissue[i,'geneid']
  
  tmp=readRDS(paste0('/data/c***/z***/projects/gtex/weights/',model,'/',tissue,'/',geneid,'.rds'))
  if(i==1){
    weights_all<-tmp
  }else{
    weights_all<-rbind(weights_all,tmp)
  }
}

write.table(weights_all[,c('rsid','gene')],paste0('/data/g***/z***/covid19_neanderthals/1kgp/pred_exp/info/snp_list_',pheno,'_asso'),sep='\t',quote = F,row.names = F)














#------------predict expression in Neanderthal------------

#--Neanderthal pre-process--
#extract
cmd='plink --bfile /data/g***/z***/covid19_neanderthals/neanderthals/AltaNea --chr 3 --from-bp 42859651 --to-bp 48909024 --make-bed --out /data/g***/z***/covid19_neanderthals/neanderthals/seg/seg_3mb'
system(cmd,wait = T)
#map to rsid
bim<-read.table('/data/g***/z***/covid19_neanderthals/neanderthals/seg/seg_3mb.bim',header = F,stringsAsFactors = F)
bim$id=seq(1,nrow(bim))
bim$chr_bp=paste0(bim$V1,'_',bim$V4)

bim<-merge(bim,weights_all[!duplicated(weights_all$rsid),c('rsid','chr_bp','ref_allele','counted_allele')],by='chr_bp',all.x=T)

bim<-bim[order(bim$id),]
write.table(bim[,c('V1','rsid','V3','V4','V5','V6')],'/data/g***/z***/covid19_neanderthals/neanderthals/seg/seg_3mb.bim',quote = F,sep='\t',row.names = F,col.names = F)
cmd=paste0('plink --bfile /data/g***/z***/covid19_neanderthals/neanderthals/seg/seg_3mb --extract /data/g***/z***/covid19_neanderthals/1kgp/pred_exp/info/snp_list_',pheno,'_asso --make-bed --out /data/g***/z***/covid19_neanderthals/neanderthals/seg/seg_3mb_rsid')
system(cmd,wait = T)
cmd='plink --bfile /data/g***/z***/covid19_neanderthals/neanderthals/seg/seg_3mb_rsid --recode A --out /data/g***/z***/covid19_neanderthals/neanderthals/seg/seg_3mb_rsid'
system(cmd,wait = T)

#load nea dosage
nea_dosage_all<-read.table('/data/g***/z***/covid19_neanderthals/neanderthals/seg/seg_3mb_rsid.raw',header = T,stringsAsFactors = F)
nea_dosage_all<-nea_dosage_all[,-seq(2,6)]
nea_allele<-read.table('/data/g***/z***/covid19_neanderthals/neanderthals/seg/seg_3mb_rsid.bim',stringsAsFactors = F)
nea_allele<-nea_allele[,c('V2','V6')]
colnames(nea_allele)<-c('rsid','nea_allele')
colnames(nea_dosage_all)[-1]<-sapply(colnames(nea_dosage_all)[-1], function(x) strsplit(x,"[_]")[[1]][1])

#--kgp pre-process--
#extract snps
cmd='plink --bfile /data/g***/z***/covid19_neanderthals/1kgp/geno/chr3_seg_5mb_maf0.001 --extract /data/g***/z***/covid19_neanderthals/1kgp/pred_exp/info/snp_list --recode A --out /data/g***/z***/covid19_neanderthals/1kgp/pred_exp/dosage'
system(cmd,wait = T)

#load dosage
kgp_dosage_all<-read.table('/data/g***/z***/covid19_neanderthals/1kgp/pred_exp/dosage.raw',header = T,stringsAsFactors = F)
kgp_dosage_all<-kgp_dosage_all[,-seq(2,6)]
kgp_allele<-data.frame(rsid=sapply(colnames(kgp_dosage_all)[-1], function(x) strsplit(x,"[_]")[[1]][1]),kgp_allele=sapply(colnames(kgp_dosage_all)[-1], function(x) strsplit(x,"[_]")[[1]][2]),stringsAsFactors = F)
colnames(kgp_dosage_all)[-1]<-kgp_allele$rsid

#load kgp available snp
kgp_available_allele<-read.table('/data/g***/z***/covid19_neanderthals/1kgp/geno/chr3_seg_5mb_maf0.001.bim',stringsAsFactors = F)

#sample info annotation
sample_anno<-read.table('/data/g***/z***/covid19_neanderthals/1kgp/info/igsr_samples.tsv',header = T,stringsAsFactors = F,sep='\t')
sample_anno<-sample_anno[,c('Sample.name','Population.code','Superpopulation.code')]


dir.create(paste0('/data/g***/z***/covid19_neanderthals/neanderthals/pred_exp/tissue_gene/',pheno))

i=1
#--pred exp--
for (i in 1:nrow(gene_tissue)){
  print(i)
  
  tissue=gene_tissue[i,'tissue']
  geneid=gene_tissue[i,'geneid']
  genename=gene_tissue[i,'gene_name']
  
  #load weight
  #weights=readRDS(paste0('/data/g***/z***/covid19_neanderthals/gtex/weights/st/',tissue,'/',geneid,'.rds'))
  weights=readRDS(paste0('/data/c***/z***/projects/gtex/weights/',model,'/',tissue,'/',geneid,'.rds'))
  
  
  #rm kgp unavailable snps
  weights<-weights[which(weights$rsid %in% kgp_available_allele$V2),]
  
  #rm nea unavailable snps
  nea_allele_tmp<-nea_allele[which(nea_allele$rsid %in% weights$rsid),]
  nea_allele_tmp<-nea_allele_tmp[!duplicated(nea_allele_tmp$rsid),]
  
  if(nrow(nea_allele_tmp)<=1){next}
  weights<-weights[sapply(nea_allele_tmp$rsid, function(x) which(weights$rsid==x)),]
  
  #harmonize nea allele
  weights$nea_allele=nea_allele_tmp$nea_allele
  #remove a/t c/g
  weights<-weights[which(!(paste0(weights$ref_allele,weights$counted_allele) %in% c('AT','TA','CG','GC'))),]
  weights$nea_allele_flip=sapply(weights$nea_allele,function(x)     switch (
    x,
    "A" = "T",
    "C" = "G",
    "T" = "A",
    "G" = "C",
    return(NA)
  ))
  weights$nea_counted_allele<-ifelse(weights$counted_allele==weights$nea_allele,weights$nea_allele,weights$nea_allele_flip)
  
  #effect nea
  weights$effect_nea=ifelse(weights$counted_allele==weights$nea_counted_allele,weights$weight,weights$weight*-1)
  #pred exp nea
  nea_exp=sum(2*weights$effect_nea)
  nea_pred_exp<-data.frame(FID='Neanderthal',pred_exp=nea_exp,stringsAsFactors = F)
  colnames(nea_pred_exp)[2]<-paste0(tissue,':',genename)
  
  
  #kgp
  if(length(which(kgp_allele$rsid %in% weights$rsid))<2){next}
  kgp_allele_tmp<-kgp_allele[which(kgp_allele$rsid %in% weights$rsid),]
  
  kgp_dosage_tmp<-kgp_dosage_all[,c('FID',kgp_allele_tmp$rsid)]
  
  weights<-weights[sapply(kgp_allele_tmp$rsid, function(x) which(weights$rsid==x)),]
  weights$effect_kgp<-ifelse(weights$counted_allele==kgp_allele_tmp$kgp_allele,weights$weight,weights$weight*-1)
  kgp_dosage_tmp[,paste0(tissue,':',genename)]<-as.matrix(kgp_dosage_tmp[,-1]) %*% weights$effect_kgp
  kgp_dosage_tmp<-kgp_dosage_tmp[,c('FID',paste0(tissue,':',genename))]
  
  #sample annotation
  pred_exp<-merge(kgp_dosage_tmp,sample_anno,by=1)
  
  #bw
  bw=0.08
  #if(length(unique(pred_exp[,2]))>3){bw=0.05}
  #if(length(unique(pred_exp[,2]))>5){bw=0.1}
  
  #plot
  df<-pred_exp
  df$population=df$Superpopulation.code
  colnames(df)[2]<-'Predicted_expression'
  png(paste0('/data/g***/z***/covid19_neanderthals/neanderthals/pred_exp/tissue_gene/',pheno,'/',genename,':',tissue,'.png'),width = 1000,height = 800,res=225)
  p<-ggplot(df, aes(x=Predicted_expression, color=population, fill=population)) +
    geom_histogram(aes(y=..density..), position="identity", alpha=0,bins = 15)+
    xlim(min(df[,2],nea_exp)-0.15,max(df[,2],nea_exp)+0.15)+
    geom_density(alpha=0.1,bw=bw)+
    theme_classic()+
    geom_vline(xintercept = nea_exp, linetype="dotted", 
               color = "red", size=1.5)+
    xlab(paste0('predicted expression of ',genename,' in ',tissue))+
    theme(legend.position="top",
          axis.title.x = element_text(size=9))
  print(p)
  dev.off()
  
  
  
  #merge kgp nea and output
  pred_exp<-rbind(kgp_dosage_tmp,nea_pred_exp)
  
  if(i==1){
    pred_exp_output=pred_exp
  }else{
    pred_exp_output=merge(pred_exp_output,pred_exp)
  }
  
}

#write.table(pred_exp_output,'/data/g***/z***/covid19_neanderthals/neanderthals/pred_exp/predixcan_egger_sig.txt',quote = F,sep='\t',row.names = F)
















