
#JTI  
#gtex v8
args<-as.numeric(commandArgs(TRUE))

library('glmnet')

#set up subjobs
run_id<-1
run_list<-list()
#for (i in c(13,31,49)){ #13,31,49  brain,liver,blood
for (i in 1:49){  #tissues
  #for (i in c(seq(1,49))[-c(13,31,49)]){ #1:49
  for (j in 1:6){  #subjobs per tissue
    for (k in 1:3){
      run_list[[run_id]]<-c(i,j,k)
      run_id=run_id+1
    }
  }
}

run_i<-run_list[[args[1]]]
print(run_i)

region = switch(run_i[3],
                '1'='seg_only',
                '2'='seg_100k',
                '3'='seg_500k')
folder = region

#tissue
tissue_list<-dir('/data/c*/z*/data/gtex/exp/v8/weighted')
tissue<-tissue_list[run_i[1]]

#mkdir
dir.create(paste0('/data/g*/z*/covid19_neanderthals/gtex/weights/',folder,'/'))
dir.create(paste0('/data/g*/z*/covid19_neanderthals/gtex/weights/',folder,'/',tissue))

#get gene list
#load annotation
anno<-read.table('/data/c*/z*/anno/gencode/37/gencode.v32.GRCh37.txt',header = T,stringsAsFactors = F)
anno<-anno[which(anno$chr=='chr3' & anno$left>43859651 & anno$right<47909024),]



#generate genotype data
dosage_raw<-read.table(paste0('/data/g*/z*/covid19_neanderthals/gtex/geno/raw/',region,'.raw'),header = T,stringsAsFactors = F)
dosage_raw<-dosage_raw[,-seq(2,6)]
colnames(dosage_raw)[-1]<-sapply(colnames(dosage_raw)[-1], function(x) strsplit(x,"[_]")[[1]][1])
dosage<-dosage_raw
colnames(dosage)[1]<-'sampleid'
dosage$sampleid=sub('GTEX-','GTEX.',dosage$sampleid)

#post imp imp
dosage[,-1]<-sapply(dosage[,-1], function(x) ifelse(is.na(x),mean(x,na.rm = T),x))

#allele info
snp_info<-read.table(paste0('/data/g*/z*/covid19_neanderthals/gtex/geno/raw/',region,'.bim'),header = F,stringsAsFactors = F)
snp_info$chr_bp<-paste0(snp_info$V1,'_',snp_info$V4)
colnames(snp_info)<-c('chr','rsid','drop','pos','counted_allele','ref_allele','chr_bp')
snp_info_all<-snp_info[,c('rsid','chr_bp','ref_allele','counted_allele')]

#exp
exp_list<-sub('....$','',dir(paste0('/data/c*/z*/data/gtex/exp/v8/weighted/',tissue)))

gene_list<-intersect(exp_list,anno$geneid)

#sub job start and end id
i_start=(run_i[2]-1)*10+1
i_end=min(run_i[2]*10,length(gene_list))

#map info
map_info = read.table('/data/c*/z*/data/gtex/geno/v8/liftover/to_hg19_mapped.bed',header = F,stringsAsFactors = F)
map_info = map_info[map_info$V1=='chr3',]
colnames(map_info) = c('chr','pos1','pos2','gtexid')
map_info$chr_bp = paste0(sub('chr','',map_info$chr),'_',map_info$pos1)

#---------------------------
#i=which(gene_list=="ENSG00000173473")

if (i_start<length(gene_list)){
  for (i in i_start:i_end){
    print(paste0('INFO processing ',gene_list[i],' ...'))
    
    if(file.exists(paste0('/data/g*/z*/covid19_neanderthals/gtex/weights/',folder,'/',tissue,'/',gene_list[i],'.rds'))){next}
    
    #overlap with gene region
    snp_info_1m<-readRDS(paste0('/data/c*/z*/data/gtex/geno/v8/gene/info_1m/',gene_list[i],'.rds'))
    snp_info_1m = merge(snp_info_1m, map_info, by = 'chr_bp')
    mapped_pos = which(snp_info_all$rsid %in% snp_info_1m$gtexid)
    
    if(length(mapped_pos)==0){next}
    
    snp_info=snp_info_all[mapped_pos,]
    geno<-dosage[,c(1,mapped_pos)]
    
    #load expression levels
    exp<-readRDS(paste0('/data/c*/z*/data/gtex/exp/v8/weighted/',tissue,'/',gene_list[i],'.rds'));colnames(exp)[3]<-'exp_w' 
    
    #merge with dhs weights
    dhs_weight_gene<-read.table(paste0('/data/g*/z*/covid19_neanderthals/jti/dhs/similarity/',tissue,'/table/covid19_segment_QN.txt'),header = T,stringsAsFactors = F)
    dhs_weight<-data.frame(tissue=colnames(dhs_weight_gene[-1]),dhs_w=as.numeric(dhs_weight_gene[1,-1]))
    dhs_weight[which(dhs_weight[,2]<0),2]<-0
    
    exp<-merge(exp,dhs_weight,by='tissue')  
    #col: tissue sampleid exp exp_w dhs_w
    
    #-----model training-----
    #fit single tissue model to get proper window size and lambda range for each gene
    exp_st<-exp[exp$tissue==tissue,] #load single tissue data
    d_st<-merge(exp_st,geno,by='sampleid')
    
    if(nrow(d_st)<10){next}
    
    set.seed(as.numeric(sub('^....','',gene_list[i])))
    fit<-cv.glmnet(x=as.matrix(d_st[,6:ncol(d_st)]),y=as.matrix(d_st[,'exp']), nfolds = 5,keep = T,alpha=0.5,nlambda=50,pmax=200) #pmax=200
    
    #find proper window size and set up for downstream analysis
    #The cross tissue lambda range was set around the best single tissue lambda
    lambda_list<-fit$lambda[max(1,which(fit$lambda==fit$lambda.min)-20):min(length(fit$lambda),which(fit$lambda==fit$lambda.min)+20)]
    
    
    #dataframe for grid search and performance collection
    r_matrix<-as.data.frame(matrix(data=NA,nrow=(4*4),ncol=6))
    i_loop<-1
    for (exp_power in c(1,4,16,64)){
      for (dhs_power in c(1,4,16,64)){
        r_matrix[i_loop,1:2]<-c(exp_power,dhs_power)
        i_loop=i_loop+1
      }
    }
    colnames(r_matrix)<-c('exp_power','dhs_power','lambda','window','r_test','p_test')
    r_matrix$window='segment'
    
    #sample id map for each fold (fold=5)
    d_tmp<-merge(exp,geno,by='sampleid')
    sample_all<-unique(d_tmp$sampleid)
    set.seed(as.numeric(sub('^....','',gene_list[1]))+1)
    id_map<-data.frame(sampleid=sample_all,id_map=sample(rep(seq(1,5),ceiling(length(sample_all)/5)),length(sample_all)),stringsAsFactors = F)
    
    #beta list of each hyper-parameter pairs
    beta_list<-list()
    
    #go through each of the hyper-parameter pairs
    for (j in 1:nrow(r_matrix)){
      print(j)
      exp_power=r_matrix[j,1];dhs_power=r_matrix[j,2]
      exp$w<-(exp$exp_w)^exp_power*(exp$dhs_w)^dhs_power
      
      #merge
      d<-merge(exp,geno,by='sampleid')
      
      #rm tissues with low similarity levels (w<0.1)
      d<-d[d[,'w']>=0.1,]; d<-d[!is.na(d[,1]),]
      
      #id map
      d<-merge(id_map,d,by='sampleid')
      print(dim(d))
      
      #target tissue position (the performance will be only estimated in the target tissue)
      tt_pos<-which(d$tissue==tissue)
      
      #cross-tissue weighted elastic net
      set.seed(as.numeric(sub('^....','',gene_list[1]))+2)
      ans<-try(cv.glmnet(x=as.matrix(d[,8:ncol(d)]),y=as.matrix(d[,'exp']),weights=d[,'w'],foldid=d[,'id_map'],lambda=lambda_list,keep = T,pmax=200)) #pmax=200
      if ('try-error' %in% class(ans)){
        r_matrix[j,'r_test']<-0
      }else{
        #correlation between pred and obs expression levels
        cor_ans<-cor.test(ans$fit.preval[tt_pos,which(ans$lambda==ans$lambda.min)],d[tt_pos,'exp'])
        
        r_matrix[j,'r_test']<-cor_ans$estimate
        r_matrix[j,'p_test']<-cor_ans$p.value
        r_matrix[j,'lambda']<-ans$lambda.min
        
        #collect the weights
        beta_list[[j]]<-as.numeric(ans$glmnet.fit$beta[,which(ans$lambda==ans$lambda.min)])
        
      }
    }
    
    #find the best hyperparameters with the largest r_test
    best_row=which.max(r_matrix[,'r_test'])[1]
    r_test=r_matrix[best_row,'r_test']
    p_test=r_matrix[best_row,'p_test']
    
    #---output---
    if (r_test>0){  #p_test<0.05 & r_test>0.1
      snp_info$gene<-gene_list[i]
      snp_info$r2<-r_test^2
      snp_info$p<-p_test
      snp_info$lambda<-r_matrix[best_row,'lambda']
      snp_info$weight<-beta_list[[best_row]]
      snp_info<-snp_info[,c(5,1:4,9,6:8)]
      nocv<-snp_info[snp_info$weight!=0,]
      if (nrow(nocv)>0){
        out_path<-paste0('/data/g*/z*/covid19_neanderthals/gtex/weights/',folder,'/',tissue,'/',gene_list[i],'.rds') 
        saveRDS(nocv,file=out_path)
      }
    }
    
  }
}else{
  print('out of range')
}






