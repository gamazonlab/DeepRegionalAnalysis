
#cross tissue XT-SCAN
#gtex v8

args<-as.numeric(commandArgs(TRUE))
library('glmnet')

#set up subjobs
run_id<-1
run_list<-list()
#for (i in c(13,31,49)){ #13,31,49  brain,liver,blood
for (i in c(32,49)){  #tissues  32 lung
#for (i in c(seq(1,49))[-c(13,31,49)]){ #1:49
  for (j in 1:50){  #subjobs per tissue
    run_list[[run_id]]<-c(i,j)
    run_id=run_id+1
  }
}

run_i<-run_list[[args]]
print(run_i)

folder='xt'

#tissue
tissue_list<-dir('/data/c*/z*/data/gtex/exp/v8/weighted')
tissue<-tissue_list[run_i[1]]

#mkdir
dir.create(paste0('/data/c*/z*/projects/gtex/weights/',folder,'/'))
dir.create(paste0('/data/c*/z*/projects/gtex/weights/',folder,'/',tissue))

#get gene list
#load annotation
anno<-read.table('/data/c*/z*/anno/gencode/37/gencode.v32.GRCh37.txt',header = T,stringsAsFactors = F)
anno<-anno[which(anno$chr=='chr3' & anno$left>(45859651-1e6) & anno$right<(45909024+1e6)),]

#load tissue level dhs similarity (when gene level dhs similarity is not available)
dhs_weight_tissue<-read.table(paste0('/data/c*/z*/projects/cross_tissue/encode/simularity_tissue/',tissue),header = T,stringsAsFactors = F)
#get gene level dhs similarity gene list
dhs_gene_list<-dir(paste0('/data/c*/z*/projects/cross_tissue/encode/simularity/',tissue,'/gene/'))

#gene list
exp_list<-sub('....$','',dir(paste0('/data/c*/z*/data/gtex/exp/v8/weighted/',tissue))) #expression 
geno_list<-sub('....$','',dir('/data/c*/z*/data/gtex/geno/v8/gene/dosage_1m/')) #dosage genotype
gene_list<-intersect(exp_list,geno_list)
gene_list<-intersect(anno$geneid,gene_list)


#---------------------------
#i=which(gene_list=='ENSG00000240654')

if(length(gene_list)>=run_i[2]){
  gene_id = gene_list[run_i[2]]
  
  print(paste0('INFO processing ',gene_id,' ...'))
  
  #load expression levels
  exp<-readRDS(paste0('/data/c*/z*/data/gtex/exp/v8/weighted/',tissue,'/',gene_id,'.rds'));colnames(exp)[3]<-'exp_w' 
  
  #merge with dhs weights
  dhs_weight_gene<-read.table(paste0('/data/g*/z*/covid19_neanderthals/jti/dhs/similarity/',tissue,'/table/covid19_segment_QN.txt'),header = T,stringsAsFactors = F)
  dhs_weight<-data.frame(tissue=colnames(dhs_weight_gene[-1]),dhs_w=as.numeric(dhs_weight_gene[1,-1]))
  dhs_weight[which(dhs_weight[,2]<0),2]<-0
  
  exp<-merge(exp,dhs_weight,by='tissue')  
  
  #-----data prepare-----
  #load geno snpinfo
  #1m
  geno_1m<-readRDS(paste0('/data/c*/z*/data/gtex/geno/v8/gene/dosage_1m/',gene_id,'.rds'))
  geno_1m[,1]<-paste0('GTEX.',geno_1m[,1])
  snp_info_1m<-readRDS(paste0('/data/c*/z*/data/gtex/geno/v8/gene/info_1m/',gene_id,'.rds'))
  snp_info_1m$pos = sapply(snp_info_1m$chr_bp, function(x) as.numeric(strsplit(x,"[_]")[[1]][2]))
  
  cis_left = anno[which(anno$geneid==gene_id),'left'] - 1e6
  cis_right = anno[which(anno$geneid==gene_id),'right'] + 1e6
  
  output = data.frame(seed = seq(1,100), r = NA, p = NA)
  
  for(k in 1:101){
    print(k)
    set.seed(k+i)
    start_pos = sample(cis_left:(cis_right-49373),1)
    end_pos = start_pos + 49373
    
    if(k == 1){
      start_pos = 45859651
      end_pos = 45909024
    }
    
    snp_info = snp_info_1m[which(snp_info_1m$pos > start_pos & snp_info_1m$pos < end_pos),]
    if(nrow(snp_info)<2){next}
    geno = geno_1m[,c('sampleid',snp_info$rsid)]
    
    #-----model training-----
    #fit single tissue model to get proper window size and lambda range for each gene
    exp_st<-exp[exp$tissue==tissue,] #load single tissue data
    d_st<-merge(exp_st,geno,by='sampleid')
    if(nrow(d_st)<10){next}
    
    set.seed(as.numeric(sub('^....','',gene_id)))
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
      #print(j)
      exp_power=r_matrix[j,1];dhs_power=r_matrix[j,2]
      exp$w<-(exp$exp_w)^exp_power*(exp$dhs_w)^dhs_power
      
      #merge
      d<-merge(exp,geno,by='sampleid')
      
      #rm tissues with low similarity levels (w<0.1)
      d<-d[d[,'w']>=0.1,]; d<-d[!is.na(d[,1]),]
      
      #id map
      d<-merge(id_map,d,by='sampleid')
      #print(dim(d))
      
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
    
    output[k, 'r'] = r_test
    output[k, 'p'] = p_test
    
  }
  
  
  output$tissue = tissue
  
  
  write.table(output,paste0('/data/g*/z*/covid19_neanderthals/random_segment/raw/',tissue,'_',gene_id,'.txt'),quote = F,sep='\t',row.names = F)
  
  
}









