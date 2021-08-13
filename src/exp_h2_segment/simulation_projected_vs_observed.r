

# #------VCF file (from 1000 gp) processing------

# extract a list of EUR sample list
# awk '{print $1}' /data/coxvgi/z*/data/kgp/phase3/chr/chr16_eur_nondup.fam > /data/g*/z*/covid19_neanderthals/ld_matrix_projected/genotype_simulation/1000gp/ceu.txt
#
# extract a list of common variants used
# awk '{print $2}' /data/g*/z*/covid19_neanderthals/ld_matrix_projected/geno/chr22.bim > /data/g*/z*/covid19_neanderthals/ld_matrix_projected/genotype_simulation/1000gp/variants.txt
# 
# 
# ml GCC OpenMPI VCFtools/0.1.15
# 
# vcftools --gzvcf /data/coxvgi/z*/data/kgp/phase3/chr/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
# --keep /data/g*/z*/covid19_neanderthals/ld_matrix_projected/genotype_simulation/1000gp/ceu.txt \
# --maf 0.01 \
# --chr 22 \
# --snps /data/g*/z*/covid19_neanderthals/ld_matrix_projected/genotype_simulation/1000gp/variants.txt \
# --from-bp 20000000 \
# --to-bp 25000000 \
# --recode \
# --out /data/g*/z*/covid19_neanderthals/ld_matrix_projected/genotype_simulation/1000gp/chr22



#args=as.numeric(commandArgs(TRUE)) 

library('sim1000G')
library('MASS')
library('Matrix')
library('RiskPortfolios')
.libPaths("~/R/rlib-3.4_new")
library('data.table')
library('pdist')
library('corpcor')
library('Matrix')

library(foreach)
library(doParallel)

segment_length = as.integer(50000)
sample_size = 500

i = 1 #FOR TEST

#set up parallele
registerDoParallel(cores = detectCores())

#generate random segments on chr 22
all_geno = foreach(i = seq(1,500,50)) %dopar% {
   
  #segment start and end position
  set.seed(i)
  start_point = as.integer(sample(seq(20000000,24000000,by=1000000),1,replace = F))
  end_point = as.integer(start_point + segment_length)
  
  #load vcf file
  vcf_file = '/data/g*/z*/covid19_neanderthals/ld_matrix_projected/genotype_simulation/1000gp/chr22.recode.vcf'
  vcf = readVCF(vcf_file, maxNumberOfVariants = 300 , 
                min_maf = 0.05, max_maf = 1,
                region_start = start_point, region_end = end_point)
  set.seed(i)
  
  #simulate genetype for 50k samples based on the empirical data (1000gp EUR)
  startSimulation(vcf, totalNumberOfIndividuals = 50000)
  ids = generateUnrelatedIndividuals(50000)
  
  genotype = retrieveGenotypes(ids)
  all_geno = as.data.frame(genotype)
  
  #scale to ~N(0,1)
  all_geno = apply(all_geno, MARGIN = 2, function(x) scale(x)) #scale to ~N(0,1)
  all_geno
}


#output datafrome
output = data.frame(seed = seq(1,500))


asso = function(x){
  ans = lm(outcome ~ x)
  as.numeric(ans$coefficients[2])
}

h2_est = function(cov_matrix, beta_gwas){
  #Moore-Penrose pseudoinverse
  inverse_c<-pseudoinverse(cov_matrix)
  #q=rank(C)  the maximum number of linearly independent vectors in rows or columns
  q=rankMatrix(inverse_c)
  q=as.numeric(q)
  h2 = (sample_size * t(beta_gwas) %*% inverse_c %*% beta_gwas - q) / (sample_size - q)
}


#500 times simulation
for(i in 1:500){
  print(i)
  
  #load all genotype
  set.seed(i)
  all = all_geno[[sample(length(all_geno),1)]]
  #take a random sample
  set.seed(i)
  sampled = all[sample(nrow(all),sample_size,replace = F),]
  
  #estimate cov matrix
  projected = covEstimation(as.matrix(sampled), control = list(type = 'cor')) 
  observed = cov(as.matrix(sampled))
  true = cov(as.matrix(all))
  
  output[i,'sample_size'] = sample_size
  output[i,'segment_length'] = segment_length
  output[i,'n_snps'] = ncol(all)
  output[i,'distance_true_observed'] <- sqrt(sum(as.numeric(true - observed))^2)
  output[i,'distance_true_projected'] <- sqrt(sum(as.numeric(true - projected))^2)
  
  for (h2 in c(0.01,0.02,0.03)){

    #define a sample set for beta (GWAS) estimation
    set.seed(i+2021)
    df_beta_est = all[sample(nrow(all),sample_size,replace = F),]
    #define a causal snp
    set.seed(i+100)
    causal_snp_pos = sample(ncol(df_beta_est),1)
    beta = sqrt(h2*1/var(df_beta_est[,causal_snp_pos])) #var(y) = 1
    #generate the outcome
    set.seed(i)
    outcome = beta*df_beta_est[,causal_snp_pos] + rnorm(nrow(df_beta_est),0,sqrt(1-h2))
    
    #estimate beta (GWAS)
    beta_gwas = apply(df_beta_est, MARGIN = 2, function(x) asso(x))
    
    #estimate h2 for the segment using beta_gwas and LD-matrix
    output[i,paste0('h2_',h2,'_projected')] = h2_est(projected, beta_gwas)
    output[i,paste0('h2_',h2,'_observed')] = h2_est(observed, beta_gwas)
    
  }
  
}

stopImplicitCluster()


write.table(output,paste0('/data/g*/z*/covid19_neanderthals/ld_matrix_projected/results/projected_observed_h2.txt'),quote = F,sep='\t',row.names = F)



#---generate the figure----

library(ggplot2)

#load the results (the same df a few lines above)
df = read.table(paste0('/data/g*/z*/covid19_neanderthals/ld_matrix_projected/results/projected_observed_h2.txt'),header = T,stringsAsFactors = F)


p = list()
i=1
for(i in 1:3){
  h2 = c(0.01,0.02,0.03)[i]
  tmp = data.frame(h2 = c(df[,paste0('h2_',h2,'_projected')],
                          df[,paste0('h2_',h2,'_observed')]),
                   group = c(rep('Projected',nrow(df)), rep('Observed',nrow(df))))
  tmp$group = as.factor(tmp$group)
  p[[letters[i]]] = ggplot(tmp, aes(x=group, y=h2)) + 
    
    geom_violin(data=tmp, aes(x=group, y=h2, fill=group), position=position_dodge(0.6),bw=0.005)+
    geom_boxplot(width=0.2,data=tmp, aes(x=group, y=h2, fill=group), position=position_dodge(0.6))+
    
    #geom_boxplot(outlier.shape = NA)+
    ylim(0,0.15)+
    theme(axis.title.x = element_blank(),
          panel.grid =element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position = 'none')+
    ylab(bquote({h^2})) + 
    geom_hline(yintercept = h2,linetype = 'dashed',color = 'red', size=1)+
    labs(tage=letters[i])
  
}


library(grid)
lay_out = function(...) {    
  x <- list(...)
  n <- max(sapply(x, function(x) max(x[[2]])))
  p <- max(sapply(x, function(x) max(x[[3]])))
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n, p)))    
  
  for (i in seq_len(length(x))) {
    print(x[[i]][[1]], vp = grid::viewport(layout.pos.row = x[[i]][[2]], 
                                           layout.pos.col = x[[i]][[3]]))
  }
} 


pdf('/data/g*/z*/covid19_neanderthals/ld_matrix_projected/results/h2_comparison.pdf',height = 5,width = 8)
output<-lay_out(list(p$a,1,1),
                list(p$b,1,2),
                list(p$c,1,3))
print(output)

dev.off()






