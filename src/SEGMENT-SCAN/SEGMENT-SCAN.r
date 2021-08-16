
#---SEGMENT-SCAN---

#local h2 estimation
#full cis region v.s. segment

library('RiskPortfolios')
library('MASS')
library('Matrix')
library('pdist')
library('corpcor')
library('optparse')

option_list = list(
  make_option("--eqtl_path", action="store", default=NA, type='character',
              help="Path to dataframe of summary statistics [required]"),
  make_option("--genotype_path", action="store", default=NA, type='character',
              help="Path to genotype file (dosage format) [required]"),
  make_option("--res_path", action="store", default=NA, type='character',
              help="Path for the result [required]"),
  
  make_option("--start_pos", action="store", default=NA, type='integer',
              help="segment position"),
  make_option("--stop_pos", action="store", default=NA, type='integer',
              help="segment position"),
  
  
  make_option("--n_cases", action="store", default=NA, type='integer',
              help="n cases (gwas)"),
  make_option("--n_controls", action="store", default=NA, type='integer',
              help="n controls (gwas)")
  
)

opt = parse_args(OptionParser(option_list=option_list))
eqtl_path = opt$eqtl_path
genotype_path = opt$genotype_path
res_path = opt$res_path
start_pos = opt$start_pos
stop_pos = opt$stop_pos
n_cases = opt$n_cases
n_controls = opt$n_controls

#load eqtl and alpha
eqtl = read.table(eqtl_path, header = T,stringsAsFactors = F)

#load genotype in dosage
dosage = read.table(genotype_path, header = T,stringsAsFactors = F)

#in segment or not
eqtl$in_segment<-ifelse((eqtl$pos>start_pos & eqtl$pos<stop_pos),1,0)

#effective sample size
ascertained_prevalence=n_cases/(n_cases+n_controls)
N_effective = 4*ascertained_prevalence*(1-ascertained_prevalence)*(n_cases+n_controls)

#output df
res = data.frame()



#---estimate h2 for the full cis region of the gene---

full<-eqtl

#estimate snp-snp var-cov matrix (projected)
full_c <- covEstimation(as.matrix(dosage), control = list(type = 'cor'))

#Moore-Penrose pseudoinverse
inverse_full_c<-ginv(full_c)

#q=rank(C)  the maximum number of linearly independent vectors in rows or columns
q=rankMatrix(inverse_full_c)
q=as.numeric(q)

#pseudo gwas beta
alpha = full[1,'alpha']
alpha_se = full[1,'alpha_se']
full$beta_gwas = alpha * full$beta_eqtl

#pseudo gwas se
full$se_gwas = ((alpha_se)^2 * (full$se_eqtl)^2 +
                  (alpha_se)^2 * (full$beta_eqtl)^2 +
                  (alpha)^2 * (full$se_eqtl)^2)^0.5

#n
res[1,'full_n_effective'] = n = N_effective
#p
res[1,'full_n_variants'] = p = nrow(full)
#h2
res[1,'full_h2_local'] = (t(full$beta_gwas) %*% inverse_full_c %*% full$beta_gwas - q/n) / (n-q) *n


#---estimate h2 for the segment (mediated by expression of the gene)---

segment<-full[full$in_segment==1,]

#estimate snp-snp var-cov matrix (projected)
segment_c<-full_c[which(full$in_segment==1),which(full$in_segment==1)]

#Moore-Penrose pseudoinverse
inverse_segment_c<-ginv(segment_c)

#q=rank(C)  the maximum number of linearly independent vectors in rows or columns
q=rankMatrix(inverse_segment_c)
q=as.numeric(q)

#pseudo gwas beta
segment$beta_gwas = alpha * segment$beta_eqtl

#pseudo gwas se
segment$se_gwas = ((alpha_se)^2 * (segment$se_eqtl)^2 +
                     (alpha_se)^2 * (segment$beta_eqtl)^2 +
                     (alpha)^2 * (segment$se_eqtl)^2)^0.5

#n
res[1,'segment_n_effective']=n=N_effective
#p
res[1,'segment_n_variants']=p=nrow(segment)
#h2
res[1,'segment_h2_local']=(t(segment$beta_gwas) %*% inverse_segment_c %*% segment$beta_gwas - q/n) / (n-q) *n


res$pi_c=res$segment_h2_local/res$full_h2_local

res$concentration_ratio=res$segment_h2_local*res$full_n_variants/res$segment_n_variants/res$full_h2_local

write.table(res,res_path,quote = F,sep='\t',row.names = F)







