## DeepRegionalAnalysis [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/gamazonlab/DeepRegionalAnalysis/blob/master/LICENSE)  

## Reference

Integrative transcriptomic, evolutionary, and causal inference framework for region-level analysis: Application to COVID-19

Dan Zhou<sup>1, 2</sup> and Eric R. Gamazon<sup>1, 2, 3, 4, 5</sup>

<sup>1</sup>Division of Genetic Medicine, Department of Medicine, Vanderbilt University Medical Center, Nashville, TN, USA  
<sup>2</sup>Vanderbit Genetics Institute, Vanderbilt University Medical Center, Nashville, TN, USA  
<sup>3</sup>Data Science Institute, Vanderbilt University Medical Center, Nashville, TN  
<sup>4</sup>Clare Hall, University of Cambridge, Cambridge, United Kingdom  
<sup>5</sup>MRC Epidemiology Unit, University of Cambridge, Cambridge, United Kingdom  

Send correspondence to:
Eric R. Gamazon <ericgamazon@gmail.com>

Code being maintained by:
Dan Zhou <zdangm@gmail.com>



## Example of running SEGMENT-SCAN to estimate expression-mediated h2 local.

Rscript ./SEGMENT-SCAN.r \
--eqtl_path ${example_dir}/df.txt \
--genotype_path ${example_dir}/genotype.txt \
--res_path ${example_dir}/res.txt \
--start_pos=45859651 \
--stop_pos=45909024 \
--n_cases=1610 \
--n_controls=2205

#options  
eqtl_path  Data frame contains eqtl effect and the effect of gene on trait. See the example file (df.txt) for the format. In the example file, "alpha" and "alpha_se" were the estimated effect size and its se of the effect of gene on trait (estimated from MR), respectively.
genotype_path: Genotype in dosage. See the example file (genotype.txt) for the format.  
res_path: Path for the results.  
start_pos: start position for the segment.  
end_pos: end position for the segment.  
n_cases: N of cases in the GWAS.  
n_controls: N of controls in the GWAS.  

#results  
full_n_variants: N of variants in the full cis region of the gene.  
full_h2_local: Estimated expression-mediated h2 for the full cis region.  
segment_n_variants: N of variants in the segment.  
segment_h2_local: Estimated expression-mediated h2 for the segment.  
pi_c: proportion of expression-mediated causal effect explained.  
concentration_ratio: concentration of expression-mediated causal effect. segment_h2_local * full_n_variants / segment_n_variants / full_h2_local






