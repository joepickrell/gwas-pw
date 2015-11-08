gwas-pw
=====

gwas-pw is a tool for jointly analysing two genome-wide association studies (GWAS). The basic setup is that you have performed two GWAS and want to identify loci that influence both traits. Instead of using two P-value thresholds to identify variants that influence both traits, the algorithm learns reasonable thresholds from the data. 

###Dependencies###
gwas-pw depends on:

- the [GNU Scientific Library](http://www.gnu.org/software/gsl/)

- the [Boost Libraries](http://www.boost.org)

###Quick Start###
The most up-to-date release is: version 0.21. See ["Releases"](https://github.com/joepickrell/gwas-pw/releases) above.
After downloading gwas-pw-0.21.tar.gz at the link above, run:

>tar -xvf gwas-pw-0.21.tar.gz

>cd gwas-pw-0.21

>./configure

>make

This will create an executable file called gwas-pw in the src directory. The most common compilation error is that the configure script cannot find Boost or GSL. You may have to tell the script explicitly where to find them. For example, on OS X using macports, installations go to the non-standard path /opt/local/lib. To configure in this case, replace the above configure step with:

>./configure LDFLAGS=-L/opt/local/lib

Example data is available in the example_data/ directory. To ensure that gwas-pw is working, run:

>gwas-pw -i example_data/aam_height_example.gz -bed example_data/all_fourier_ls.bed -phenos AAM HEIGHT


###Input file format###
The input file must have the following columns (in any order, they will be identified by the header). Rows must be sorted by chromosomal position:

1. SNPID: A string with a SNP identifier
2. CHR: chromosome
3. POS: position
4. Z_[pheno1]: the signed Z score measuring the evidence for association to phenotype 1 at the SNP
5. V_[pheno1]: the variance in the effect size estimate with phenotype 1 at this SNP
6. Z_[pheno2]: the signed Z score measuring the evidence for association to phenotype 2 at the SNP (note that the allele chosen to set the sign should be identical for the two phenotypes)
7. V_[pheno2]: the variance in the effect size estimate with phenotype 2 at this SNP

Note the [pheno1] and [pheno2] will be supplied by you at the command line.

###Output file format###

There are three output files:

1. [output].segbfs.gz contains a line for each segment of the genome. The columns are:

-chunk: the internal numerical identifer for the segment
-NSNP: the number of SNPs in the segment
-chr: chromosome 
-st: star position
-sp: end position 
-max_abs_Z_[pheno1]: the maximum absolute value of the Z-score for phenotype 1 in the region 
-max_abs_Z_[pheno2]: the maximum absolute value of the Z-score for phenotype 2 in the region 
-logBF_1: ln(regional Bayes factor supporting model 1 [association only to phenotype 1] versus the null)
-logBF_2: ln(regional Bayes factor supporting model 2 [association only to phenotype 2] versus the null) 
-logBF_3: ln(regional Bayes factor supporting model 3 [shared association to both phenotypes] versus the null) 
-logBF_4: ln(regional Bayes factor supporting model 3 [two distinct associations, one to each phenotype] versus the null) 
-pi_1: prior on model 1 
-pi_2: prior on model 2
-pi_3: prior on model 3 
-pi_4: prior on model 4 
-PPA_1: posterior probability of model 1 
-PPA_2: posterior probability of model 2 
-PPA_3: posterior probability of model 3 
-PPA_4: posterior probability of model 4


