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

-[output].segbfs.gz contains a line for each segment of the genome. The columns are:

1. chunk: the internal numerical identifer for the segment
2. NSNP: the number of SNPs in the segment
3. chr: chromosome 
4. st: star position
5. sp: end position 
6. max_abs_Z_[pheno1]: the maximum absolute value of the Z-score for phenotype 1 in the region 
7. max_abs_Z_[pheno2]: the maximum absolute value of the Z-score for phenotype 2 in the region 
8. logBF_1: ln(regional Bayes factor supporting model 1 [association only to phenotype 1] versus the null)
9. logBF_2: ln(regional Bayes factor supporting model 2 [association only to phenotype 2] versus the null) 
10. logBF_3: ln(regional Bayes factor supporting model 3 [shared association to both phenotypes] versus the null) 
11. logBF_4: ln(regional Bayes factor supporting model 3 [two distinct associations, one to each phenotype] versus the null) 
12. pi_1: prior on model 1 
13. pi_2: prior on model 2
14. pi_3: prior on model 3 
15. pi_4: prior on model 4 
16. PPA_1: posterior probability of model 1 
17. PPA_2: posterior probability of model 2 
18. PPA_3: posterior probability of model 3 
19. PPA_4: posterior probability of model 4

-[output].bfs.gz contains a line for each SNP in the genome. The columns are:

1. id: SNP identifier
2. chr: chromosome 
3: pos: position
4. logBF_1: ln(Bayes factor measure the suppport for model 1 at the SNP) 
5. logBF_2: ln(Bayes factor measure the suppport for model 2 at the SNP) 
6. logBF_3: ln(Bayes factor measure the suppport for model 3 at the SNP) 
7. Z_[pheno1]: Z-score for association to phenotype 1
8. V_[pheno1]: variance in the effect size estimate for phenotype 1
9. Z_[pheno2]: Z-score for association to phenotype 2
10. V_[pheno2]: variance in the effect size estimate for phenotype 2 
11. pi_1: prior on this SNP being the causal one under model 1 
12. pi_2: prior on this SNP being the causal one under model 2 
13. pi_3: prior on this SNP being the causal one under model 3 
14. PPA_1: posterior probability that this SNP is the causal one under model 1
15. PPA_2: posterior probability that this SNP is the causal one under model 2 
16. PPA_3: posterior probability that this SNP is the causal one under model 3
17. chunk: the internal numerical identifer for the segment this SNP falls in

-[output].MLE contains the estimated regional prior probabilites of each model (same as in [output].segbfs.gz)

###Options###

-i [file name] name of the input file, in the format described above

-phenos [string] [string] names of the phenotypes, such the the Z scores are in columns labeled Z_[pheno1] and Z_[pheno2]

-o [string] stem for names of output files

-bed [file name] gwas-pw splits the genome into approximately independent blocks. To input these blocks from a .bed file, use this option. We recommend using the bed files available from https://bitbucket.org/nygcresearch/ldetect-data

-noprint don't print the Bayes factors

-k [integer] as an alternative to spliting the genome into blocks based on the bed file, input the number of SNPs per block. If neither -k or -bed is specified, this defaults to blocks of 5,000 SNPs

-cor [float] if the two GWAS were performed using overlapping cohorts, use this flag to specify the expected correlation in summary statistics under the null (defaults to zero)
