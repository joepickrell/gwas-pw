gwas-pw
=====

gwas-pw is a tool for jointly analysing two genome-wide association studies (GWAS). The basic setup is that you have performed two GWAS and want to identify loci that influence both traits. Instead of using two P-value thresholds to identify variants that influence both traits, the algorithm learns the appropriate thresholds from the data. 

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

