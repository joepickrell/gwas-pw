
bin_PROGRAMS = gwas-pw test
DISTCHECK_CONFIGURE_FLAGS=LDFLAGS=-L/opt/local/lib CXXFLAGS=-I/opt/local/include
gwas_pw_SOURCES = CmdLine.cpp gwas-pw.cpp gzstream.cpp SNP_PW.cpp SNPs_PW.cpp gwaspw_params.cpp LDmatrix.cpp
test_SOURCES = LDmatrix.cpp test.cpp gzstream.cpp gwaspw_params.cpp SNP_PW.cpp SNPs_PW.cpp
noinst_HEADERS = LDmatrix.h gzstream.h CmdLine.h gwaspw_params.h SNP_PW.h SNPs_PW.h
