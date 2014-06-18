
bin_PROGRAMS = gwas-pw test
DISTCHECK_CONFIGURE_FLAGS=LDFLAGS=-L/opt/local/lib CXXFLAGS=-I/opt/local/include
gwas_pw_SOURCES = CmdLine.cpp gwas-pw.cpp gzstream.cpp SNP_PW.cpp SNPs_PW.cpp fgwas_params.cpp LDmatrix.cpp
test_SOURCES = test.cpp gzstream.cpp fgwas_params.cpp SNP_PW.cpp
noinst_HEADERS = gzstream.h CmdLine.h fgwas_params.h SNP_PW.h SNPs_PW.h
