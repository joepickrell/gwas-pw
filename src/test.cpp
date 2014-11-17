/*
 * test.cpp
 *
 *  Created on: Jan 14, 2013
 *      Author: pickrell
 */

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>


#include "SNPs_PW.h"
#include "gwaspw_params.h"
#include "LDmatrix.h"
using namespace std;

int main(){

	Fgwas_params p;
	p.infile = "/Users/jkpickrell/projects/gwas-pw/t.gz";
	p.pairwise = true;
	p.pheno1 = "HDL";
	p.pheno2 = "LDL";
	//p.finemap = true;
	p.cor = 0;
	p.overlap = true;
	p.ldfile = "/Users/jkpickrell/projects/gwas-pw/covariance_matrix/all_chr1_ld";
	p.K = 2;
	p.Nhap = 700;
	SNPs_PW s(&p);

	//pair<double, double> Ri = ld.get_R(d[i].pos, d[j].pos);
	s.llk();
	//cout << s.d[0].BF1+s.d[1].BF2_C(&s.d[0], 0.0235308159933, 0, 0.0745772458041) << "\n";
	//cout <<s.d[0].BF1+s.d[1].BF2 << "\n";
	//cout << s.d[1].BF2 << " BF2 "<< s.d[1].BF2_C(&s.d[0], 0.009, 0, 0.081)  <<" BF2_C\n";
	/*
	boost::numeric::ublas::compressed_matrix<double> m (3, 3, 3 * 3);
	for (unsigned i = 0; i < m.size1 (); ++ i)
		for (unsigned j = 0; j < m.size2 (); ++ j)
			m (i, j) = 3 * i + j;
			//cout << 3*i+j << " "<< m(i,j)<< "\n";
	std::cout << m << std::endl;
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			cout << m(i,j) << "\n";
		}
	}

	vector<bool> an;
	vector<int> dists;
	vector<vector<pair<int, int> > > dmodels;
	SNP_PW s1("rs1", "chr1" , 10000, 3, 10, 0.001, 0.002, an, dists, dmodels, 0.5, 0.5);
	SNP_PW s2("rs2", "chr1" , 20000, 3, 8, 0.001, 0.0025, an, dists, dmodels, 0.5, 0.5);
	cout << s1.BF1 << " "<< s1.BF2 << " "<< s1.BF3 << "\n";

	cout << s1.BF2_C(&s2, 0.5, 0.3, 0.5) << "\n";
	*/
	//vector<int> tmp;
	//tmp.push_back(4123849); tmp.push_back(4131181);
	//LDmatrix ldm("/Users/jkpickrell/projects/gwas-pw/covariance_matrix/all_ldfiles", "chr1", tmp, 100);
	//pair<double, double> test = ldm.get_R(4123849, 4131181);
	//cout << test.first << " "<< test.second << "\n";
	//cout << ldm.get_ld(106048793, 106048793) << "\n" << ldm.get_ld(106049818, 106048793) << "\n";
	//vector<double> m;
	//m.push_back(-2.0);
	//m.push_back(-3.0);
	//s.V1 =  0.0002380952;
	//s.V2 =  0.0002380952;
	//s.W = 0.3;
	//s.Z1 = 10.0;
	//s.Z2 = -10.0;
	//cout << s.calc_logBF1(0) << " "<< s.calc_logBF1(-0.1) << " "<< s.calc_logBF1(0.1) << " "<< s.calc_logBF1(0.2)<< " "<< s.calc_logBF1(-0.2) <<"\n";
	//cout << s.calc_logBF2(0) << " "<< s.calc_logBF2(-0.1) << " "<< s.calc_logBF2(0.1) << " "<< s.calc_logBF2(0.2) << " "<< s.calc_logBF2(-0.2)<< "\n";
	//cout << s.calc_logBF3(0) << " "<< s.calc_logBF3(-0.1) << " "<< s.calc_logBF3(0.1) << " "<< s.calc_logBF3(0.2)<< " "<< s.calc_logBF3(-0.2)<< "\n";
	//Fgwas_params p;
	//p.wannot.push_back("wgEncodeDukeDnase8988T");
	//p.wannot.push_back("wgEncodeDukeDnaseAoSMC");
	//p.cc = false;
	//p.K = 2000;
	//p.finemap = false;
	//p.dannot.push_back("tssdist");
	//p.segannot.push_back("tssdist");
	//p.loquant = 0.33;
	//p.hiquant = 0.67;
	//p.distmodels.push_back("dist_model");
	//string test = "1340.5";
	//cout << atoi(test.c_str()); cout.flush();
	//p.infile = "/Users/pickrell/Documents/workspace/GWAS/test/hb.wpos.wf.COMBINED.wdist.gz";
	//p.infile = "/Users/pickrell/Documents/workspace/GWAS/test/hb_testin.gz";
	//SNPs s(&p);
	//s.GSL_optim();
	//s.cross10(false);

	return 0;
}

