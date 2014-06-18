/*
 * test.cpp
 *
 *  Created on: Jan 14, 2013
 *      Author: pickrell
 */

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>


#include "SNP_PW.h"
#include "fgwas_params.h"
#include "LDmatrix.h"
using namespace std;
using namespace boost::numeric::ublas;

int main(){
	compressed_matrix<double> m (3, 3, 3 * 3);
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
	SNP_PW s;
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

