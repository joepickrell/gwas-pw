/*
 * test.cpp
 *
 *  Created on: Jan 14, 2013
 *      Author: pickrell
 */

#include "SNP_PW.h"
#include "fgwas_params.h"
using namespace std;

int main(){

	SNP_PW s;
	vector<double> m;
	m.push_back(0.5);
	m.push_back(-3);
	vector<vector<double> > C;
	vector<double> tmp, tmp2;
	tmp.push_back(0.4);
	tmp.push_back(0.1);
	tmp2.push_back(0.1);
	tmp2.push_back(0.4);
	C.push_back(tmp);
	C.push_back(tmp2);

	cout << s.ln_MVN(m, C) << "\n";
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

