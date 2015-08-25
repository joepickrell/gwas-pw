/*
 * fgwas_params.cpp
 *
 *  Created on: Jun 17, 2013
 *      Author: pickrell
 */

#include "gwaspw_params.h"
using namespace std;

Fgwas_params::Fgwas_params(){
	K = 5000;
	V.push_back(0.01); V.push_back(0.1); V.push_back(0.5);
	print = true;
	zformat = true;
	wannot.clear();
	dannot.clear();
	distmodels.clear();
	segannot.clear();
	outstem = "gwas-pw";
	dropchr = false;
	cc = false;
	finemap = false;
	ridge_penalty = 0.2;
	xv = false;
	onlyp = false;
	cond = false;
	testcond_annot = "";
	pairwise = true;
	pheno1 = "NA";
	pheno2 = "NA";
	burnin = 3000;
	nsamp = 30000;
	overlap = false;
	sampfreq = 10;
	rev = false;
	MCMC_gauss_SD = sqrt(1.0/5.0);
	//alpha_prop = 50.0;
	alpha_prior.push_back(2);
	for (int i = 0; i < 4; i++) {
		alpha_prior.push_back(-2);
	}
	cor = 0;
	Nhap = 0;
	bedseg = false;
	numberedseg = false;
	segment_bedfile = "";
	MCMC = false;
}

void Fgwas_params::print_stdout(){
	cout <<"\n\n";
	cout << ":::Parameter settings::::\n";
	cout << ":: Input file: "<< infile << "\n";
	cout << ":: Output stem: "<< outstem << "\n";
	cout << ":: Phenotype 1: "<< pheno1 << "\n";
	cout << ":: Phenotype 2: "<< pheno2 << "\n";
	if (!bedseg) cout << ":: K: " << K << "\n";
	else cout << ":: Segment bedfile: "<< segment_bedfile << "\n";
	cout << ":: V: ";
	for (int i = 0; i < V.size(); i++) cout << V[i] << " ";
	cout << "\n";
	cout << ":: Fine-mapping?: ";
	if (finemap) cout << "yes\n";
	else cout << "no\n";
	cout << ":: Numbered segments?: ";
	if (numberedseg) cout << "yes\n";
	else cout << "no\n";
	cout << ":: Print: ";
	if (print) cout << "yes\n";
	else cout <<"no\n";
	cout << ":: MCMC: ";
	if (MCMC) cout << "yes\n";
	else cout <<"no\n";
	cout <<":: Correlation: "<< cor << "\n";
	if (overlap){
		cout <<":: Overlapping cohorts: yes\n";
		cout <<":: LD file: "<< ldfile<< "\n";
		cout <<":: Nhap: "<< Nhap << "\n";
	}
	else{
		cout << ":: Overlapping cohorts: no\n";
	}
	//cout << ":: SNP annotations:";
	//for (vector<string>::iterator it = wannot.begin(); it != wannot.end(); it++) cout << " "<< *it; cout << "\n";
	//cout << ":: Distance models:";
	//for (int i = 0; i < dannot.size(); i++)	cout << " " << dannot[i]<< ":" << distmodels[i];  cout << "\n";
	//cout << ":: Segment annotation (low quantile, high quantile):";
	//for (vector<string>::iterator it = segannot.begin(); it != segannot.end(); it++) cout << " "<< *it<<  " ("<< loquant << " "<< hiquant << ")"; cout << "\n";
	//cout << ":: Conditional analysis of: " << testcond_annot << "\n";
	cout << ":::::::::::::::::::::::::\n";
	cout <<"\n\n";
}
