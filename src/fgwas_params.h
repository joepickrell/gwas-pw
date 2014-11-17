/*
 * fgwas_params.h
 *
 *  Created on: Jun 17, 2013
 *      Author: pickrell
 */

#ifndef FGWAS_PARAMS_H_
#define FGWAS_PARAMS_H_


#include <string>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <sys/stat.h>
#include "gzstream.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <boost/algorithm/string.hpp>
#include "CmdLine.h"
#include <assert.h>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/tokenizer.hpp>
using namespace std;

class Fgwas_params{
public:
	Fgwas_params();
	int K; //block size
	vector<double> V; //prior variances
	bool print, zformat;
	string infile, outstem, ldfile;
	bool cc;
	//annotation lists
	vector<string> wannot, dannot, distmodels, segannot;
	double loquant, hiquant;

	//drop chromosomes
	bool dropchr;
	string chrtodrop;
	void print_stdout();
	bool finemap;
	double ridge_penalty;
	bool xv;
	bool onlyp;
	bool cond;
	string testcond_annot;
	bool pairwise;
	string pheno1, pheno2; //for pairwise
	int burnin;
	int nsamp;
	int sampfreq;
	vector<double> alpha_prior;
	int seed;
	double MCMC_gauss_SD;
	double cor;
	bool overlap;
	bool MCMC, bedseg;
	string segment_bedfile;
	int Nhap;
	bool rev;
	bool numberedseg;
};


#endif /* FGWAS_PARAMS_H_ */
