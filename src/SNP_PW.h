//SNP_PW.h

#ifndef SNP_PW_H_
#define SNP_PW_H_

#include <string>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <sys/stat.h>
#include "gzstream.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <boost/algorithm/string.hpp>
#include "CmdLine.h"
using namespace std;

class SNP_PW{
public:
	SNP_PW();
	SNP_PW(string, string, int, double, double, double, double, vector<bool>, vector<int>, vector<vector<pair<int, int> > >, vector<double>, double); //for pairwise
	string id;
	string chr;
	int pos;

	double Z1; // Z-score
	double Z2;

	double f, N1, N2;
	double V1, V2;
	vector<double> W; //prior variances
	double BF1; // Bayes factor for configuration [1,0]
	double BF2; // [0,1]
	double BF3; // [1,1]

	//corrected BFs including other SNP in LD
	double BF1_C_ind(SNP_PW*, double, pair<double, double>, double, double);
	double BF1_C(SNP_PW*, double, pair<double, double>, double);
	double BF2_C_ind(SNP_PW*, double, pair<double, double>, double, double);
	double BF2_C(SNP_PW*, double, pair<double, double>, double);
	pair<pair<double, double>, pair<double, double> > condZ(SNP_PW*, pair<double, double>, double);
	//double BF3_C(SNP_PW*, double, pair<double, double>);

	//return betas
	double get_beta1();
	double get_beta2();



	int chunknumber;
	float dens;
	vector<bool> annot;
	vector<float> annot_weight;
	vector<int> dists;
	bool condannot;
	void append_distannots(vector<vector<pair<int, int> > >); // convert distances to annotations according to distance models
	double approx_v2();
	double approx_v1();
	double sumlog(double, double);
	double calc_logBF1(double);
	double calc_logBF1_ind(double, double);
	double calc_logBF2(double);
	double calc_logBF2_ind(double, double);
	double calc_logBF3(double);
	double calc_logBF3_ind(double, double);
	int nannot;
	double get_x(vector<double>);
	//double get_x_cond(vector<double>, double);
	double ln_MVN(vector<double>, vector<vector<double> >);
};


#endif /* SNP_PW_H_ */
