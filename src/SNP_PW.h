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
	SNP_PW(string, string, int, double, double, double, double, vector<bool>, vector<int>, vector<vector<pair<int, int> > >, double, double); //for pairwise
	string id;
	string chr;
	int pos;

	double Z1; // Z-score
	double Z2;

	double f, N1, N2;
	double V1, V2, W;
	double BF1; // Bayes factor
	double BF2; // second BF for pairwise
	double BF3; // third BF for both

	//corrected BFs including other SNP in LD
	double BF1_C(double, double, double, double, double);
	double BF2_C(double, double, double, double, double);
	double BF3_C(double, double, double, double, double);

	int chunknumber;
	float dens;
	vector<bool> annot;
	vector<float> annot_weight;
	vector<int> dists;
	bool condannot;
	void append_distannots(vector<vector<pair<int, int> > >); // convert distances to annotations according to distance models
	double approx_v2();
	double approx_v1();
	double calc_logBF1(double);
	double calc_logBF2(double);
	double calc_logBF3(double);
	int nannot;
	double get_x(vector<double>);
	//double get_x_cond(vector<double>, double);
	double ln_MVN(vector<double>, vector<vector<double> >);
};


#endif /* SNP_PW_H_ */
