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
	SNP_PW(string, string, int, double, double, vector<bool>, vector<int>, vector<vector<pair<int, int> > >); //for pairwise
	string id;
	string chr;
	int pos;

	double Z; // Z-score
	double BF; // Bayes factor
	double BF2; // second BF for pairwise
	int chunknumber;
	float dens;
	vector<bool> annot;
	vector<float> annot_weight;
	vector<int> dists;
	bool condannot;
	void append_distannots(vector<vector<pair<int, int> > >); // convert distances to annotations according to distance models
	int nannot;
	double get_x(vector<double>);
	double get_x_cond(vector<double>, double);
};


#endif /* SNP_PW_H_ */
