/*
 * SNP_PW.cpp
 *
 */
#include "SNP_PW.h"
using namespace std;

SNP_PW::SNP_PW(){

}



SNP_PW::SNP_PW(string rs, string c, int p, double b1, double b2, vector<bool> an, vector<int> ds, vector<vector<pair<int, int> > > dmodels){
	//for pairwise
	id = rs;
	chr = c;
	pos = p;
	BF = b1;
	BF2 = b2;
	for (vector<bool>::iterator it = an.begin(); it != an.end(); it++) {
		annot.push_back(*it);
		if (*it) annot_weight.push_back(1.0);
		else annot_weight.push_back(0.0);
	}

	//distribute weights
	float s = 0;
	for (vector<float>::iterator it = annot_weight.begin(); it != annot_weight.end(); it++) s += *it;
	if (s > 0) for (int i = 0; i < annot_weight.size(); i++) annot_weight[i]  = annot_weight[i]/s;

	for (vector<int>::iterator it = ds.begin(); it != ds.end(); it++) dists.push_back(*it);
	// append distance annotations
	append_distannots(dmodels);
	nannot = annot.size();
}


void SNP_PW::append_distannots(vector<vector<pair<int, int> > > dmodels){
	for (int i = 0; i < dists.size(); i++){
		int dist = dists[i];
		bool found = false;
		vector<pair<int, int> > model = dmodels[i];
		for (vector<pair<int, int> >::iterator it = model.begin(); it != model.end(); it++){
			int st = it->first;
			int sp = it->second;
			if (dist >= st && dist < sp) {
				if (found){
					cerr << "ERROR: SNP "<< id << " is in more than one distance bin for distance measure number "<< i << "\n";
					exit(1);
				}
				annot.push_back(true);
				annot_weight.push_back(1.0);
				found = true;
			}
			else {
				annot.push_back(false);
				annot_weight.push_back(0.0);
			}
		}
	}
}

double SNP_PW::get_x(vector<double> lambda){
	if (lambda.size() != nannot){
		cerr << "ERROR: SNP "<< id << ". Lambda has "<< lambda.size()<< " entries. nannot is " << nannot << "\n";
		exit(1);
	}
	if (nannot == 0) return 0;
	double toreturn = 0;
	for (int i = 0; i < nannot; i++) {
		if (annot[i]) toreturn += lambda[i];
	}
	return toreturn;
}

double SNP_PW::get_x_cond(vector<double> lambda, double lambdac){
	if (lambda.size() != nannot){
		cerr << "ERROR: SNP "<< id << ". Lambda has "<< lambda.size()<< " entries. nannot is " << nannot << "\n";
		exit(1);
	}
	if (nannot == 0) return 0;
	double toreturn = 0;
	for (int i = 0; i < nannot; i++) {
		if (annot[i]) toreturn += lambda[i];
	}
	if (condannot) toreturn+= lambdac;
	return toreturn;
}

