/*
 * LDmatrix.h
 *
 *  Created on: Jun 18, 2014
 *      Author: jkpickrell
 */

#ifndef LDMATRIX_H_
#define LDMATRIX_H_


#include "gwaspw_params.h"
using namespace std;

class LDmatrix{
public:
	LDmatrix();
	~LDmatrix();
	LDmatrix(string, string, vector<int>, int);
	void print();
	void process_infilelist(string);
	void read_matrix();
	string chrom;
	int minpos, maxpos, Nhap;
	double get_ld(int, int);
	pair<double, double> get_R(int, int);
	vector<double> get_hapfreqs(double, double, double, int, int);
	vector<vector<double> > get_cov(vector<double>);
	vector<double> get_delta(vector<double>);
private:
	boost::numeric::ublas::compressed_matrix<double>*m;
	map<int, int> pos2index;
	map<int, int> index2pos;
	set<int> pos2keep;
	vector<string> infilelist;

};


#endif /* LDMATRIX_H_ */
