/*
 * LDmatrix.h
 *
 *  Created on: Jun 18, 2014
 *      Author: jkpickrell
 */

#ifndef LDMATRIX_H_
#define LDMATRIX_H_


#include "fgwas_params.h"
using namespace std;

class LDmatrix{
public:
	LDmatrix();
	LDmatrix(string, vector<int>);
	void print();
	void process_infilelist(string);
private:
	boost::numeric::ublas::compressed_matrix<double> m;
	map<int, int> pos2index;
	map<int, int> index2pos;
	set<int> pos2keep;
	int minpos, maxpos;
	vector<string> infilelist;
};


#endif /* LDMATRIX_H_ */
