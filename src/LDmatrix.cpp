/*
 * LDmatrix.cpp
 *
 *  Created on: Jun 18, 2014
 *      Author: jkpickrell
 */

#include "LDmatrix.h"

LDmatrix::LDmatrix(){

}


LDmatrix::LDmatrix(string infiles, vector<int> pos){

	// give it the positions to calculate LD for
	// also list of input files in a single file
	// each file in infile should have name like [chr].[start].[stop].gz

	pos2keep.clear();
	pos2index.clear();
	index2pos.clear();
	m.clear();
	int index = 0;
	int previous = pos[0];
	for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++){
		pos2keep.insert(*it);
		pos2index.insert(make_pair(*it, index));
		index2pos.insert(make_pair(index, *it));
		if (index > 0) assert(*it > previous);
		index++;
		previous = *it;
	}
	minpos = index2pos[0];
	maxpos = index2pos[pos[pos.size()-1] ];
	int Nsnp = pos.size();
	m.resize(Nsnp, Nsnp, Nsnp*800);



}

void LDmatrix::process_infilelist(string list){
	infilelist.clear();

	ifstream in(list.c_str());
	vector<string> line;
	struct stat stFileInfo;
	int intStat;
	string st, buf;

	intStat = stat(infile.c_str(), &stFileInfo);
	if (intStat !=0){
		std::cerr<< "ERROR: cannot open file " << infile << "\n";
		exit(1);
	}
    while(getline(in, st)){
    	buf.clear();
    	stringstream ss(st);
    	line.clear();
    	while (ss>> buf){
    		line.push_back(buf);
    	}
    	infilelist.push_back(line[0]);
    }
}
