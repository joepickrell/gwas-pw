/*
 * LDmatrix.cpp
 *
 *  Created on: Jun 18, 2014
 *      Author: jkpickrell
 */

#include "LDmatrix.h"

LDmatrix::LDmatrix(){

}


LDmatrix::LDmatrix(string infiles, string chr, vector<int> pos){

	// give it the positions to calculate LD for
	// also list of input files in a single file
	// each file in infile should have name like [chr].[start].[stop].gz
	//cout << "here\n";
	pos2keep.clear();
	pos2index.clear();
	index2pos.clear();
	//m.clear();
	//cout << "here2\n"; cout.flush();
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
	chrom = chr;
	minpos = index2pos[0];
	maxpos = index2pos[pos.size()-1];
	assert(minpos < maxpos);
	int Nsnp = pos.size();
	//cout << "here3\n"; cout.flush();
	m = new boost::numeric::ublas::compressed_matrix<double>(Nsnp, Nsnp, Nsnp*800);
	//cout << "here4\n"; cout.flush();
	process_infilelist(infiles);
	read_matrix();
	//cout << *m << "\n";

}

void LDmatrix::process_infilelist(string list){
	infilelist.clear();
	ifstream in(list.c_str());
	vector<string> line;
	struct stat stFileInfo;
	int intStat;
	string st, buf;

	intStat = stat(list.c_str(), &stFileInfo);
	if (intStat !=0){
		std::cerr<< "ERROR: cannot open file " << list << "\n";
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

void LDmatrix::read_matrix(){
	for (vector<string>::iterator it = infilelist.begin(); it != infilelist.end(); it++){
		//cout << *it << "\n";
        typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
        boost::char_separator<char> sep("/");
        tokenizer tokens(*it, sep);
        vector<string> splitfile;
        for (tokenizer::iterator tok_iter = tokens.begin();  tok_iter != tokens.end(); ++tok_iter){
        	string tmp = *tok_iter;
        	splitfile.push_back(tmp);
        }
       // for (vector<string>::iterator it2 = splitfile.begin(); it2 != splitfile.end(); it2++) cout << *it2 << " ";
       // cout <<  "\n";
        string filename = splitfile.at(splitfile.size()-1);
        boost::char_separator<char> sep2(".");
        tokenizer tokens2(filename, sep2);
        vector<string> fileinfo;
        for (tokenizer::iterator tok_iter = tokens2.begin();  tok_iter != tokens2.end(); ++tok_iter){
        	string tmp = *tok_iter;
        	fileinfo.push_back(tmp);
        }
        string tmpchr = fileinfo[0];
        int startpos = atoi(fileinfo[1].c_str());
        int endpos = atoi(fileinfo[2].c_str());
        // does the file match the chromosome [chr].[startpos].[endpos].gz
        if (tmpchr != chrom) continue;
       // cout << minpos << " "<< maxpos << "\n";cout.flush();
        if ( (startpos < minpos && endpos >= minpos) || (startpos >=minpos && startpos < maxpos)){
        	cout << "Reading LD from "<<*it << "\n";
            //cout << startpos << " "<< endpos << " here\n";
        	igzstream in(it->c_str());
        	vector<string> line;
        	struct stat stFileInfo;
        	int intStat;
        	string st, buf;

        	intStat = stat(it->c_str(), &stFileInfo);
        	if (intStat !=0){
        		std::cerr<< "ERROR: cannot open file " << *it << "\n";
        		exit(1);
        	}
            while(getline(in, st)){
            	buf.clear();
            	stringstream ss(st);
            	line.clear();
            	while (ss>> buf){
            		line.push_back(buf);
            	}
            	int pos1 = atoi(line[2].c_str());
            	int pos2 = atoi(line[3].c_str());
            	double D = atof(line[7].c_str());
            	if (pos2index.find(pos1) != pos2index.end() && pos2index.find(pos2) != pos2index.end()){
                	int index1 = pos2index[pos1];
                	int index2 = pos2index[pos2];
                	(*m)(index1, index2) = D;
            	}

            }
        }
	}
}

double LDmatrix::get_ld(int p1, int p2){
	assert(pos2index.find(p1) != pos2index.end());
	assert(pos2index.find(p2) != pos2index.end());
   	if (pos2index.find(p1) != pos2index.end() && pos2index.find(p2) != pos2index.end()){
   		int index1 = pos2index[p1];
   		int index2 = pos2index[p2];
   		//only storing one of (i, j) and (j, i), so need to check both
   		double t1 = (*m)(index1, index2);
   		double t2 = (*m)(index2, index1);
   		if (index1 == index2 and t1 < 1e-8){
   			cerr<< "ERROR: cannot find position "<< p1 << " in covariance matrix\n";
   			exit(1);
   		}
   		if (fabs(t1)>fabs(t2)) return t1;
   		else return t2;

   	}
   	else{
   		cout << "ERROR: looking for "<< p1 << " and "<< p2 <<" in LDmatrix, don't exist";
   		exit(1);
   	}
}
