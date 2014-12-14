/*
 * LDmatrix.cpp
 *
 *  Created on: Jun 18, 2014
 *      Author: jkpickrell
 */

#include "LDmatrix.h"

LDmatrix::LDmatrix(){

}
LDmatrix::~LDmatrix(){
	delete m;
}

LDmatrix::LDmatrix(string infiles, string chr, vector<int> pos, int N){

	// give it the positions to calculate LD for
	// also list of input files in a single file
	// each file in infile should have name like [chr].[start].[stop].gz
	//cout << "here\n";
	pos2keep.clear();
	pos2index.clear();
	index2pos.clear();
	Nhap = N;
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
	assert(minpos <= maxpos);
	int Nsnp = pos.size();
	//cout << "here3\n"; cout.flush();
	m = new boost::numeric::ublas::compressed_matrix<double>(Nsnp, Nsnp, Nsnp*2000);
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
	int n_nonzero = 0;
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
		if (fileinfo.size() < 3){
			cerr << "ERROR: file name of "<<*it<<" don't have [chr].[st].[sp].gz\n";
			exit(1);
		}
        	string tmpchr = fileinfo.at(0);
        	int startpos = atoi(fileinfo.at(1).c_str());
        	int endpos = atoi(fileinfo.at(2).c_str());
        	// does the file match the chromosome [chr].[startpos].[endpos].gz
        	if (tmpchr != chrom) continue;
       		// cout << minpos << " "<< maxpos << "\n";cout.flush();
        	if ( (startpos < minpos && endpos >= minpos) || (startpos >=minpos && startpos < maxpos)){
        		cout << "Reading LD from "<<*it << "..."; cout.flush();
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
            		int pos1 = atoi(line.at(2).c_str());
            		int pos2 = atoi(line.at(3).c_str());
            		double D = atof(line.at(7).c_str());
            		if (pos2index.find(pos1) != pos2index.end() && pos2index.find(pos2) != pos2index.end()){
                		//cout << pos1 << " "<< pos2 << " "<< D << " "; cout.flush();
				int index1 = pos2index[pos1];
                		int index2 = pos2index[pos2];
                		(*m)(index1, index2) = D;
				//cout << "added "<< n_nonzero <<"\n"; cout.flush();
				n_nonzero++;
            		}

            	}
	    	cout << "done\n"; cout.flush();
            }
	
	}
}

pair<double, double> LDmatrix::get_R(int p1, int p2){
	// return R = D/V1
	// and the approximate variance in R.
	// var(R) = \delta \sigma \delta^T, see multivariate delta method in Agresti (2002)

	double D = get_ld(p1, p2);
	double V1 = get_ld(p1, p1);
	double V2 = get_ld(p2, p2);
	double R = D/V1;
	vector<double> hapfreqs = get_hapfreqs(V1, V2, D, p1, p2);
	vector<vector<double> > C = get_cov(hapfreqs);
	vector<double> delta = get_delta(hapfreqs);

	vector<double> dC;
	for (int i = 0; i < 3; i++){
		double tmp =0;
		for (int j = 0; j < 3; j++){
			tmp+= delta[j]*C[i][j];
		}
		dC.push_back(tmp);
	}
	double VR = dC[0]*delta[0] + dC[1]*delta[1]+dC[2]*delta[2];
	//cout << "in LD "<< R << " " << VR << "\n";
	return(make_pair(R, VR));
}

vector<double> LDmatrix::get_delta(vector<double> fs){
	vector<double> toreturn;
	double p = fs[0]+fs[1];
	double q = 1-p;
	double t1 = (fs[1]*q*q - fs[2]*p*p)/ (p*p*q*q);
	double t2 = (-fs[2]*p*p - fs[0]*q*q) / (p*p*q*q);
	double t3 = 1/(p-1);
	toreturn.push_back(t1);toreturn.push_back(t2);toreturn.push_back(t3);
	//cout << "delta "<< t1 << " "<< t2 << " "<< t3 << "\n";
	return toreturn;
}
vector<vector<double> > LDmatrix::get_cov(vector<double> fs){
	vector<vector<double> > toreturn;
	for (int i = 0; i < 3; i++){
		vector<double> tmp;
		for (int i =0 ; i < 3; i++) tmp.push_back(0.0);
		toreturn.push_back(tmp);
	}
	//for (int i = 0; i < 3; i++){
	//	cout << fs[i] << " ";
	//}
	//cout << " f_in_cov\n";
	for (int i = 0; i < 3; i++){
		for (int j = i; j< 3; j++){
			if (i == j) toreturn[i][j] = fs[i] *(1-fs[i])/(double)Nhap;
			else{
				toreturn[i][j] = -1* fs[i]*fs[j] / (double)Nhap;
				toreturn[j][i] = -1* fs[i]*fs[j] / (double)Nhap;
			}
		}
	}
	//for (int i = 0; i < 3; i++){
	//	for (int j = 0; j < 3; j++){
	//		cout << toreturn[i][j] << " ";
	//	}
	//	cout << "\n";
	//}
	return(toreturn);
}


vector<double> LDmatrix::get_hapfreqs(double V1, double V2, double D, int pos1, int pos2){
	double tmpp1_1 = (1.0 + sqrt(1-4.0*V1))/2.0;
	double tmpp1_2 = (1.0 - sqrt(1-4.0*V1))/2.0;

	double tmpp2_1 = (1.0 + sqrt(1-4.0*V2))/2.0;
	double tmpp2_2 = (1.0 - sqrt(1-4.0*V2))/2.0;

	double f11, f10, f01, f00, p1, p2;
	double tmpf11_1 = D+ tmpp1_1* tmpp2_1;
	double tmpf11_2 = D+ tmpp1_2* tmpp2_1;

	//cout << tmpf11_1 << " "<< tmpf11_2 << " "<< tmpp1_1 << " "<< tmpp1_2 << " "<< tmpp2_1 << " "<< tmpp2_2 << "\n";
	if (tmpf11_1 < tmpp1_1 and tmpf11_1 < tmpp2_1 and tmpp1_1 + tmpp2_1 - tmpf11_1 < 1){
		f11 = tmpf11_1;
		p1 = tmpp1_1;
		p2 = tmpp2_1;
		f10 = p1- f11;
		f01 = p2 - f11;
		f00 = 1-f11-f10-f01;
	}
	else if(tmpf11_2 < tmpp1_2 and tmpf11_2 < tmpp2_1 and tmpp1_2 + tmpp2_2 - tmpf11_2 < 1){
		f11 = tmpf11_2;
		p1 = tmpp1_2;
		p2 = tmpp2_1;
		f10 = p1 -f11;
		f01 = p2 - f11;
		f00 = 1-f11-f10-f01;
	}
	// rare pathological cases give negative haplotype frequencies, set those to zero
	if (f11< 0 or f10 < 0 or f01 < 0 or f00 >1 or f11 >1 or f10 > 1 or f01 >1 or f00 < 0){
		cerr << "WARNING: trouble getting haplotype frequencies from V1: "<<V1 <<" ("<< pos1<<") V2: "<< V2 << " ("<< pos2 << ") D: "<< D<< " (negative haplotype fs)\n";
		if (f11 < 0) f11 = 0;
		if (f10 < 0) f10 = 0;
		if (f01 < 0) f01 = 0;
		if (f00 < 0) f00 = 0;
		if (f11 > 1) f11 = 1;
		if (f10 > 1) f10 = 1;
		if (f01 > 1) f01 = 1;
		if (f00 > 1) f00 = 1;
		 //renormalize
		double sum = f11+f10+f01+f00;
		f11 = f11/sum;
		f10 = f10/sum;
		f01 = f01/sum;
		f00 = f00/sum;
		//exit(1);
	}
	vector<double> toreturn;
	//cout << p1 << " "<< p2 << "\n";
	toreturn.push_back(f11); toreturn.push_back(f10); toreturn.push_back(f01); toreturn.push_back(f00);
	//for (vector<double>::iterator it = toreturn.begin(); it != toreturn.end(); it++) cout << *it << " ";
	//cout << " hap\n";
	return(toreturn);

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
