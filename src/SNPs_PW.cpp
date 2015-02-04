/*
 * SNPs.cpp
 *
 *  Created on: Jan 14, 2013
 *      Author: pickrell
 */

#include "SNPs_PW.h"
using namespace std;

SNPs_PW::SNPs_PW(){

}


SNPs_PW::SNPs_PW(Fgwas_params *p){
	params = p;
	params->print_stdout();
	precomputed = false;
	//read distance models
	for (vector<string>::iterator it = params->distmodels.begin(); it != params->distmodels.end(); it++) dmodels.push_back( read_dmodel(*it));

	//read input file
	if (params->pairwise) load_snps_pw(params->infile, params->wannot, params->dannot, params->segannot);
	else{
		cerr<< "ERROR: in pairwise, parameter file says this isn't a pairwise run\n";
		exit(1);
	}
	//make segments
	if (params->finemap || params->numberedseg) make_segments_finemap();
	else{
		make_chrsegments();
		if (params->bedseg) make_segments(params->segment_bedfile);
		else make_segments(params->K);
	}

	//double-check input quality
	check_input();



	//initialize
	snppri.clear();
	pi.clear();
	pi.push_back(1);pi.push_back(1);pi.push_back(1);pi.push_back(1);pi.push_back(1);
	init_segpriors();
    phi = (1+sqrt(5))/2;
    resphi = 2-phi;
	for (int i = 0; i < d.size(); i++){
		vector<double> sp;
		sp.push_back(1.0);sp.push_back(1.0);sp.push_back(1.0);
		snppri.push_back(sp);
		snppost.push_back(1.0);
	}

	nannot = annotnames.size();
	nsegannot = segannotnames.size();
	for (int i = 0; i < nannot; i++)	lambdas.push_back(0);
	for (int i = 0; i < nsegannot; i++) seglambdas.push_back(0);
	set_priors();
}

void SNPs_PW::check_input(){

	double meansize = 0.0;
	double meannsnp = 0.0;
	int minnsnp = 1000000000;
	int maxsnp = 0;
	int nseg = segments.size();
	for (vector<pair<int, int> >::iterator it = segments.begin(); it != segments.end(); it++){

		int st = it->first;
		int sp = it->second;
		int toadd = d[sp-1].pos- d[st].pos;
		int l = sp-st;
		if (l < minnsnp) minnsnp = l;
		if (l > maxsnp) maxsnp = l;
		meannsnp += (double) l / (double) nseg;

		meansize += ((double) toadd/ 1000000.0) / (double) nseg;

		int prevpos = d[st].pos;
		string prevchr = d[st].chr;
		for (int i= st+1; i < sp; i++){
			string testchr = d[i].chr;
			int testpos = d[i].pos;
			if (testchr == prevchr and prevpos >= testpos){  //test that each segment is only a single chromosome, is ordered
				cerr<< "ERROR: SNPs out of order\nChromosome "<<testchr << ". Position "<< prevpos << " seen before "<< testpos<< "\n";
				exit(1);
			}
			prevpos = testpos;
			prevchr = testchr;
		}
	}

	cout << "Number of segments: "<< segments.size()<< "\nMean segment size: "<< meansize<< " Mb ("<< meannsnp << " SNPs)\n";
	cout << "Minimum number of SNPs/segment: "<< minnsnp << "\n";
	cout << "Maximum number of SNPs/segment: "<< maxsnp << "\n";
	if (meansize > 10.0){
		cout << "\n****\n**** WARNING: mean segment size is over 10Mb, this often causes convergence problems (in human data). Consider reducing window size (using -k).\n****\n\n"; cout.flush();
	}
}


map<string, vector<pair< int, int> > > SNPs_PW::read_bedfile(string bedfile){
	map<string, vector<pair< int, int> > > toreturn;
	ifstream in(bedfile.c_str());
	vector<string> line;
	struct stat stFileInfo;
	int intStat;
	string st, buf;

	intStat = stat(bedfile.c_str(), &stFileInfo);
	if (intStat !=0){
		std::cerr<< "ERROR: cannot open file " << bedfile << "\n";
		exit(1);
	}
    while(getline(in, st)){
    	buf.clear();
    	stringstream ss(st);
    	line.clear();
    	while (ss>> buf){
    		line.push_back(buf);
    	}
    	string chr = line[0];
    	int start = atoi(line[1].c_str());
    	int stop = atoi(line[2].c_str());
    	if (toreturn.find(chr) == toreturn.end()) {
    		vector<pair<int, int> > tmp;
    		toreturn.insert(make_pair(chr, tmp));
    	}
    	toreturn[chr].push_back(make_pair(start, stop));

    }

	// check to make sure all entries are ordered

	for (map<string, vector<pair<int, int> > >::iterator it = toreturn.begin(); it != toreturn.end(); it++){
		vector<pair<int, int> > tmp = it->second;
		pair<int, int> prev = tmp[0];
		if (prev.second < prev.first) {
				cerr<<"ERROR: start and stop positions in bed file out of order: "<< it->first << " "<< prev.first<< " "<< prev.second << "\n";
				exit(1);
		}
		for (int i = 1; i < tmp.size(); i++){
			pair<int, int> current = tmp[i];
			//cout << current.first << " "<< current.second << " "<< prev.second << "\n";
			if (current.second < current.first) {
					cerr<<"ERROR: start and stop positions in bed file out of order: "<< it->first << " "<< current.first<< " "<< current.second << "\n";
					exit(1);
			}
			if (current.first < prev.second) {
						cerr<<"ERROR: bed file out of order: "<< it->first << " "<< prev.first<< " "<< prev.second << "\n" << it->first << " "<< current.first<< " "<< current.second << "\n";
						exit(1);
			}
			prev =tmp[i];
		}
	}
	return toreturn;

}

void SNPs_PW::make_segments(string bedfile){
	segments.clear();
	map<string, vector<pair<int, int> > > bedsegs = read_bedfile(bedfile);
	for (int i = 0; i < chrnames.size(); i++){
		string tmpchr = chrnames[i];
		//cout << tmpchr << "\n";
		if (bedsegs.find(tmpchr) == bedsegs.end()){
			cerr << "ERROR: chromsome "<< tmpchr << " not found in .bed file\n";
			exit(1);
		}
		vector<pair<int, int> > intervals = bedsegs[tmpchr];
		pair<int, int> chromosome = chrsegments[i];
		pair<int, int> currentseg = intervals[0];
		int j = chromosome.first;
		int start = j;
		int intervalindex = 0;
		while (j < chromosome.second){
			int jpos = d[j].pos;

			if (jpos < currentseg.first){
				//cout<< jpos << " "<< currentseg.first << " "<< currentseg.second << " here1\n";
				cerr << "ERROR: current segment is "<< currentseg.first << " "<< currentseg.second << ", position is "<< jpos << "\n";
				exit(1);
			}
			else if (jpos > currentseg.second and j == start){
				//cout<< jpos << " "<< currentseg.first << " "<< currentseg.second << " here2\n";
				intervalindex++;
				if (intervalindex >= intervals.size()){
					cerr << "ERROR: position "<< jpos << " is outside range of bed file\n";
					exit(1);
				}
				currentseg = intervals[intervalindex];
			}
			else if (jpos >=currentseg.first and jpos < currentseg.second) {
				//cout<< jpos << " "<< currentseg.first << " "<< currentseg.second << " here3\n";
				j = j+1;
			}
			else if (jpos >= currentseg.second and j != start){
				//cout<< jpos << " "<< currentseg.first << " "<< currentseg.second << " here4\n";
				//cout << start << " "<< j << " adding\n";
				segments.push_back(make_pair(start, j));
				start = j;
				intervalindex++;
				if (intervalindex >= intervals.size()){
					cerr << "ERROR: position "<< jpos << " is outside range of bed file\n";
					exit(1);
				}
				currentseg = intervals[intervalindex];
				//j = j+1;
			}
			else{
				cerr << "ERROR: missing something when reading bed? "<< jpos << " "<<  currentseg.first << " "<< currentseg.second << "\n";
				exit(1);
			}


		}
		//cout << start << " "<< j << " done adding\n";
		segments.push_back(make_pair(start, j));
	}
}
vector<pair<int, int> > SNPs_PW::read_dmodel(string infile){
	vector<pair<int, int> > toreturn;
	ifstream in(infile.c_str());
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
    	int start = atoi(line[0].c_str());
    	int stop = atoi(line[1].c_str());
    	assert(stop < start);
    	if (stop < start){
    		cerr <<"ERROR: in distance model "<< infile<< " , " << start << " is after "<< stop << "\n";
    		exit(1);
    	}
    	toreturn.push_back(make_pair(start, stop));
    }
	return toreturn;
}


void SNPs_PW::init_segpriors(){
	//segannot.clear();
	segpriors.clear();
	alpha.clear();
	for (int i = 0; i < 5; i++) {
		alpha.push_back(params->alpha_prior[i]);
	}
	//alpha.push_back(10);alpha.push_back(10);alpha.push_back(10);alpha.push_back(10);alpha.push_back(10);
	//segannotnames.clear();
	//vector<double> segmeans;
	//vector<double> means2sort;

	for (vector<pair<int, int> >::iterator it = segments.begin(); it != segments.end(); it++){
		vector<double> sp;
		sp.push_back(0.2);sp.push_back(0.2);sp.push_back(0.2);sp.push_back(0.2);sp.push_back(0.2);
		segpriors.push_back(sp);
		/*
		if (params->segannot.size() > 0){
			double segmean = 0.0;
			int total = 0;
			for (int i = it->first; i < it->second; i++){
				segmean += d[i].dens;
				total ++;
			}
			segmeans.push_back( segmean / (double) total);
			means2sort.push_back( segmean / (double) total);
		}
		*/
	}


	/*
	if (params->segannot.size() < 1) return;

	segannotnames.push_back(params->segannot[0] +"_lo");
	segannotnames.push_back(params->segannot[0] +"_hi");
	sort(means2sort.begin(), means2sort.end());
	double locutoff = means2sort[floor( (double) means2sort.size()* params->loquant ) ];
	double hicutoff = means2sort[floor( (double) means2sort.size()* params->hiquant ) ];
	for (int i = 0; i < segments.size(); i++){
		vector<bool> annots;
		double m = segmeans[i];
		if (m <= locutoff) annots.push_back(true);
		else annots.push_back(false);
		if (m >= hicutoff) annots.push_back(true);
		else annots.push_back(false);
		segannot.push_back(annots);
	}
	*/
}

void SNPs_PW::load_snps_pw(string infile, vector<string> annot, vector<string> dannot, vector<string> segannot){
	igzstream in(infile.c_str()); //only gzipped files
	vector<string> line;
	struct stat stFileInfo;
	int intStat;
	string st, buf;

	intStat = stat(infile.c_str(), &stFileInfo);
	if (intStat !=0){
		std::cerr<< "ERROR: cannot open file " << infile << "\n";
		exit(1);
	}

	// read header
	getline(in, st);
	buf.clear();
	stringstream ss(st);
	line.clear();
	while (ss>> buf){
		line.push_back(buf);
	}

	//make a map of header to index
	map<string, int> header_index;
	for (int i = 0; i < line.size(); i++){
		header_index[line[i]] = i;
	}
	// get the indices of the annotations
	vector<int> annot_index;
   	for (vector<string>::iterator it = annot.begin(); it != annot.end(); it++){
   		int i = 0;
   		bool found = false;
   		while (i < line.size() and !found){
   			if (line[i] == *it) {
   				annot_index.push_back(i);
   				found = true;
   			}
   			i++;
   		}
   		if (!found){
   			cerr << "ERROR: cannot find annotation "<< *it << "\n";
   			exit(1);
   		}
    	annotnames.push_back(*it);
   	}

   	// get indices of distance annotations
   	vector<int> dannot_index;
   	for (int j = 0; j < dannot.size(); j++){

   		string jname = dannot[j];
   		int i = 0;
   		bool found = false;
   		while (i < line.size() and !found){
   			if (line[i] == jname) {
   				dannot_index.push_back(i);
   				found = true;
   			}
   			i++;
   		}
   		if (!found){
   			cerr << "ERROR: cannot find distance annotation "<< jname << "\n";
   			exit(1);
   		}
   		append_dannotnames(jname, dmodels[j]);
   	}
   	// get indices for segannot
   	int segannotindex;
   	if (segannot.size() > 0){
   		if (header_index.find(segannot[0]) == header_index.end()){
   			cerr << "ERROR: cannot find segment annotation "<< segannot[0] << "\n";
   			exit(1);
   		}
   		segannotindex = header_index[segannot[0]];
   	}
   	// get indices for the rs, maf, chr, pos, N, Ncase, Ncontrol,
   	int rsindex, chrindex, posindex, segnumberindex, condindex, Z1index, Z2index, V1index, V2index;

   	if (header_index.find("SNPID") == header_index.end()){
   		cerr << "ERROR: cannot find SNPID in header\n";
   		exit(1);
   	}
   	else rsindex = header_index["SNPID"];

 	if (header_index.find("CHR") == header_index.end()){
   		cerr << "ERROR: cannot find CHR in header\n";
   		exit(1);
   	}
   	else chrindex = header_index["CHR"];

 	if (header_index.find("POS") == header_index.end()){
   		cerr << "ERROR: cannot find POS in header\n";
   		exit(1);
   	}
   	else posindex = header_index["POS"];

 	if (header_index.find("Z_"+params->pheno1) == header_index.end()){
 		cerr << "ERROR: cannot find Z_"+params->pheno1+" in header\n";
 		exit(1);
 	}
	else Z1index = header_index["Z_"+params->pheno1];

	if (header_index.find("V_"+params->pheno1) == header_index.end()){
 		cerr << "ERROR: cannot find V_"+params->pheno1+" in header\n";
 		exit(1);
 	}
	else V1index = header_index["V_"+params->pheno1];

	if (header_index.find("Z_"+params->pheno2) == header_index.end()){
 		cerr << "ERROR: cannot find Z_"+params->pheno2+" in header\n";
 		exit(1);
 	}
	else Z2index = header_index["Z_"+params->pheno2];
	if (header_index.find("V_"+params->pheno2) == header_index.end()){
 		cerr << "ERROR: cannot find V_"+params->pheno2+" in header\n";
 		exit(1);
 	}
	else V2index = header_index["V_"+params->pheno2];

	if (header_index.find("SEGNUMBER") == header_index.end() && (params->finemap || params->numberedseg)){
   		cerr << "ERROR: cannot find SEGNUMBER in header\n";
   		exit(1);
   	}
   	else segnumberindex = header_index["SEGNUMBER"];

	if (params->cond && header_index.find(params->testcond_annot) == header_index.end()){
		cerr << "ERROR: cannot find annotation "<< params->testcond_annot << "\n";
		exit(1);
	}
	else if (params->cond) condindex = header_index[params->testcond_annot];
	string oldchr = "NA";
    while(getline(in, st)){
    	buf.clear();
    	stringstream ss(st);
    	line.clear();
    	while (ss>> buf){
    		line.push_back(buf);
    	}
    	string rs = line[rsindex];
    	string chr = line[chrindex];

    	if (line[Z1index] == "NA" || line[Z2index]==  "NA") continue;
    	double Z1 = atof(line[Z1index].c_str());
    	double Z2 = atof(line[Z2index].c_str());
    	double V1 = atof(line[V1index].c_str());
    	double V2 = atof(line[V2index].c_str());
    	assert(V1 >0);
    	assert(V2 >0);
    	if (chr != oldchr) oldchr = chr;
    	if (params->dropchr and chr == params->chrtodrop) continue;
    	int pos = atoi(line[posindex].c_str());
    	vector<bool> an;
    	vector<int> dists;
    	for (vector<int>::iterator it = annot_index.begin(); it != annot_index.end(); it++){
    		assert(line[*it] == "1" || line[*it] == "0");
    		if (line[*it] == "1") an.push_back(true);
    		else if (line[*it] == "0") an.push_back(false);
    		else{
    			cerr << "ERROR: only 0 and 1 allowed for annotations, found "<< line[*it] <<"\n";
    			exit(1);
    		}
    	}
    	for (vector<int>::iterator it = dannot_index.begin(); it != dannot_index.end(); it++){
    		dists.push_back( atoi(line[*it].c_str()));
    	}

    	SNP_PW s(rs, chr , pos, Z1, Z2, V1, V2, an, dists, dmodels, params->V, params->cor);
    	if (params->finemap || params->numberedseg){
    		int snumber = atoi(line[segnumberindex].c_str());
    		s.chunknumber = snumber;
    	}
    	if (segannot.size() > 0) s.dens = atof(line[segannotindex].c_str());


    	if (params->cond){
    		if (line[condindex] == "1") s.condannot =true;
    		else if (line[condindex] == "0") s.condannot = false;
    		else{
    			cerr << "ERROR: only 0 and 1 allowed for annotations, found "<< line[condindex] <<"\n";
    			exit(1);
    		}
    	}
    	d.push_back(s);
    }
}


void SNPs_PW::append_dannotnames(string name, vector<pair<int, int> > model){
	for (vector<pair<int, int> >::iterator it = model.begin(); it != model.end(); it++){
		stringstream ss;
		ss << name << "_" << it->first << "_"<< it->second;
		annotnames.push_back(ss.str());
	}
}


vector<pair< pair<int, int>, pair<double, double> > > SNPs_PW::get_cis(){
	vector<pair<pair<int, int>, pair<double, double> > > toreturn;
	vector<double> startalphas;
	for (vector<double>::iterator it = alpha.begin(); it != alpha.end(); it++) startalphas.push_back(*it);
	int sti = 0;
	if (params->finemap) sti = 1;
	for (int i = sti; i < alpha.size(); i++){
		toreturn.push_back(get_cis_alpha(i));
		alpha[i] = startalphas[i];
		set_priors();
	}
	return toreturn;
}



pair< pair<int, int>, pair<double, double> > SNPs_PW::get_cis_alpha(int which){
	pair< pair<int, int>, pair<double, double> > toreturn;
	double tau = 0.001;
	double startlk = llk();
	double thold = startlk - 2;

	double min = -20.0;
	double max = 20.0;
	double test = alpha[which];
	//cout <<  startlk << " "<< thold << " "<< test<< "\n";
	//for (int i = 0; i < 5; i++) cout << alpha[i] <<  " "<< pi[i] << " "<<  i << "\n";
	if (test > max) max = test+20.0;
	if (test < min) min = test-20.0;
	//if (max < 0) max = 20.0;
	//if (min > 0) min = -20.0;

	//upper
	double hi;
	double lo;
	int convhi = 1;
	int convlo = 1;
	alpha[which] = max;
	set_priors();

	if (llk() > thold) hi = pi[which];
	else{
		alpha[which] = test;
		set_priors();
		double start = (test+max)/2;

		convhi = golden_section_alpha_ci(test, start, max, tau, which, thold);
		hi = pi[which];
	}
	cout << hi << " "<< llk() << " hi\n";
	//for (int i = 0; i < 5; i++) cout << alpha[i] <<  " "<< pi[i] << " "<<  i << "\n";
	//lower

	alpha[which] = test;
	set_priors();
	//cout << llk()<< " " << startlk << " "<< thold << "\n";
	alpha[which] = min;
	set_priors();
	if (llk() > thold) lo = pi[which];
	else{
		alpha[which] = test;
		set_priors();
		double start = (test+min)/2;
		convlo = golden_section_alpha_ci(min, start, test, tau, which, thold);
		lo = pi[which];
	}
	cout << lo << " "<< llk() << " lo\n";
	pair<int, int> conv = make_pair(convlo, convhi);
	pair<double, double> ci = make_pair(lo, hi);
	return make_pair(conv, ci);
}


/*

pair< pair<int, int>, pair<double, double> > SNPs::get_cis_condlambda(){
	pair< pair<int, int>, pair<double, double> > toreturn;
	double startlk = llk();
	double thold = startlk - 2;
	cout <<  startlk << " "<< thold << "\n";
	double min = -20.0;
	double max = 20.0;
	double test = condlambda;
	if (test > max) max = test+20.0;
	if (test < min) min = test-20.0;
	if (max < 0) max = 20.0;
	if (min > 0) min = -20.0;
	double start = (test+max)/2;
	double tau = 0.001;
	int convhi = golden_section_condlambda_ci(test, start, max, tau, thold);
	double hi = condlambda;
	cout << hi << " "<< llk() << " hi\n";
	start = (test+min)/2;
	int convlo = golden_section_condlambda_ci(min, start, test, tau, thold);
	double lo = condlambda;
	cout << lo << " "<< llk() << " lo\n";
	pair<int, int> conv = make_pair(convlo, convhi);
	pair<double, double> ci = make_pair(lo, hi);
	return make_pair(conv, ci);
}
*/

/*
pair< pair<int, int>, pair<double, double> > SNPs::get_cis_segpi(){
	pair< pair<int, int>, pair<double, double> > toreturn;
	double startsegpi = segpi;
	double startlk = llk();
	double thold = startlk - 2;
	//cout <<  startlk << " "<< thold << "\n";
	double test =  log(startsegpi) - log(1-startsegpi);

	double min = -10.0;
	double max = 10.0;
	int nit = 0;
	double start = (test+max)/2;
	double tau = 0.001;
	int convhi = golden_section_segpi_ci(test, start, max, tau, thold,  &nit);
	double hi = segpi;
	cout << hi << " "<< llk() << " hi\n";
	start = (test+min)/2;
	nit = 0;
	int convlo = golden_section_segpi_ci(min, start, test, tau, thold, &nit);
	double lo = segpi;
	cout << lo << " "<< llk() << " lo\n";
	pair<int, int> conv = make_pair(convlo, convhi);
	pair<double, double> ci = make_pair(lo, hi);
	return make_pair(conv, ci);
}

*/

/*
pair< pair<int, int>, pair<double, double> > SNPs::get_cis_lambda(int which){
	pair< pair<int, int>, pair<double, double> > toreturn;

	double startlk = llk();
	double thold = startlk - 2;
	cout <<  startlk << " "<< thold << "\n";
	double min = -20.0;
	double max = 20.0;
	double test = lambdas[which];
	if (test > max) max = test+20.0;
	if (test < min) min = test-20.0;
	if (max < 0) max = 20.0;
	if (min > 0) min = -20.0;
	double start = (test+max)/2;
	double tau = 0.001;
	int convhi = golden_section_lambda_ci(test, start, max, tau, which, thold);
	double hi = lambdas[which];
	cout << hi << " "<< llk() << " hi\n";
	start = (test+min)/2;
	int convlo = golden_section_lambda_ci(min, start, test, tau, which, thold);
	double lo = lambdas[which];
	cout << lo << " "<< llk() << " lo\n";
	pair<int, int> conv = make_pair(convlo, convhi);
	pair<double, double> ci = make_pair(lo, hi);
	return make_pair(conv, ci);
}
*/

/*
pair< pair<int, int>, pair<double, double> > SNPs::get_cis_seglambda(int which){
	pair< pair<int, int>, pair<double, double> > toreturn;

	double startlk = llk();
	double thold = startlk - 2;
	cout <<  startlk << " "<< thold << "\n";
	double min = -20.0;
	double max = 20.0;
	double test = seglambdas[which];
	if (test > max) max = test+20.0;
	if (test < min) min = test-20.0;
	if (max < 0) max = 20.0;
	if (min > 0) min = -20.0;
	double start = (test+max)/2;
	double tau = 0.001;
	int convhi = golden_section_seglambda_ci(test, start, max, tau, which, thold);
	double hi = seglambdas[which];
	cout << hi << " "<< llk() << " hi\n";
	start = (test+min)/2;
	int convlo = golden_section_seglambda_ci(min, start, test, tau, which, thold);
	double lo = seglambdas[which];
	cout << lo << " "<< llk() << " lo\n";
	pair<int, int> conv = make_pair(convlo, convhi);
	pair<double, double> ci = make_pair(lo, hi);
	return make_pair(conv, ci);
}
*/

void SNPs_PW::print(){
	cout << "rs chr pos BF1 BF2";
	for (int i =0; i < nannot; i++) cout << " "<< annotnames[i];
	cout << "\n";
	for (vector<SNP_PW>::iterator it = d.begin(); it != d.end(); it++){
		cout << it->id << " "<< it->chr << " "<< it->pos << " "<< it->BF1 <<  " "<< it->BF2 << " "<< it->BF3;
		for (int i = 0; i < nannot; i++) cout << " "<< it->annot[i];
		cout << "\n";
	}
}


void SNPs_PW::print(string outfile, string outfile2){
	ogzstream out(outfile.c_str());
	ogzstream out2(outfile2.c_str());
	out << "id chr pos logBF_1 logBF_2 logBF_3 Z_"<< params->pheno1<<" V_"<<params->pheno1<<" Z_"<< params->pheno2 << " V_"<<params->pheno2<<" pi_1 pi_2 pi_3 PPA_1 PPA_2 PPA_3 chunk";
	out2 << "chunk NSNP chr st sp max_abs_Z_"<< params->pheno1<<" max_abs_Z_"<<params->pheno2<<" logBF_1 logBF_2 logBF_3 logBF_4 pi_1 pi_2 pi_3 pi_4 PPA_1 PPA_2 PPA_3 PPA_4";
	for (vector<string>::iterator it = annotnames.begin(); it != annotnames.end(); it++) out << " "<< *it;
	out << "\n";
	for (vector<string>::iterator it = segannotnames.begin(); it != segannotnames.end(); it++) out2 << " "<< *it;
	out2 << "\n";
	int segnum = 0;
	for (vector<pair<int, int> >::iterator it = segments.begin(); it != segments.end(); it++){
		int stindex = it->first;
		int spindex = it->second;
		out2 << segnum << " " << spindex-stindex << " "<< d[stindex].chr << " "<< d[stindex].pos << " "<< d[spindex-1].pos << " ";
		vector<double> segbfs = get_segbfs(segnum);
		if (params->finemap){
			double tmp = segbfs[0];
			for (int i = 0; i < 4; i++) segbfs[i] = segbfs[i] - tmp;
		}
		double maxZ1 = 0;
		double maxZ2 = 0;
		double sumlbf1 = -1000;
		double sumlbf2 = -1000;
		double sumlbf3 = -1000;
		for (int i = stindex; i < spindex; i++){
			double Z1 = fabs(d[i].Z1);
			double Z2 = fabs(d[i].Z2);
			sumlbf1 = sumlog(sumlbf1, d[i].BF1);
			sumlbf2 = sumlog(sumlbf2, d[i].BF2);
			sumlbf3 = sumlog(sumlbf3, d[i].BF3);
			if (Z1> maxZ1) maxZ1 = Z1;
			if (Z2> maxZ2) maxZ2 = Z2;
		}
		vector<double> lpostodds;
		vector<double> lpriorodds;
		for (int i = 0; i < 4; i++) lpriorodds.push_back(log(pi[i+1]/pi[0]));
		if (params->finemap){
			for (int i = 0; i < 4; i++) lpriorodds[i] = log(pi[i+1]/pi[1]);
		}
		for (int i = 0; i < 4; i++) lpostodds.push_back(segbfs[i]+lpriorodds[i]);
		vector<double> segPPA;
		for (int i = 0; i < 4; i++){
			double num = lpostodds[i];
			double denom = 0;
			for (int i = 0; i <4 ; i++) denom = sumlog(denom,  lpostodds[i]);
			if (params->finemap){
				denom = 0;
				for (int i = 1; i < 4; i++) denom = sumlog(denom, lpostodds[i]);
			}
			segPPA.push_back(exp(num - denom));
		}

		out2 << maxZ1<< " "<< maxZ2<< " "<< segbfs[0] <<" "<< segbfs[1]<< " "<< segbfs[2]<< " "<< segbfs[3] << " " << pi[1]<< " "<< pi[2]<< " "<< pi[3]<< " "<< pi[4] << " "<< segPPA[0]<< " "<< segPPA[1]<< " "<< segPPA[2]<< " "<< segPPA[3];
		for (int i = 0; i < nsegannot; i++) out2 << " "<< segannot[segnum][i];
		out2 << "\n";
		for (int i =stindex ; i < spindex; i++){
			out << d[i].id << " "<< d[i].chr << " "<< d[i].pos << " "<< d[i].BF1 << " "<< d[i].BF2 << " " << d[i].BF3<< " "<< d[i].Z1 <<  " " << d[i].V1 << " "<< d[i].Z2 <<  " " << d[i].V2<< " "<< snppri[i][0]<< " "<< snppri[i][1]<< " "<< snppri[i][2] << " " << exp(d[i].BF1- sumlbf1) << " "<< exp(d[i].BF2- sumlbf2) << " "<< exp(d[i].BF3- sumlbf3) << " "<< segnum;
			for (int j = 0; j < annotnames.size(); j++) out << " "<< d[i].annot[j];
			out << "\n";
		}
		segnum++;
	}
}


/*
double SNPs::cross10(bool penalize){
	//do 10-fold cross validation
	//
	// split segments into groups
	// L = 1/10* \sum_i L*(i)
	// where L*(i) is the likelihood of data in group i after optimizing model without it
	//
	double toreturn =0;
	vector< set<int> > split10 = make_cross10();
	vector<double> Lstar;
	for (vector<set<int> >::iterator it = split10.begin(); it != split10.end(); it++){
		GSL_xv_optim(*it, penalize);
		double tmpl = 0;
		for (set<int>::iterator it2 = it->begin(); it2 != it->end(); it2++) tmpl += llk(*it2);
		//cout << tmpl << "\n";
		Lstar.push_back(tmpl);
		toreturn += tmpl;
	}
	toreturn = toreturn /10.0;
	return toreturn;
}
*/

/*
vector<set<int> > SNPs::make_cross10(){
	vector<set<int> > toreturn;
	int nper = floor((double) segments.size() / 10.0);
	for (int i = 0; i < 10; i++){
		set<int> tmp;
		for (int j = i*nper; j < i*nper+nper; j++) tmp.insert(j);
		toreturn.push_back(tmp);
	}
	//for (vector<set<int> >::iterator it = toreturn.begin(); it != toreturn.end(); it++){
	//	for (set<int>::iterator it2 = it->begin(); it2 != it->end(); it2++){
	//		cout << *it2 << " ";
	//	}
	//	cout << "\n";
	//}

	return toreturn;

}

*/

void SNPs_PW::make_segments(int size){
	segments.clear();
	int counter = 0;
	for (vector<pair<int, int> >::iterator it = chrsegments.begin(); it != chrsegments.end(); it++){
		int starti = it->first;
		int endi = it->second;
		int length = endi-starti;
		int bestmod = length % size;
		int bestsize = size;
		for (int i = size -20; i < size+20; i++){
			if (i < 1) continue;
			int test = length % i;
			if (test < bestmod){
				bestmod = test;
				bestsize = i;
			}
		}
		int nseg = length/bestsize;
		for (int i = 0; i < nseg; i++){
			int sstart = starti+i*bestsize;
			int send = starti+i*bestsize+bestsize;
			if (i > (nseg-2))	send = endi;
			segments.push_back(make_pair(sstart, send));
			for (int i = sstart ; i < send ; i++) d[i].chunknumber = counter;
			counter++;
		}
	}
}


void SNPs_PW::make_segments_finemap(){
	//cout << "here\n"; cout.flush();
	segments.clear();
	int sstart = 0;
	int wseg = d[0].chunknumber;
	for (int i = 1; i < d.size(); i++){
		int testseg = d[i].chunknumber;
		if (testseg < wseg){
			cerr<< "ERROR: segment number "<< testseg << " occurs after "<< wseg << ". For numbered segments/finemapping, order the input file by SEGNUMBER.\n";
			exit(1);
		}
		if (testseg > wseg) {
			segments.push_back(make_pair(sstart, i));
			sstart = i;
			wseg = testseg;
		}
	}
	segments.push_back(make_pair(sstart, d.size()));
}

void SNPs_PW::make_chrsegments(){
	chrsegments.clear();
	chrnames.clear();
	int i = 0;
	int start = i;
	int startpos = d[i].pos;
	string startchr = d[i].chr;
	while (i < d.size()){
		int tmppos = d[i].pos;
		string tmpchr = d[i].chr;
		if (tmpchr != startchr){
			int end = i;
			chrnames.push_back(startchr);
			chrsegments.push_back(make_pair(start, end));
			start = i;
			startpos = d[i].pos;
			startchr = d[i].chr;
		}
		i++;
	}
	int end = i;
	chrnames.push_back(startchr);
	chrsegments.push_back(make_pair(start, end));
}

void SNPs_PW::print_segments(){
	for (vector<pair<int, int> >::iterator it = segments.begin(); it != segments.end(); it++){
		cout << it->first << " "<< it->second << "\n";
	}

}

void SNPs_PW::print_chrsegments(){
	for (vector<pair<int, int> >::iterator it = chrsegments.begin(); it != chrsegments.end(); it++){
		cout << it->first << " "<< it->second << "\n";
	}

}


void SNPs_PW::set_priors(){

	set_segpriors(); // right now only segment priors, constant across segments
	for (int i = 0; i < segments.size(); i++) set_priors(i);
}


/*
void SNPs::set_priors_cond(){
	set_segpriors();
	for (int i = 0; i < segments.size(); i++) set_priors_cond(i);
}
*/


void SNPs_PW::set_segpriors(){
	// 5 models, prior on model i = exp(a_i)/ sum_j(exp)a_j))
	assert (alpha.size()==5);
	vector<double> segp;
	double s = 0;
	int sti = 0;
	if (params->finemap) sti = 1;
	//cout << sti << " sti\n"; cout.flush();
	for (int i = sti; i < alpha.size(); i++) {
		segp.push_back(exp(alpha[i]));
		s+= exp(alpha[i]);
	}
	for (int i = sti; i < alpha.size(); i++) pi[i] = segp[i-sti]/s;

	//for (int i = 0; i < pi.size(); i++){
	//		cout << pi[i] << " ok\n";
	//}
	// priors now constant across segments
	/*
	for (int i = 0; i < segments.size(); i++){

		double logitprior = log(segpi) - log(1-segpi);
		for (int j = 0; j < nsegannot; j++){
			if (segannot[i][j]) logitprior += seglambdas[j];
		}
		double prior = 1.0/  (1.0 + exp(-logitprior));

		segpriors[i] = prior;
	}
	*/
}




void SNPs_PW::set_priors(int which){
	pair<int, int> seg = segments[which];
	int st = seg.first;
	int sp = seg.second;

	double sumxs = d[st].get_x(lambdas);
	vector<double> tmps;
	tmps.push_back(sumxs);

	for (int i = st+1; i < sp ;i++) {
		double tmpx = d[i].get_x(lambdas);
		tmps.push_back(tmpx);
		sumxs = sumlog(sumxs, tmpx);
	}

	for (int i = st; i <  sp ; i++) {
		//doing this in log space
		snppri.at(i).at(0)= tmps.at(i-st) - sumxs;
		snppri.at(i).at(1)= tmps.at(i-st) - sumxs;
		snppri.at(i).at(2)= tmps.at(i-st) - sumxs;
		if (!isfinite(snppri[i][0]) || !isfinite(snppri[i][1]) || !isfinite(snppri[i][2]) ){
			cerr << "ERROR: prior for SNP "<< i << " is " << snppri[i][0] <<" " << snppri[i][1]<< " "<< snppri[i][2] <<"\n";
			exit(1);
		}
	}
}


/*
void SNPs::set_priors_cond(int which){
	pair<int, int> seg = segments[which];
	int st = seg.first;
	int sp = seg.second;
	double sumxs = 0;
	for (int i = st; i < sp ;i++) {
		//cout << i << " "<< d[i].get_x(lambdas) << "\n";
		//double tmpx = exp(d[i].get_x(lambdas));
		double tmpx = d[i].get_x_cond(lambdas, condlambda);
		snppri[i] = tmpx;
		//sumxs += tmpx;
		sumxs = sumlog(sumxs, tmpx);
	}
	for (int i = st; i <  sp ; i++) {
		//snppri[i] = snppri[i]/sumxs;
		//doing this in log space
		snppri[i] = snppri[i] - sumxs;
		if (!isfinite(snppri[i])){
			cerr << "ERROR: prior for SNP "<< i << " is " << snppri[i] << "\n";
			exit(1);
		}
		//if (snppri[i] < DBL_MIN) snppri[i] = DBL_MIN;
		//cout << i << " "<< snppri[i] << "\n";
	}
}

*/
double SNPs_PW::llk(){
	if (!precomputed) seg_toadd.clear();
	double toreturn = -1000;
	for (int i = 0; i < segments.size(); i++) toreturn += llk(i);
	data_llk = toreturn;
	if (!precomputed) precomputed = true;

	// add on prior on proportions (for MAP estimator)
	// this is primarily for numerical stability
	double pi_lk = lnvecdens(alpha, params->alpha_prior);
	toreturn = toreturn+pi_lk;
	return toreturn;
}

/*
double SNPs::llk_xv(set<int> skip, bool penalize){
	double toreturn = 0;
	for (int i = 0; i < segments.size(); i++) {
		if (skip.find(i) != skip.end()) continue;
		toreturn += llk(i);
	}
	if (penalize){
		for (vector<double>::iterator it = lambdas.begin(); it != lambdas.end(); it++) toreturn -= params->ridge_penalty * *it * *it;
		for (vector<double>::iterator it = seglambdas.begin(); it != seglambdas.end(); it++) toreturn -= params->ridge_penalty * *it * *it;

	}//data_llk = toreturn;
	return toreturn;
}
*/
/*
double SNPs::llk_ridge(){
	double toreturn = 0;
	for (int i = 0; i < segments.size(); i++) toreturn += llk(i);
	for (vector<double>::iterator it = lambdas.begin(); it != lambdas.end(); it++) toreturn -= params->ridge_penalty * *it * *it;
	for (vector<double>::iterator it = seglambdas.begin(); it != seglambdas.end(); it++) toreturn -= params->ridge_penalty * *it * *it;
	data_llk = toreturn;
	return toreturn;
}
*/
void SNPs_PW::get_all_condZ(){
	string outfile = params->outstem +".condZ.gz";
	ogzstream out(outfile.c_str());
	out << "SNP Z_1 V_1 Z_2 V_2 BF_1 BF_2 conditional_on R CZ_1 CV_1 CZ_2 CV_2 CBF_2\n";
	for (int i = 0; i < segments.size(); i++) get_condZ(i, out);
}
void SNPs_PW::get_condZ(int which, ogzstream & out){
	// conditional Z scores for all pairs of SNPs in a segment
	//
	if (!params->overlap){
		cerr<< "ERROR: calling conditional Z scores but no LD matrix\n";
		exit(1);
	}
	else{
		pair<int, int> seg = segments[which];
		int st = seg.first;
		int sp = seg.second;

		vector<int> poss;
		for (int i = st; i < sp ; i++) poss.push_back(d[i].pos);
		LDmatrix ld(params->ldfile, d[st].chr, poss, params->Nhap);
		for (int i = st; i < sp ; i++){
			for (int j = i+1; j < sp ; j++){
				double tmpVi = ld.get_ld(d[i].pos, d[i].pos);
				double tmpVj = ld.get_ld(d[j].pos, d[j].pos);


				double  VarR_i = tmpVi/tmpVj; //this is pi (1-pi) / (pj(1-pj))
				double VarR_j = tmpVj/tmpVi; //this is pj (1-pj)/ (pi (1-pi))

				pair<double, double> Ri = ld.get_R(d[i].pos, d[j].pos); //this is D/Var( V_i )
				pair<double, double> Rj = ld.get_R(d[j].pos, d[i].pos); //this is D/Var( V_j );

				//conditional effects of SNP j
				pair<pair<double, double>, pair<double, double> > ef1 = d[j].condZ(&d[i], Rj, VarR_j);
				out << d[j].id << " "<< d[j].Z1 << " " << d[j].V1 << " " << d[j].Z2 << " "<< d[j].V2 << " "<< d[j].BF1<< " "<< d[j].BF2 << " "<< d[i].id << " " << Rj.first << " " << ef1.first.first/sqrt(ef1.first.second) << " "<< ef1.first.second << " "<<ef1.second.first/sqrt(ef1.second.second) << " "<< ef1.second.second << " "<< d[j].BF2_C(&d[i], params->cor, Rj, VarR_j) <<"\n";

				//conditional effects of SNP i
				pair<pair<double, double>, pair<double, double> > ef2 = d[i].condZ(&d[j], Ri, VarR_i);
				out << d[i].id << " "<< d[i].Z1 << " " << d[i].V1 << " " << d[i].Z2 << " "<< d[i].V2 << " "<< d[i].BF1 << " "<< d[i].BF2 << " "<< d[j].id << " "<< Ri.first << " " <<ef2.first.first/sqrt(ef2.first.second) << " "<< ef2.first.second << " "<<ef2.second.first/sqrt(ef2.second.second) << " "<< ef2.second.second << " "<< d[i].BF2_C(&d[j], params->cor, Ri, VarR_i)<< "\n";


			}
		}
	}
}


double SNPs_PW::llk(int which){
	//
	// M0: no SNPs affect any trait
	// M1: 1 SNP affects trait 1
	// M2: 1 SNP affects trait 2
	// M3: 1 SNP affects both trait
	// M4: 2 SNPs, 1 affects each trait
	//
	// logLK ~ log(P(M0)+ P(M1)*LK(M1)/LK(M0)  + P(M2)*LK(M2)/LK(M0) + P(M3)*LK(M3)/LK(M0)+ P(M4)*LK(M4)/LK(M0))
	//
	// LK(M1)/LK(M0) = sum(pi * BF1)
	// ...
	//
	// only complication is when this is a single cohort, then need to compute M4 worrying about LD
	//
	double toreturn = 0;
	if (precomputed){
		double m1 = seg_toadd.at(which).at(0);
		double m2 = seg_toadd.at(which).at(1);
		double m3 = seg_toadd.at(which).at(2);
		double m4 = seg_toadd.at(which).at(3);


		m1 = m1+log(pi[1]);
		m2 = m2+log(pi[2]);
		m3 = m3+log(pi[3]);
		m4 = m4+log(pi[4]);

		if (params->finemap){
			double tmp = 1+ exp(m1-m3)+exp(m2-m3)+exp(m4-m3);
			toreturn = m3+log(tmp);
			return toreturn;
		}

		double m0 = log(pi[0]);
		double tmp = 1+exp(m0-m3)+ exp(m1-m3)+exp(m2-m3)+exp(m4-m3);
		toreturn = m3+log(tmp);
		//cout << m0 << " "<< m1 << " "<< m2 << " "<< m3 << " "<< m4 << " " << toreturn << "\n";
		return toreturn;
	}
	pair<int, int> seg = segments[which];
	int st = seg.first;
	int sp = seg.second;

	//initialize to ~0 (in log space)
	double m1 = -1000;
	double m2 = -1000;
	double m3 = -1000;
	double m4 = -1000;


	//int counter = 0;
	// if these are non-overlapping cohorts
	if (!params->overlap){
		for (int i = st; i < sp ; i++){

			//term 1: one associated SNP for pheno 1
			double tmp2add1 = snppri.at(i).at(0)+ d[i].BF1;
			m1 = sumlog(m1, tmp2add1);


			//term 2: one associated SNP for pheno 2
			double tmp2add2 = snppri.at(i).at(1)+ d[i].BF2;
			m2 = sumlog(m2, tmp2add2);

			//term 3: one associated SNP, both phenos
			double tmp2add3 = snppri.at(i).at(2)+ d[i].BF3;
			m3 = sumlog(m3, tmp2add3);

			//term 4: two associated SNPs, both phenos
			double tmp2add4 = -10000;
			for (int j = i+1; j < sp ; j++){


				double tmp2_4 = snppri.at(i).at(0)+snppri.at(j).at(1)+d[i].BF1+d[j].BF2;
				double tmp2_42 = snppri.at(i).at(1)+snppri.at(j).at(0)+d[j].BF1+d[i].BF2;
				//cout << d[i].id << " "<< d[j].id << " "<< tmp2_4 << " "<< d[i].BF1+d[j].BF2 << " "<< tmp2_42 <<  " "<< d[j].BF1+d[i].BF2 << "\n";
				tmp2add4 = sumlog(tmp2add4, tmp2_4);

				tmp2add4 = sumlog(tmp2add4, tmp2_42);

			}
			m4 = sumlog(m4, tmp2add4);
		}
	}
	else{
		vector<int> poss;
		for (int i = st; i < sp ; i++) {
			//cout << d[i].pos << "\n";
			poss.push_back(d[i].pos);
		}
		//cout << "readld\n"; cout.flush();
		LDmatrix ld(params->ldfile, d[st].chr, poss, params->Nhap);
		//cout << "done\n"; cout.flush();
		for (int i = st; i < sp ; i++){

			//term 1: one associated SNP for pheno 1
			double tmp2add1 = snppri.at(i).at(0)+ d[i].BF1;
			m1 = sumlog(m1, tmp2add1);


			//term 2: one associated SNP for pheno 2
			double tmp2add2 = snppri.at(i).at(1)+ d[i].BF2;
			m2 = sumlog(m2, tmp2add2);

			//term 3: one associated SNP, both phenos
			double tmp2add3 = snppri.at(i).at(2)+ d[i].BF3;
			m3 = sumlog(m3, tmp2add3);

			//term 4: two associated SNPs, both phenos
			double tmp2add4 = -10000;
			for (int j = i+1; j < sp ; j++){
				//double D = ld.get_ld(d[i].pos, d[j].pos);
				double tmpVi = ld.get_ld(d[i].pos, d[i].pos);
				double tmpVj = ld.get_ld(d[j].pos, d[j].pos);

				double  VarR_i = tmpVi/tmpVj; //this is pi (1-pi) / (pj(1-pj))
				double VarR_j = tmpVj/tmpVi; //this is pj (1-pj)/ (pi (1-pi))



				pair<double, double> Ri = ld.get_R(d[i].pos, d[j].pos); //this is D/Var( V_i )
				pair<double, double> Rj = ld.get_R(d[j].pos, d[i].pos); //this is D/Var( V_j );




				double BFj_ci, tmp2_4, BFi_cj, tmp2_42;
				if (!params->rev){
					//SNP i affects pheno 1, BF at SNP j for pheno 2 conditional on SNPi
					BFj_ci = d[j].BF2_C(&d[i], params->cor, Rj, VarR_j);
					tmp2_4 = snppri.at(i).at(0)+snppri.at(j).at(1)+d[i].BF1+BFj_ci;
					//if (BFj_ci > d[j].calc_logBF2(params->cor)+3){
					//	cout << d[j].id << " "<< d[i].id << " "<< Rj.first << " "<< BFj_ci << " "<< d[j].calc_logBF2(params->cor) << "\n";
					//}
					//SNP j affects pheno 1, BF at SNP i for pheno 2 conditional on SNPj
					BFi_cj = d[i].BF2_C(&d[j],params->cor, Ri, VarR_i);
					tmp2_42 = snppri.at(i).at(1)+snppri.at(j).at(0)+d[j].BF1+ BFi_cj;
				}
				else{
					BFj_ci = d[j].BF1_C(&d[i], params->cor, Rj, VarR_j);
					tmp2_4 = snppri.at(i).at(0)+snppri.at(j).at(1)+d[i].BF2+BFj_ci;

					BFi_cj = d[i].BF1_C(&d[j],params->cor, Ri, VarR_i);
					tmp2_42 = snppri.at(i).at(1)+snppri.at(j).at(0)+d[j].BF2+ BFi_cj;
				}


				tmp2add4 = sumlog(tmp2add4, tmp2_4);

				tmp2add4 = sumlog(tmp2add4, tmp2_42);

			}
			m4 = sumlog(m4, tmp2add4);
		}
	}
	vector<double> toadd;
	toadd.push_back(m1);toadd.push_back(m2);toadd.push_back(m3);toadd.push_back(m4);
	seg_toadd.push_back(toadd);


	//cout << m1 << " "<< m2 << " "<< m3 << " "<< m4<< "\n";
	m1 = m1+log(pi[1]);
	m2 = m2+log(pi[2]);
	m3 = m3+log(pi[3]);
	m4 = m4+log(pi[4]);

	if (params->finemap){
		double tmp = 1+ exp(m1-m3)+exp(m2-m3)+exp(m4-m3);
		cout << alpha[4] << " "<< pi[4] << " "<< m1 << " "<< m2 << " "<< m3 << " "<< m4 << " " << toreturn << " "<< which << "\n";
		toreturn = m3+log(tmp);
		return toreturn;
	}
	double m0 = log(pi[0]);
	double tmp = 1+exp(m0-m3)+ exp(m1-m3)+exp(m2-m3)+exp(m4-m3);
	toreturn = m3+log(tmp);
	cout << alpha[4] << " "<< pi[4] << " "<< m0 << " "<< m1 << " "<< m2 << " "<< m3 << " "<< m4 << " " << toreturn << " "<< which << "\n";
	return toreturn;
}

double SNPs_PW::sumlog(double logx, double logy){
        if (logx > logy) return logx + log(1 + exp(logy-logx));
        else return logy + log(1 + exp(logx-logy));
}

vector<double> SNPs_PW::get_segbfs(int i){
	vector<double> toreturn;
	if (!precomputed) llk();
	for (vector<double>::iterator it = seg_toadd[i].begin(); it != seg_toadd[i].end(); it++) toreturn.push_back(*it);
	return toreturn;
}

void SNPs_PW::MCMC(gsl_rng *r){
	llk();
	double spi1 = 0;
	double spi2 = 0;
	double spi3 = 0;
	double spi4 = 0;
	int ns = 0;
	string outMCMC = params->outstem+".MCMC";
	ofstream outr(outMCMC.c_str());
	outr << "i ";
	if (!params->finemap) outr << "pi_0 ";
	outr << "pi_1 pi_2 pi_3 pi_4 lk\n";
	for (int i = 0; i < params->burnin; i++) {
		MCMC_update(r);
		if (i % params->sampfreq ==0) {
			cout << i << " ";
			outr << "#" << i << " ";
			if (!params->finemap){
				cout << pi[0] << " ";
				outr << pi[0] << " ";
			}
			cout <<  pi[1]<< " "<<pi[2] << " "<< pi[3] << " "<< pi[4]<< " "<< data_llk << "\n";
			outr <<  pi[1]<< " "<<pi[2] << " "<< pi[3] << " "<< pi[4]<< " "<< data_llk << "\n";
		}
	}
	int nsamp = 0;
	int naccept = 0;
	for (int i = 0; i < params->nsamp; i++){
		naccept+= MCMC_update(r);
		nsamp++;
		if (i % params->sampfreq ==0) {
			cout << i << " ";
			outr << i << " ";
			if (!params->finemap){
				cout << pi[0] << " ";
				outr << pi[0] << " ";
			}
			cout <<  pi[1]<< " "<<pi[2] << " "<< pi[3] << " "<< pi[4]<< " "<< data_llk << "\n";
			outr <<  pi[1]<< " "<<pi[2] << " "<< pi[3] << " "<< pi[4]<< " "<< data_llk << "\n";
			spi1 += pi[1];
			spi2 += pi[2];
			spi3 += pi[3];
			spi4 += pi[4];
			ns++;
		}

	}
	outr << "# " << (double) naccept / (double) nsamp << " :acceptance probability\n";
	outr << "# " << spi1/ (double) ns << " : mean pi_1\n";
	outr << "# " << spi2/ (double) ns << " : mean pi_2\n";
	outr << "# " << spi3/ (double) ns << " : mean pi_3\n";
	outr << "# " << spi4/ (double) ns << " : mean pi_4\n";
	pi[1] = spi1/ (double) ns;
	pi[2] = spi2/ (double) ns;
	pi[3] = spi3/ (double) ns;
	pi[4] = spi4/ (double) ns;
}

int SNPs_PW::MCMC_update(gsl_rng *r){
	// return 1 if update accepted, 0 otw
	//save old values of alpha, pi
	vector<double> tmpalpha;
	vector<double> tmppi;
	double tmpllk = data_llk;
	for (int i = 0; i < alpha.size(); i++) {
		tmpalpha.push_back(alpha[i]);
		tmppi.push_back(pi[i]);
	}

	//propose new alpha, reset priors
	vector<double> newalpha = propose_alpha(r);
	for (int i = 0; i < 5; i++) alpha[i] = newalpha[i];
	set_priors();
	llk();

	//prior on initial alphas
	double tmppi_lk = lnvecdens(tmpalpha, params->alpha_prior);
	//prior on new alphas
	double pi_lk = lnvecdens(alpha, params->alpha_prior);


	//compute acceptance prob
	double num = data_llk+ pi_lk;
	double denom = tmpllk+ tmppi_lk;
	double diff = num - denom;
	double acceptance = exp(diff);
	if (acceptance > 1) return 1;

	//generate random unif
	double unif = gsl_rng_uniform(r);
	if (unif < acceptance) return 1;
	else{
		for (int i = 0; i < 5; i++) alpha[i] = tmpalpha[i];
		set_priors();
		llk();
	}
	return 0;


}

double SNPs_PW::lnvecdens(vector<double> current, vector<double> prior){
	vector<double> tmp;
	for (int i = 0; i < 5; i++) tmp.push_back( lndgauss(current[i]-prior[i], 3));
	double toreturn = tmp[0];
	for (int i = 1; i < 5; i++) toreturn+= tmp[i];
	return toreturn;
}

double SNPs_PW::lndgauss(double dif, double se){
        double toreturn = 0;
        toreturn += -log (se * sqrt(2.0*M_PI));
        toreturn += -(dif*dif) /(2*se*se);
        return toreturn;
}

double SNPs_PW::dirichlet_lndens(vector<double> alphas, vector<double> thetas){
	size_t K = alphas.size();
	double t[5];
	double a[5];
	for (int i = 0; i < 5; i++){
		t[i] = thetas[i];
		a[i] = alphas[i];
	}

	double toreturn = gsl_ran_dirichlet_lnpdf(K, a, t);
	return toreturn;
}

vector<double> SNPs_PW::propose_alpha(gsl_rng *r){
	//Propose new alphas, normally distributed around old alphas
	vector<double> toreturn;
	for (int i = 0; i < 5; i++){
		double tmp = gsl_ran_gaussian(r, params->MCMC_gauss_SD);
		toreturn.push_back(alpha[i]+ tmp);
	}
	return toreturn;
	//Dirichlet (had problems with hitting the edge of the simplex)
	/*
	vector<double> toreturn;
	double t[5];
	double a[5] = {alpha[0]* segments.size(), alpha[1]*segments.size(), alpha[2]*segments.size(), alpha[3]*segments.size(), alpha[4]*segments.size()};
	size_t K = 5;
	gsl_ran_dirichlet(r, K, a, t);
	for (int i =0 ; i < 5; i++) toreturn.push_back(t[i] * segments.size());
	return toreturn;
	*/
}
/*
void SNPs::set_post(){
	for (int i = 0; i < segments.size(); i++) set_post(i);
}

*/

/*
void SNPs::set_post(int which){
	pair<int, int> seg = segments[which];
	int st = seg.first;
	int sp = seg.second;
	double seglk = llk(which);
	for (int i = st; i < sp; i++){
		double num = log(segpi) + d[i].BF + snppri[i];
		//double num = log(segpi) + d[i].BF + log(snppri[i]);
		double lpost = num - seglk;
		snppost[i] = exp(lpost);
	}
}
*/
/*
void SNPs::print_segprobs(string outfile){
	ogzstream out(outfile.c_str());
	for (vector<pair<int, int> >::iterator it = segments.begin(); it != segments.end(); it++){
		int st = it->first;
		int sp = it->second;
		double total = 0;
		for (int i = st ; i < sp ; i++){
			total += snppost[i];
		}
		out << d[st].chr << " "<< d[st].pos << " "<< d[sp].pos << " " <<  total << "\n";
	}
}
*/

/*
void SNPs::optimize_segpi(){
	double start = log(segpi) - log(1-segpi);
	double min = -10.0;
	double max = 10.0;
	double tau = 0.001;
	golden_section_segpi(min, start, max, tau);
	cout << segpi << "\n";
}
*/

/*
void SNPs::optimize_condlambda(){
	double start = 0;
	double min = -20.0;
	double max = 20.0;
	double tau = 0.001;
	golden_section_condlambda(min, start, max, tau);
	cout << condlambda << "\n";
}
*/


void SNPs_PW::GSL_optim(){
	//cout << "here2\n"; cout.flush();
	if (params->finemap){
		GSL_optim_fine();
		return;
	}
	//cout << "here3\n"; cout.flush();
	int nparam = alpha.size()-1;
	size_t iter = 0;
	double size;
    int status;
    const gsl_multimin_fminimizer_type *T =
    		gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s;
    gsl_vector *x;
    gsl_vector *ss;
    gsl_multimin_function lm;
    lm.n = nparam;
    lm.f = &GSL_llk;
    struct GSL_params p;
    p.d = this;
    lm.params = &p;
    //
    // initialize parameters
    //
    x = gsl_vector_alloc (nparam);
    for (int i = 0; i < nparam; i++)   gsl_vector_set(x, i, alpha[i+1]);

    // set initial step sizes to 1
    ss = gsl_vector_alloc(nparam);
    gsl_vector_set_all(ss, 1.0);
    s = gsl_multimin_fminimizer_alloc (T, nparam);

    gsl_multimin_fminimizer_set (s, &lm, x, ss);
    do
     {
             iter++;
             status = gsl_multimin_fminimizer_iterate (s);

             if (status){
                     printf ("error: %s\n", gsl_strerror (status));
                     break;
             }
             size = gsl_multimin_fminimizer_size(s);
             status = gsl_multimin_test_size (size, 0.001);
             //cout << iter << " "<< iter %10 << "\n";
             if (iter % 20 < 1 || iter < 20){
            	 cout <<"iteration: "<< iter << " "<< alpha[0]<< " "<< alpha[1]<< " "<< alpha[2]<< " "<< alpha[3]<< " "<< alpha[4];
            	 cout << " "<< s->fval << " "<< size <<  "\n";
             }

     }
     while (status == GSL_CONTINUE && iter <5000);
     if (iter > 4999) {
             cerr << "WARNING: failed to converge. Maybe try -mcmc to check parameter estimates.\n";
             //exit(1);
     }
     //segpi = 1.0 /  (1.0 + exp (- gsl_vector_get(s->x, 0)));
     for (int i = 0; i <  nparam; i++) {
    	 alpha[i+1] = gsl_vector_get(s->x, i);
    	 //cout << "finished alpha "<< i+1 << " "<< alpha[i+1] << "\n";
     }


     gsl_multimin_fminimizer_free (s);
     gsl_vector_free (x);
     gsl_vector_free(ss);
}


void SNPs_PW::GSL_optim_fine(){

	int nparam = alpha.size()-2;
	size_t iter = 0;
	double size;
    int status;
    const gsl_multimin_fminimizer_type *T =
    		gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s;
    gsl_vector *x;
    gsl_vector *ss;
    gsl_multimin_function lm;
    lm.n = nparam;
    lm.f = &GSL_llk_fine;
    struct GSL_params p;
    p.d = this;
    lm.params = &p;
    //
    // initialize parameters
    //

    x = gsl_vector_alloc (nparam);
    for (int i = 0; i < nparam; i++)   gsl_vector_set(x, i, alpha[i+2]);

    // set initial step sizes to 1
    ss = gsl_vector_alloc(nparam);
    gsl_vector_set_all(ss, 1.0);
    s = gsl_multimin_fminimizer_alloc (T, nparam);
   // cout << "here3\n"; cout.flush();
    gsl_multimin_fminimizer_set (s, &lm, x, ss);

    do
     {
             iter++;
             status = gsl_multimin_fminimizer_iterate (s);

             if (status){
                     printf ("error: %s\n", gsl_strerror (status));
                     break;
             }
             size = gsl_multimin_fminimizer_size(s);
             status = gsl_multimin_test_size (size, 0.001);
             //cout << iter << " "<< iter %10 << "\n";
             if (iter % 20 < 1 || iter < 20){
            	 cout <<"iteration: "<< iter << " "<< alpha[0]<< " "<< alpha[1]<< " "<< alpha[2]<< " "<< alpha[3]<< " "<< alpha[4];
            	 cout << " "<< s->fval << " "<< size <<  "\n";
             }

     }
     while (status == GSL_CONTINUE && iter <5000);
     if (iter > 4999) {
             cerr << "WARNING: failed to converge. Maybe try -mcmc to check parameter estimates.\n";
             //exit(1);
     }
     //segpi = 1.0 /  (1.0 + exp (- gsl_vector_get(s->x, 0)));
     for (int i = 0; i <  nparam; i++) {
    	 alpha[i+2] = gsl_vector_get(s->x, i);
    	 //cout << "finished alpha "<< i+1 << " "<< alpha[i+1] << "\n";
     }


     gsl_multimin_fminimizer_free (s);
     gsl_vector_free (x);
     gsl_vector_free(ss);
}

/*
void SNPs::GSL_xv_optim(set<int> toskip, bool penalize){
	int nparam = nannot+nsegannot+1;
	size_t iter = 0;
	double size;
    int status;
    const gsl_multimin_fminimizer_type *T =
    		gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s;
    gsl_vector *x;
    gsl_vector *ss;
    gsl_multimin_function lm;
    lm.n = nparam;
    lm.f = &GSL_llk_xv;
    struct GSL_params p;
    p.d = this;
    p.toskip = toskip;
    p.penalize = penalize;
    lm.params = &p;
    //cout << llk()<< "\n"; cout.flush();
    //
    // initialize parameters
    //
    x = gsl_vector_alloc (nparam);
    gsl_vector_set(x, 0, log(segpi) - log(1-segpi));
    for (int i = 0; i < nparam-1; i++)   gsl_vector_set(x, i+1, 1);

    // set initial step sizes to 1
    ss = gsl_vector_alloc(nparam);
    gsl_vector_set_all(ss, 1.0);
    s = gsl_multimin_fminimizer_alloc (T, nparam);

    gsl_multimin_fminimizer_set (s, &lm, x, ss);
    do
     {
             iter++;
             status = gsl_multimin_fminimizer_iterate (s);

             if (status){
                     printf ("error: %s\n", gsl_strerror (status));
                     break;
             }
             size = gsl_multimin_fminimizer_size(s);
             status = gsl_multimin_test_size (size, 0.001);
             //cout << iter << " "<< iter %10 << "\n";
             if (iter % 20 < 1 || iter < 2){
            	 cout <<"iteration: "<< iter << " "<< segpi;
            	 for (int i = 0; i < nsegannot; i++) cout <<  " "<< seglambdas[i];
            	 for (int i = 0; i < nannot; i++) cout << " "<< lambdas[i];
            	 cout << " "<< s->fval << " "<< size <<  "\n";
             }

     }
     while (status == GSL_CONTINUE && iter <5000);
     if (iter > 4999) {
             cerr << "WARNING: failed to converge\n";
             //exit(1);
     }
     segpi = 1.0 /  (1.0 + exp (- gsl_vector_get(s->x, 0)));
     for (int i = 0; i < nsegannot; i++) seglambdas[i] = gsl_vector_get(s->x, i+1);
     for (int i = 0; i < nannot; i++) lambdas[i] = gsl_vector_get(s->x, i+nsegannot+1);


     gsl_multimin_fminimizer_free (s);
     gsl_vector_free (x);
     gsl_vector_free(ss);
}




void SNPs::GSL_optim_fine(){
	int nparam = nannot;
	if (nparam < 1) return;
	size_t iter = 0;
	double size;
    int status;
    const gsl_multimin_fminimizer_type *T =
    		gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s;
    gsl_vector *x;
    gsl_vector *ss;
    gsl_multimin_function lm;
    lm.n = nparam;
    lm.f = &GSL_llk_fine;
    struct GSL_params p;
    p.d = this;
    lm.params = &p;
    //cout << llk()<< "\n"; cout.flush();
    //
    // initialize parameters
    //

    x = gsl_vector_alloc (nparam);
    for (int i = 0; i < nannot; i++)   gsl_vector_set(x, i, 1);

    // set initial step sizes to 1
    ss = gsl_vector_alloc(nparam);
    gsl_vector_set_all(ss, 1.0);

    s = gsl_multimin_fminimizer_alloc (T, nparam);
    //cout << "here o1\n"; cout.flush();
    gsl_multimin_fminimizer_set (s, &lm, x, ss);
    //cout << "here2\n"; cout.flush();
    do
     {
             iter++;
             status = gsl_multimin_fminimizer_iterate (s);

             if (status){
                     printf ("error: %s\n", gsl_strerror (status));
                     break;
             }
             size = gsl_multimin_fminimizer_size(s);
             status = gsl_multimin_test_size (size, 0.001);
             //cout << iter << " "<< iter %10 << "\n";
             if (iter % 20 < 1 || iter < 2){
            	 cout <<"iteration: "<< iter;
            	 for (int i = 0; i < nannot; i++) cout << " "<< lambdas[i];
            	 cout << " "<< s->fval << " "<< size <<  "\n";
             }

     }
     while (status == GSL_CONTINUE && iter <5000);
     if (iter > 4999) {
             cerr << "WARNING: failed to converge\n";
             //exit(1);
     }
     for (int i = 0; i < nannot; i++) lambdas[i] = gsl_vector_get(s->x, i);


     gsl_multimin_fminimizer_free (s);
     gsl_vector_free (x);
     gsl_vector_free(ss);
}



void SNPs::GSL_optim_ridge_fine(){
	int nparam = nannot;
	if (nparam < 1) return;
	size_t iter = 0;
	double size;
    int status;
    const gsl_multimin_fminimizer_type *T =
    		gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s;
    gsl_vector *x;
    gsl_vector *ss;
    gsl_multimin_function lm;
    lm.n = nparam;
    lm.f = &GSL_llk_ridge_fine;
    struct GSL_params p;
    p.d = this;
    lm.params = &p;
    //cout << llk()<< "\n"; cout.flush();
    //
    // initialize parameters
    //

    x = gsl_vector_alloc (nparam);
    for (int i = 0; i < nannot; i++)   gsl_vector_set(x, i, 1);

    // set initial step sizes to 1
    ss = gsl_vector_alloc(nparam);
    gsl_vector_set_all(ss, 1.0);

    s = gsl_multimin_fminimizer_alloc (T, nparam);
    //cout << "here o1\n"; cout.flush();
    gsl_multimin_fminimizer_set (s, &lm, x, ss);
    //cout << "here2\n"; cout.flush();
    do
     {
             iter++;
             status = gsl_multimin_fminimizer_iterate (s);

             if (status){
                     printf ("error: %s\n", gsl_strerror (status));
                     break;
             }
             size = gsl_multimin_fminimizer_size(s);
             status = gsl_multimin_test_size (size, 0.001);
             //cout << iter << " "<< iter %10 << "\n";
             if (iter % 20 < 1 || iter < 2){
            	 cout <<"iteration: "<< iter;
            	 for (int i = 0; i < nannot; i++) cout << " "<< lambdas[i];
            	 cout << " "<< s->fval << " "<< size <<  "\n";
             }

     }
     while (status == GSL_CONTINUE && iter <5000);
     if (iter > 4999) {
             cerr << "WARNING: failed to converge\n";
             //exit(1);
     }
     for (int i = 0; i < nannot; i++) lambdas[i] = gsl_vector_get(s->x, i);


     gsl_multimin_fminimizer_free (s);
     gsl_vector_free (x);
     gsl_vector_free(ss);
}

void SNPs::GSL_optim_ridge(){
	if (params->finemap){
		GSL_optim_ridge_fine();
		return;
	}
	int nparam = nannot+nsegannot+1;
	size_t iter = 0;
	double size;
    int status;
    const gsl_multimin_fminimizer_type *T =
    		gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s;
    gsl_vector *x;
    gsl_vector *ss;
    gsl_multimin_function lm;
    lm.n = nparam;
    lm.f = &GSL_llk_ridge;
    struct GSL_params p;
    p.d = this;
    lm.params = &p;
    //cout << llk()<< "\n"; cout.flush();
    //
    // initialize parameters
    //
    x = gsl_vector_alloc (nparam);
    gsl_vector_set(x, 0, log(segpi) - log(1-segpi));
    for (int i = 0; i < nparam-1; i++)   gsl_vector_set(x, i+1, 1);

    // set initial step sizes to 1
    ss = gsl_vector_alloc(nparam);
    gsl_vector_set_all(ss, 1.0);
    s = gsl_multimin_fminimizer_alloc (T, nparam);

    gsl_multimin_fminimizer_set (s, &lm, x, ss);
    do
     {
             iter++;
             status = gsl_multimin_fminimizer_iterate (s);

             if (status){
                     printf ("error: %s\n", gsl_strerror (status));
                     break;
             }
             size = gsl_multimin_fminimizer_size(s);
             status = gsl_multimin_test_size (size, 0.001);
             //cout << iter << " "<< iter %10 << "\n";
             if (iter % 20 < 1 || iter < 2){
            	 cout <<"iteration: "<< iter << " "<< segpi;
            	 for (int i = 0; i < nsegannot; i++) cout <<  " "<< seglambdas[i];
            	 for (int i = 0; i < nannot; i++) cout << " "<< lambdas[i];
            	 cout << " "<< s->fval << " "<< size <<  "\n";
             }

     }
     while (status == GSL_CONTINUE && iter <5000);
     if (iter > 4999) {
             cerr << "WARNING: failed to converge\n";
             //exit(1);
     }
     segpi = 1.0 /  (1.0 + exp (- gsl_vector_get(s->x, 0)));
     for (int i = 0; i < nsegannot; i++) seglambdas[i] = gsl_vector_get(s->x, i+1);
     for (int i = 0; i < nannot; i++) lambdas[i] = gsl_vector_get(s->x, i+nsegannot+1);


     gsl_multimin_fminimizer_free (s);
     gsl_vector_free (x);
     gsl_vector_free(ss);
}

int SNPs::golden_section_segpi(double min, double guess, double max, double tau){
        double x;

        if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
        else x = guess - resphi *(guess-min);
        if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
                double new_segpi = (min+max)/2;
                segpi =  1.0  / ( 1.0 + exp(-new_segpi));
                data_llk = llk();
                return 0;
        }

        segpi = 1.0  / ( 1.0 + exp(-x));
        double f_x = -llk();

        segpi =  1.0  / ( 1.0 + exp(-guess));
        double f_guess = -llk();
        cout << x << " " <<  guess << " "<< segpi << " "<< f_x << " "<< f_guess <<  "\n";
        if (f_x < f_guess){
                if ( (max-guess) > (guess-min) )        return golden_section_segpi(guess, x, max, tau);
                else return golden_section_segpi(min, x, guess, tau);
        }
        else{
                if ( (max - guess) > (guess - min)  ) return golden_section_segpi(min, guess, x, tau);
                else return golden_section_segpi(x, guess, max, tau);
        }
}


int SNPs::golden_section_condlambda(double min, double guess, double max, double tau){
        double x;
        if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
        else x = guess - resphi *(guess-min);
        if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
                double new_segpi = (min+max)/2;
                condlambda =  new_segpi;
                set_priors_cond();
                data_llk = llk();
                return 0;
        }

        condlambda= x;
        set_priors_cond();
        double f_x = -llk();

        condlambda = guess;
        set_priors_cond();
        double f_guess = -llk();

        cout << x << " " <<  guess << " "<< f_x << " "<< f_guess  << " "<< max << " "<< min << "\n";
        if (f_x < f_guess){
                if ( (max-guess) > (guess-min) )        return golden_section_condlambda(guess, x, max, tau);
                else return golden_section_condlambda(min, x, guess, tau);
        }
        else{
                if ( (max - guess) > (guess - min)  ) return golden_section_condlambda(min, guess, x, tau);
                else return golden_section_condlambda(x, guess, max, tau);
        }
}

int SNPs::golden_section_condlambda_ci(double min, double guess, double max, double tau, double target){
        double x;

        if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
        else x = guess - resphi *(guess-min);
        if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
                double new_segpi = (min+max)/2;
                condlambda =  new_segpi;
                set_priors_cond();
                data_llk = llk();
                double tmpdiff = data_llk- target;
                double d2 = tmpdiff*tmpdiff;
                if (d2 > tau) return 1;
                else return 0;
        }

        condlambda = x;
        set_priors_cond();
        double f_x = llk()-target;
        f_x = f_x*f_x;
       // double x_llk = llk();

        condlambda = guess;
        set_priors_cond();
        double f_guess = llk()-target;
        f_guess = f_guess*f_guess;
       // double guess_llk = llk();
        cout << x << " " <<  guess << " "<< f_x << " "<< f_guess  << " "<< max << " "<< min << "\n";
        if (f_x < f_guess){
                if ( (max-guess) > (guess-min) )        return golden_section_condlambda_ci(guess, x, max, tau, target);
                else return golden_section_condlambda_ci(min, x, guess, tau, target);
        }
        else{
                if ( (max - guess) > (guess - min)  ) return golden_section_condlambda_ci(min, guess, x, tau, target);
                else return golden_section_condlambda_ci(x, guess, max, tau, target);
        }
}

*/

int SNPs_PW::golden_section_alpha_ci(double min, double guess, double max, double tau, int which, double target){
        double x;

        if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
        else x = guess - resphi *(guess-min);
        if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
                double new_segpi = (min+max)/2;
                alpha[which] =  new_segpi;
                set_priors();
                data_llk = llk();
                double tmpdiff = data_llk- target;
                double d2 = tmpdiff*tmpdiff;
                if (d2 > tau) return 1;
                else return 0;
        }

        alpha[which] = x;
        set_priors();
        double l_x = llk();
        double f_x = l_x-target;
        f_x = f_x*f_x;
       // double x_llk = llk();

        alpha[which] = guess;
        set_priors();
        double l_guess = llk();
        double f_guess = l_guess-target;
        f_guess = f_guess*f_guess;
       // double guess_llk = llk();
        cout << x << " " <<  guess << " "<< f_x << " "<< f_guess  << " "<< max << " "<< min << " "<< l_x << " "<< l_guess<< " "<< target <<"\n";
        if (f_x < f_guess){
                if ( (max-guess) > (guess-min) )        return golden_section_alpha_ci(guess, x, max, tau, which, target);
                else return golden_section_alpha_ci(min, x, guess, tau, which, target);
        }
        else{
                if ( (max - guess) > (guess - min)  ) return golden_section_alpha_ci(min, guess, x, tau, which, target);
                else return golden_section_alpha_ci(x, guess, max, tau, which, target);
        }
}

/*
int SNPs::golden_section_seglambda_ci(double min, double guess, double max, double tau, int which, double target){
        double x;

        if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
        else x = guess - resphi *(guess-min);
        if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
                double new_segpi = (min+max)/2;
                seglambdas[which] =  new_segpi;
                set_priors();
                data_llk = llk();
                double tmpdiff = data_llk- target;
                double d2 = tmpdiff*tmpdiff;
                if (d2 > tau) return 1;
                else return 0;
        }

        seglambdas[which] = x;
        set_priors();
        double f_x = llk()-target;
        f_x = f_x*f_x;
        //double x_llk = llk();

        seglambdas[which] = guess;
        set_priors();
        double f_guess = llk()-target;
        f_guess = f_guess*f_guess;
        //double guess_llk = llk();
        cout << x << " " <<  guess << " "<< f_x << " "<< f_guess << " "<< max << " "<< min << "\n";
        if (f_x < f_guess){
                if ( (max-guess) > (guess-min) )        return golden_section_seglambda_ci(guess, x, max, tau, which, target);
                else return golden_section_seglambda_ci(min, x, guess, tau, which, target);
        }
        else{
                if ( (max - guess) > (guess - min)  ) return golden_section_seglambda_ci(min, guess, x, tau, which, target);
                else return golden_section_seglambda_ci(x, guess, max, tau, which, target);
        }
}



int SNPs::golden_section_segpi_ci(double min, double guess, double max, double tau, double target, int * nit){
        double x;
        *nit = *nit +1;
        if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
        else x = guess - resphi *(guess-min);
        if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
                double new_segpi = (min+max)/2;
                segpi =  1.0  / ( 1.0 + exp(-new_segpi));
                set_priors();
                data_llk = llk();
                double tmpdiff = data_llk- target;
                double d2 = tmpdiff*tmpdiff;
                if (d2 > tau) return 1;
                else return 0;
        }
        else if (*nit > 50) return 1;

        segpi = 1.0  / ( 1.0 + exp(-x));
        set_priors();
        double f_x = llk()-target;
        f_x = f_x*f_x;

        segpi =  1.0  / ( 1.0 + exp(-guess));
        set_priors();
        double f_guess = llk()- target;
        f_guess = f_guess*f_guess;

        cout << x << " " <<  guess << " "<< f_x << " "<< f_guess <<  " "<< max << " "<< min << " "<< target << "\n";
        if (f_x < f_guess){
                if ( (max-guess) > (guess-min) )        return golden_section_segpi_ci(guess, x, max, tau, target, nit);
                else return golden_section_segpi_ci(min, x, guess, tau, target, nit);
        }
        else{
                if ( (max - guess) > (guess - min)  ) return golden_section_segpi_ci(min, guess, x, tau, target, nit);
                else return golden_section_segpi_ci(x, guess, max, tau, target, nit);
        }
}


int SNPs::golden_section_l0(double min, double guess, double max, double tau){
        double x;

        if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
        else x = guess - resphi *(guess-min);
        if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
                double new_segpi = (min+max)/2;
                lambdas[0] =  new_segpi;
                data_llk = llk();
                return 0;
        }

        lambdas[0] = x;
        set_priors();
        double f_x = -llk();

        lambdas[0] = guess;
        set_priors();
        double f_guess = -llk();
        cout << x << " " <<  guess << " "<< segpi << " "<< f_x << " "<< f_guess <<  " "<< max << " "<< min << "\n";
        if (f_x < f_guess){
                if ( (max-guess) > (guess-min) )        return golden_section_l0(guess, x, max, tau);
                else return golden_section_l0(min, x, guess, tau);
        }
        else{
                if ( (max - guess) > (guess - min)  ) return golden_section_l0(min, guess, x, tau);
                else return golden_section_l0(x, guess, max, tau);
        }
}

*/
double GSL_llk(const gsl_vector *x, void *params ){
	//first set times
	int na = ((struct GSL_params *) params)->d->alpha.size()-1;

	for (int i = 0; i < na; i++){
		((struct GSL_params *) params)->d->alpha[i+1] = gsl_vector_get(x, i);
	}

	((struct GSL_params *) params)->d->set_priors();
	return -((struct GSL_params *) params)->d->llk();
}

double GSL_llk_fine(const gsl_vector *x, void *params ){
	//first set times
	int na = ((struct GSL_params *) params)->d->alpha.size()-2;

	for (int i = 0; i < na; i++){
		((struct GSL_params *) params)->d->alpha[i+2] = gsl_vector_get(x, i);
	}

	((struct GSL_params *) params)->d->set_priors();
	return -((struct GSL_params *) params)->d->llk();
}

/*
double GSL_llk_xv(const gsl_vector *x, void *params ){
	//first set times
	int na = ((struct GSL_params *) params)->d->nannot;
	int ns = ((struct GSL_params *) params)->d->nsegannot;
	set<int> toskip = ((struct GSL_params *) params)->toskip;
	bool penalize  = ((struct GSL_params *) params)->penalize;
	((struct GSL_params *) params)->d->segpi = 1.0 /  (1.0 + exp (- gsl_vector_get(x, 0)));

	for (int i = 0; i < ns; i++){
		((struct GSL_params *) params)->d->seglambdas[i] = gsl_vector_get(x, i+1);
	}

	for (int i = 0; i < na; i++){
		((struct GSL_params *) params)->d->lambdas[i] = gsl_vector_get(x, i+ns+1);
	}
	((struct GSL_params *) params)->d->set_priors();
	return -((struct GSL_params *) params)->d->llk_xv(toskip, penalize);
}


double GSL_llk_ridge(const gsl_vector *x, void *params ){
	//first set times
	int na = ((struct GSL_params *) params)->d->nannot;
	int ns = ((struct GSL_params *) params)->d->nsegannot;
	((struct GSL_params *) params)->d->segpi = 1.0 /  (1.0 + exp (- gsl_vector_get(x, 0)));

	for (int i = 0; i < ns; i++){
		((struct GSL_params *) params)->d->seglambdas[i] = gsl_vector_get(x, i+1);
	}

	for (int i = 0; i < na; i++){
		((struct GSL_params *) params)->d->lambdas[i] = gsl_vector_get(x, i+ns+1);
	}
	((struct GSL_params *) params)->d->set_priors();
	return -((struct GSL_params *) params)->d->llk_ridge();
}


double GSL_llk_fine(const gsl_vector *x, void *params ){
	//first set times
	int na = ((struct GSL_params *) params)->d->nannot;

	for (int i = 0; i < na; i++){
		((struct GSL_params *) params)->d->lambdas[i] = gsl_vector_get(x, i);
	}
	((struct GSL_params *) params)->d->set_priors();
	return -((struct GSL_params *) params)->d->llk();
}


double GSL_llk_ridge_fine(const gsl_vector *x, void *params ){
	//first set times
	int na = ((struct GSL_params *) params)->d->nannot;

	for (int i = 0; i < na; i++){
		((struct GSL_params *) params)->d->lambdas[i] = gsl_vector_get(x, i);
	}
	((struct GSL_params *) params)->d->set_priors();
	return -((struct GSL_params *) params)->d->llk_ridge();
}


*/
