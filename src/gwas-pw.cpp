/*
 * fgwas.cpp
 *
 *  Created on: Apr 2, 2013
 *      Author: pickrell
 */


#include "SNPs_PW.h"
#include "gwaspw_params.h"
using namespace std;


void printopts(){
        cout << "\ngwas-pw v0.21\n";
        cout << "by Joe Pickrell (jkpickrell@nygenome.org)\n\n";
        cout << "-i [file name] input file w/ Z-scores, variances in beta estimates\n";
        cout << "-phenos [string] [string] names of the phenotypes\n";
        cout << "-o [string] stem for names of output files\n";
        cout << "-bed [file name] read block positions from a .bed file\n";
        cout << "-noprint don't print the Bayes factors\n";
        //cout << "-ld [string] input file with LD matrices\n";
        //cout << "-w [string] which annotation(s) to use. Separate multiple annotations with plus signs\n";
        //cout << "-dists [string:string] the name of the distance annotation(s) and the file(s) containing the distance model(s)\n";
        cout << "-k [integer] block size in number of SNPs (5000)\n";
        cout << "-cor [float] expected correlation in summary statistics under the null [0]\n";
        //cout << "-dists [string:string] the name of the distance annotation(s) and the file(s) containing the distance model(s)\n";
        //cout << "-nburn [integer] iterations of burn-in (5000)\n";
        //cout << "-nsamp [integer] iterations of sampling (50000)\n";
        //cout << "-jumpsd [float] SD of normally distributed MCMC jumps (0.44)\n";
        //cout << "-prior [float] [float] [float] [float] [float] logistic normal prior on fractions (0,0,0,0,0)\n";

        //cout << "-nhap [int] number of haplotypes used in the estimatation of the covariance matrix\n";
        //cout << "-fine do fine-mapping\n";
        //cout << "-mcmc do MCMC\n";
        cout << "\n";
}

bool printcond = false;
int main(int argc, char *argv[]){
	Fgwas_params p;

    CCmdLine cmdline;
    if (cmdline.SplitLine(argc, argv) < 1){
        printopts();
        exit(1);
    }
    //get the input file
    if (cmdline.HasSwitch("-i")) p.infile = cmdline.GetArgument("-i", 0).c_str();
    else{
    	cerr << "ERROR: missing input file (-i)\n";
        printopts();
        exit(1);
    }
    //get the output file
    if (cmdline.HasSwitch("-o")) p.outstem = cmdline.GetArgument("-o", 0);
    if (cmdline.HasSwitch("-pcond")) printcond = true;
    if (cmdline.HasSwitch("-v")) {
     	p.V.clear();
     	vector<string> strs;
     	string s = cmdline.GetArgument("-v", 0);
     	boost::split(strs, s ,boost::is_any_of(","));
     	for (int i  = 0; i < strs.size(); i++) {
     		p.V.push_back( atof(strs[i].c_str()) );
     	}
     }
    //LD file
    if (cmdline.HasSwitch("-ld")) {
    	p.overlap = true;
    	p.ldfile = cmdline.GetArgument("-ld", 0);
    	if (cmdline.HasSwitch("-nhap")) p.Nhap = atoi(cmdline.GetArgument("-nhap", 0).c_str());
    	else{
    		cerr << "ERROR: inputing LD matrix, -nhap flag\n";
    	   	printopts();
    	   	exit(1);
    	}
    }
    else if (cmdline.HasSwitch("-cor")){
    	//cerr << "WARNING: including correlation, did you mean to include an LD file?\n";
    	//printopts();
    	//exit(1);
    }
    if (cmdline.HasSwitch("-rev")) p.rev = true;
    if (cmdline.HasSwitch("-numbered")) p.numberedseg = true;
    //set K
    if (cmdline.HasSwitch("-k")) p.K = atoi(cmdline.GetArgument("-k", 0).c_str());
    if (cmdline.HasSwitch("-bed")) {
    	p.bedseg = true;
    	p.segment_bedfile = cmdline.GetArgument("-bed", 0);
    }

    if (cmdline.HasSwitch("-noprint")) p.print = false;

    //names of the phenotypes, expecting header like NAME1_Z NAME1_V NAME2_Z NAME2_V
    if (cmdline.HasSwitch("-phenos")){
     	p.pairwise = true;
     	p.pheno1 = cmdline.GetArgument("-phenos", 0);
     	p.pheno2 = cmdline.GetArgument("-phenos", 1);
     }
    else{
    	cerr << "ERROR: missing phenotypes (-pheno)\n";
        printopts();
        exit(1);
    }
    /*
    if (cmdline.HasSwitch("-w")){
    	vector<string> strs;
    	string s = cmdline.GetArgument("-w", 0);
    	boost::split(strs, s ,boost::is_any_of("+"));
    	for (int i  = 0; i < strs.size(); i++) {
    		p.wannot.push_back( strs[i] );
    	}
    }
    if (cmdline.HasSwitch("-dists")){
     	vector<string> strs;
     	string s = cmdline.GetArgument("-dists", 0);
     	boost::split(strs, s ,boost::is_any_of("+"));
     	for (int i  = 0; i < strs.size(); i++) {
     		vector<string> strs2;
     		boost::split(strs2, strs[i], boost::is_any_of(":"));
     		p.dannot.push_back( strs2[0] );
     		p.distmodels.push_back(strs2[1]);
     	}
     }
    if (cmdline.HasSwitch("-drop")){
    	p.dropchr = true;
    	string s = cmdline.GetArgument("-drop", 0);
    	p.chrtodrop = s;
    }


    if (cmdline.HasSwitch("-dens")) {
    	p.segannot.push_back(cmdline.GetArgument("-dens", 0));
    	p.loquant = atof(cmdline.GetArgument("-dens", 1).c_str());
    	p.hiquant = atof(cmdline.GetArgument("-dens", 2).c_str());
    }
    */
    if (cmdline.HasSwitch("-fine")) p.finemap = true;
    if (cmdline.HasSwitch("-mcmc")) p.MCMC = true;
    if (cmdline.HasSwitch("-seed")){
    	p.seed = atoi(cmdline.GetArgument("-seed", 0).c_str());
    }
    else p.seed = unsigned( time(NULL));
    if (cmdline.HasSwitch("-nburn")){
     	p.burnin = atoi(cmdline.GetArgument("-nburn", 0).c_str());
    }
    if (cmdline.HasSwitch("-nsamp")){
      	p.nsamp = atoi(cmdline.GetArgument("-nsamp", 0).c_str());
     }
    if (cmdline.HasSwitch("-jumpsd")){
      	p.MCMC_gauss_SD = atof(cmdline.GetArgument("-jumpsd", 0).c_str());
     }
    if (cmdline.HasSwitch("-cor")){
        	p.cor = atof(cmdline.GetArgument("-cor", 0).c_str());
        	//p.overlap = true;
       }
    if (cmdline.HasSwitch("-prior")){
    	if (cmdline.GetArgumentCount("-prior") != 5) {
    		cerr << "ERROR: -prior needs 5 entries, "<< cmdline.GetArgumentCount("-prior") << " given\n";
    		exit(1);
    	}

       	p.alpha_prior[0] = atof(cmdline.GetArgument("-prior", 0).c_str());
       	p.alpha_prior[1] = atof(cmdline.GetArgument("-prior", 1).c_str());
       	p.alpha_prior[2] = atof(cmdline.GetArgument("-prior", 2).c_str());
       	p.alpha_prior[3] = atof(cmdline.GetArgument("-prior", 3).c_str());
       	p.alpha_prior[4] = atof(cmdline.GetArgument("-prior", 4).c_str());
      }


      //random number generator
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_ranlxs2;
    r = gsl_rng_alloc(T);
    int seed = (int) time(0);
    gsl_rng_set(r, p.seed);


    SNPs_PW s(&p);
    if (printcond) s.get_all_condZ();
    s.GSL_optim();
    vector<double> ml;
    for (int i = 0; i < 5; i++)ml.push_back(s.pi[i]);

    vector<pair<pair<int, int>, pair<double, double> > > cis = s.get_cis();
	string outML = p.outstem+".MLE";
	ofstream outr(outML.c_str());
	int sti = 0;
	if (p.finemap) sti = 1;
	for (int i = sti; i < 5; i++){
		outr << "pi_"<< i <<" "<< cis.at(i-sti).second.first << " "<< ml[i]<< " "<< cis.at(i-sti).second.second << "\n";
	}
	outr.close();
    if (p.MCMC) s.MCMC(r);
	if (p.print) s.print(p.outstem+".bfs.gz", p.outstem+".segbfs.gz");
	//if (p.finemap) return 0;


	return 0;
}

