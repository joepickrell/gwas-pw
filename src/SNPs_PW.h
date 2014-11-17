/*
 * SNPs.h
 *
 *  Created on: Jan 14, 2013
 *      Author: pickrell
 */

#ifndef SNPS_PW_H_
#define SNPS_PW_H_

#include "SNP_PW.h"
#include "gwaspw_params.h"
#include "LDmatrix.h"
using namespace std;

class SNPs_PW{
public:
	SNPs_PW();
	SNPs_PW(Fgwas_params *);
	Fgwas_params *params;
	vector<SNP_PW> d;
	void load_snps_pw(string, vector<string>, vector<string>, vector<string>);

	//snp annotations
	vector<vector<double> > snppri; //needs prior influencing trait 1, prior trait 2, prior both
	vector<double> snppost;
	vector<double> lambdas;
	vector<string> annotnames;
	vector<string> chrnames;
	int nannot;
	vector<vector<pair<int, int> > > dmodels; // hold the distance models
	double condlambda; //for conditional analysis

	//for the MCMC, need to set priors
	vector<double> logist_params;
	vector<double> alpha;
	vector<double> pi;
	//segment annotations
	int nsegannot;
	vector<string> segannotnames;
	vector<pair<int, int> > segments;
	vector<pair<int, int> > chrsegments;
	vector<double> seglambdas;
	vector<vector<bool> > segannot;
	vector<vector<double> > segpriors;
	void init_segpriors();
	void set_segpriors();
	map<string, vector<pair< int, int> > > read_bedfile(string);
	void make_segments(string);

	//funtions for working with distance annotations
	void append_dannotnames(string, vector<pair<int, int> >);
	vector<pair<int, int> > read_dmodel(string);

	//pre-computing likelihoods
	vector<vector<double> > seg_toadd; //contains the relative likelihoods for all 4 models
	bool precomputed;

	//MCMC
	void MCMC(gsl_rng *);
	int MCMC_update(gsl_rng *);
	vector<double> propose_alpha(gsl_rng *);
	double dirichlet_lndens(vector<double>, vector<double>);
	double lndgauss(double, double);
	double lnvecdens(vector<double>, vector<double>);
	//10-fold cross-validation
	double cross10(bool);
	vector<set<int> > make_cross10();
	double llk_xv(set<int>, bool);

	double phi, resphi;

	void print();
	void print(string, string);

	void make_segments(int);
	void make_chrsegments();
	void make_segments_finemap();
	void print_segments();
	void print_chrsegments();
	void set_priors(int);
	void set_priors();
	void set_priors_cond(int);
	void set_priors_cond();

	void set_post(int);
	void set_post();
	void GSL_optim();
	void GSL_xv_optim(set<int>, bool);
	void GSL_optim_fine();
	void GSL_optim_ridge();
	void GSL_optim_ridge_fine();
	double llk(int);
	void get_condZ(int, ogzstream &);
	void get_all_condZ();
	double llk();
	double llk_ridge();

	double data_llk;
	double sumlog(double, double);
	void print_segprobs(string);
	void optimize_segpi();
	void optimize_condlambda();
	void optimize_l0();
	vector<pair<pair<int, int>, pair<double, double> > > get_cis();
	pair<pair<int, int>, pair<double, double> > get_cis_segpi();
	pair<pair<int, int>, pair<double, double> > get_cis_condlambda();
	pair<pair<int, int>, pair<double, double> > get_cis_lambda(int);
	pair<pair<int, int>, pair<double, double> > get_cis_alpha(int);
	pair<pair<int, int>, pair<double, double> > get_cis_seglambda(int);

	int golden_section_segpi(double, double, double, double);
	int golden_section_condlambda(double, double, double, double);
	int golden_section_segpi_ci(double, double, double, double, double, int *);
	int golden_section_alpha_ci(double, double, double, double, int, double);
	int golden_section_condlambda_ci(double, double, double, double, double);
	int golden_section_seglambda_ci(double, double, double, double, int, double);
	int golden_section_l0(double, double, double, double);

	void check_input();
	vector<double> get_segbfs(int);

};


struct GSL_params{
        SNPs_PW *d;
};

extern double GSL_llk(const gsl_vector *, void *GSL_params);
extern double GSL_llk_ridge(const gsl_vector *, void *GSL_params);
extern double GSL_llk_fine(const gsl_vector *, void *GSL_params);
extern double GSL_llk_ridge_fine(const gsl_vector *, void *GSL_params);
extern double GSL_llk_xv(const gsl_vector *, void *GSL_params);

#endif /* SNPS_H_ */
