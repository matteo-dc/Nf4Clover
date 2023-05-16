#ifndef GLOBAL_HPP
#define GLOBAL_HPP

#include "aliases.hpp"
#include <string>
#include <fstream>
#include <omp.h>

// number of configurations
extern int nconfs;
// number of jackknives
extern int njacks;
// size of the cluster
extern int clust_size;
// sizes
extern coords_t size;
// confs range
extern int conf_init, conf_step;
// number of valence masses
extern int nm;
// number of equivalent bilinear mass combinations
extern int neq;
// number of equivalent Zq mass combinations
extern int neq2;
// number of bilinears
extern int nbil;
// number of Wilson parameters
extern int nr;
// number of mr combinations
extern int nmr;
// number of types
extern int ntypes;
// number of hits
extern int nhits;
// number of combos
extern int combo;
// number of lepton combos
extern int combo_lep;
//// number of momenta
//extern int moms;
// number of equivalent momenta
extern int neqmoms;
// number of colours
extern int Nc;
// number of flavours
extern int Nf;
// string action
extern string action;
// number of betas
extern int nbeta;
// beta
extern vector<double> beta;
// beta label
extern vector<string> beta_label;
// beta label
extern vector<string> volume_label;
// number of theta
extern int ntheta;
// theta label
extern vector<string> theta_label;
// kappa
extern double kappa;
// number of sea masses
extern vector<int> nm_Sea;
// max number of sea masses
extern int nm_Sea_max;
// sea masses label
extern vector<vector<int>> SeaMasses_label;
// sea masses vector
extern vector<vector<double>> SeaMasses;
// sea mass for given ensemble
extern double mu_sea;
// valence masses
extern vector<double> mass_val;
// plaquette
extern double plaquette;
// inverse lattice spacing (Gev)
extern vector<double> ainv;
// clover coefficient
extern vector<double> csw;
// Lambda QCD
extern double LambdaQCD;
// path to upper folder
extern string path_folder;
// path to global ensembles folder
extern string path_ensemble;
// analysis method
extern string analysis;
// path to specific analysis folder
extern vector<string> path_analysis;
// path to the ensemble
extern string path_to_ens;
// add suffix to analysis folder
extern string an_suffix;
// p2min and p2max values for M3 and M4 (Gev^2)
extern vector<double> p2min_M3_M4;
extern vector<double> p2max_M3_M4;
// load label
extern string load_label;

// use Sigma1 parameter
extern int UseSigma1;
// use Effective Mass parameter
extern int UseEffMass;
// N(p) vector
//    extern vector<int> Np;
// path to mom list
extern string mom_path;
// string scheme
extern string scheme;
// range for deltam_cr fit
extern int delta_tmin;
extern int delta_tmax;
// boundary conditions
//extern string BC_str;
extern string BC;
// minimum p2 value of a2p2->0 limit range
extern double p2min;
// maximum p2 value of a2p2->0 limit range
extern double p2max;
// filtered yes/no
extern vector<bool> filt_moms;
// filter threshold
extern double thresh;
// out folder for quarks
extern string out_hadr;
// out folder for leptons
extern string out_lep;
// types of lepton propagators
extern int ntypes_lep;
// create only basic
extern int only_basic;
// compute mPCAC
extern int compute_mpcac;
// bool for inte or free analysis
extern bool free_analysis;
extern bool inte_analysis;
extern bool eta_analysis;
// bool for clover analysis
extern bool clover_analysis;
// bool for QCD factorized on the RIGHT
extern int QCD_on_the_right;
// bool for subtraction with g2 boosted
extern int sub_boosted;
// bool for subtraction with p2tilde
extern int sub_ptilde;
// bool for cut-off subtraction procedure
extern int subtraction;
// bool for recomputing basic
extern bool recompute_basic;
// reference p2
extern double p2ref;
// scales for step scaling function
extern double stepfunc_min;
extern double stepfunc_max;
// load averaged quantities
extern int load_ave;
// load chiral extrapolated quantities
extern int load_chir;
// chiral fit ansatz
extern string chir_ansatz_val;
extern string chir_ansatz_sea;
#endif
