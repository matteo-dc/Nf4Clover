#ifndef OPER_HPP
#define OPER_HPP

#include "aliases.hpp"
#include "global.hpp"
#include "read.hpp"
#include "read_input.hpp"
#include <array>
#include <iostream>
#include <vector>

#ifndef EXTERN_OPER
#define EXTERN_OPER extern
#endif

enum SCHEME_t
{
  RI_MOM,
  SMOM,
  ERR
};

EXTERN_OPER SCHEME_t sch;

struct oper_t
{
  vector<array<int, 1>> linmoms; // list of momenta used for Z, relative to glb list
  vector<array<int, 3>> bilmoms; // lists of momenta used for bilinears {k,i,j}

  void read_mom_list(const string &path);

  // volume
  double V;

  // size of nm and nr
  int _nm;
  int _nr;
  int _nmr;

  // size of linmoms and bilmoms
  int _linmoms;
  int _bilmoms;

  // variables that characterize the data struct
  double _beta;
  string _beta_label;
  int _nm_Sea;
  string _SeaMasses_label;
  double _mu_sea;
  string _volume_label;
  double g2;
  double g2_tilde;
  double _csw;

  // paths
  string ensemble_name;
  string path_to_ens;
  string path_to_beta;
  string path_to_moms;
  string path_print;

  // mom lists
  vector<coords_t> mom_list;
  vector<p_t> p, p_tilde;
  vector<double> p2, p2_tilde; //, p2_tilde_eqmoms;
  vector<double> p4, p4_tilde;
  vector<bool> filt_moms;
  vector<double> Np;

  // number of momenta
  int moms;

  // mPCAC
  void compute_mPCAC(const string &suffix);

  // effective valence mass
  vvvd_t eff_mass;
  vvvd_t eff_mass_corr;
  vvvd_t read_eff_mass(const string name);
  void compute_eff_mass();

  // effective valence mass (time dependent)
  vvvvd_t eff_mass_time;
  vvvvd_t eff_mass_corr_time;
  vvvvd_t read_eff_mass_time(const string name);

  // effective slope
  vd_t effective_slope(vd_t data, vd_t M, int TH);

  // effective sea mass
  vd_t eff_mass_sea;
  vd_t read_eff_mass_sea(const string name);
  void compute_eff_mass_sea();

  // compute the basic RC estimators
  void create_basic(const int b, const int msea);

  void set_ins();

  void set_moms();

  void set_ri_mom_moms();

  void set_smom_moms();

  void clear_all();

  // allocate vectors
  void allocate();
  // allocate valarrays
  void allocate_val();
  // check allocation
  void check_allocation();

  void resize_output(oper_t out);

  vector<string> setup_read_qprop(FILE *input[]);

  void build_prop(const vvvprop_t &prop, vvvprop_t &jprop);

  void ri_mom()
  {
    compute_sigmas();
    compute_bil();
  }

  void smom() { ri_mom(); }

  // create props
  void compute_prop();

  // print averaged quantities
  void print(const string suffix);

  // print extrapolated Z
  void printZ(const string suffix);

  // load averaged quantities and extrapolated Z
  void load(const string suffix);

  // definition of sigmas: trace of propagator
  vector<vvvvd_t> sigma;
  // compute sigmas
  void compute_sigmas();
  vvvd_t compute_sigma(vvvprop_t &jprop_inv, const int iproj, const int imom);

  // definition of projected bils
  vector<jproj_t> jG;
  // compute projected bils
  void compute_bil();

  // definition of jZq
  vector<jZq_t> jZq;
  // compute Zq
  void compute_Zq();

  // definition of Z
  vector<jZbil_t> jZ;
  // definition of ZV/ZA and ZP/ZS
  vector<jZbil_t> jZVoverZA;
  vector<jZbil_t> jZPoverZS;
  vector<jZbil_t> jZAoverZV;

  // compute Zbils
  void compute_Zbil();

  // average r
  oper_t average_r(/*const bool recompute_Zbil=false*/);

  // chiral valence extrapolation
  oper_t chiral_extr();

  // linear coefficients of chiral extrapolations
  vector<vvd_t> bval;
  vector<vvd_t> bsea;

  // step scaling function
  vector<vd_t> stepfunc;

  oper_t subOa2(const int b);

  oper_t subOainf(const int b, const string &suffix);

  oper_t subOainfMartha(const int b, const string &suffix);

  // evolution to 1/a scale
  oper_t evolve(const double ainv, const double p2_evol);
  oper_t evolveToAinv(const double ainv);

  // evolution to 1/a scale of mixed eta
  oper_t evolve_mixed(double ainv);

  // filter momenta
  oper_t filter_moms();
  oper_t filter_momsZ();

  // average of equivalent momenta
  oper_t average_equiv_moms();
  oper_t average_equiv_momsZ();

  // pick spatial-democratic momenta
  oper_t pick_democratic();

  // returns the range in which x is contained
  pair<double, double> get_a2p2tilde_range(const int size, const double a2p2_ref,
                                           const size_t n = 5) const;
  pair<double, double> get_a2p2_range(const int size, const double a2p2_ref,
                                      const size_t n = 5) const;

  // interpolate to p2ref
  oper_t interpolate_to_p2ref_near_neighbour(int b);
  oper_t interpolate_to_p2ref(int b);

  // interpolate to p2=1/a2
  oper_t interpolate_to_ainv_near_neighbour(int b);
  oper_t interpolate_to_ainv(int b);

  // plot Zq and Z
  void plot(const string suffix);
  void plot(const string suffix, int b);
  // plot chiral extrapolations
  void plot_bil_chir_extr(int mom, int ins, int ibil, vd_t x, vd_t y, vd_t dy, vvd_t jpars,
                          string suffix);

  // plot sigmas
  void plot_sigmas();

  // a2p2->0 extrapolation
  oper_t a2p2_extr(int b);
  oper_t a2p2_extr_with_pole(int b);
  oper_t a2p2_extr_with_p4(int b);
  oper_t a2p2_extr_with_pole_and_p4(int b);

  // remove hadronic contamination ~1/p2
  oper_t remove_hadr_cont(double _ainv);

  // improvement of Z w/o hadr. cont.
  oper_t Z_improvement(double _ainv);
};

// valarray of oper_t struct;
using voper_t = valarray<oper_t>;
using vvoper_t = valarray<voper_t>;
using vvvoper_t = valarray<vvoper_t>;
using vvvvoper_t = valarray<vvvoper_t>;
using vvvvvoper_t = valarray<vvvvoper_t>;

// combined a2p2 extrapolation
voper_t combined_M3(voper_t in);
voper_t combined_M4(voper_t in);
voper_t combined_M5(voper_t in);
voper_t combined_M5_quadratic(voper_t in);
voper_t combined_M5_log(voper_t in);

// chiral sea extrapolation
oper_t chiral_sea_extr(voper_t in);

#endif
