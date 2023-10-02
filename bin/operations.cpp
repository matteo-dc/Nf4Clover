#include "Dirac.hpp"
#include "aliases.hpp"
#include "allocate.hpp"
#include "ave_err.hpp"
#include "evolution.hpp"
#include "fit.hpp"
#include "global.hpp"
#include "print.hpp"
#include "prop.hpp"
#include "read.hpp"
#include "sigmas.hpp"
#include "subtraction.hpp"
#include "vertices.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

// #define EXTERN_OPER
#include "operations.hpp"

using namespace std::chrono;

SCHEME_t get_scheme()
{
  if (scheme == "RI_MOM")
    return RI_MOM;
  if (scheme == "SMOM")
    return SMOM;

  return ERR;
}

void oper_t::set_ins()
{
  sigma::set_ins();
  jprop::set_ins();
  qprop::set_ins();
  gbil::set_ins();
}

void oper_t::set_moms()
{
  // read mom list
  read_mom_list(path_to_moms);
  moms = mom_list.size();
  cout << "Read: " << moms << " momenta from \"" << mom_path << "\" (BC: " << BC << ")." << endl;

  switch (get_scheme())
  {
  case RI_MOM:
    set_ri_mom_moms();
    break;
  case SMOM:
    cout << "SMOM!" << endl;
    set_smom_moms();
    break;
  case ERR:
    cout << "Invalid scheme." << endl;
    exit(0);
    break;
  }
  _linmoms = linmoms.size();
  _bilmoms = bilmoms.size();

  cout << _linmoms << " linmoms, " << _bilmoms << " bilmoms." << endl;

  print_vec(p2_tilde, path_print + "p2_tilde.txt");
  print_vec(p2, path_print + "p2.txt");
}

void oper_t::set_ri_mom_moms()
{
  linmoms.resize(moms);
  bilmoms.resize(moms);

  for (int imom = 0; imom < moms; imom++)
  {
    linmoms[imom] = {imom};
    bilmoms[imom] = {imom, imom, imom};
  }
}

void oper_t::set_smom_moms()
{
  cout << " moms not initialized for SMOM." << endl;
  exit(0);
}

void oper_t::create_basic(const int b, const int msea)
{
  //    step = "basic";
  cout << endl;
  cout << "----- basic -----" << endl << endl;

  _beta = beta[b];
  _beta_label = beta_label[b];
  _volume_label = volume_label[b];
  _nm_Sea = nm_Sea[b];
  _SeaMasses_label = to_string(SeaMasses_label[b][msea]);
  _mu_sea = SeaMasses[b][msea];
  _csw = csw[b];

  if (inte_analysis)
  {
    // e.g. /.../matteo/Nf4/
    path_ensemble = path_folder + path_analysis[0] + "/";

    // e.g. /.../matteo/Nf4_Clover/C.d.50.32/
    // if (!(strcmp(an_suffix.c_str(), "_Clover") == 0 or
    //       strcmp(an_suffix.c_str(), "_Clover2023") == 0))
    // {
    //   cout << "Suffix '_Clover' needed!" << endl;
    //   exit(0);
    // }

    path_to_beta = path_ensemble;

    // e.g. B1m
    ensemble_name = _beta_label + ".d." + _SeaMasses_label + "." + _volume_label;

    // e.g. /.../matteo/Nf4/B_b1.95/B1m/
    path_to_ens = path_to_beta + ensemble_name + "/";
  }
  else
  {
    cout << "Only interacting analysis implemented." << endl;
    exit(0);
  }

  cout << path_to_ens << endl;
  read_input(path_to_ens, ensemble_name);
  path_to_moms = path_to_ens + mom_path;

  path_print = path_to_ens + "print/";

  V = sizeV[0] * sizeV[1] * sizeV[2] * sizeV[3];

  _nm = nm;
  _nr = nr;
  _nmr = _nm * _nr;

  g2 = 6.0 / _beta;
  g2_tilde = g2 / plaquette;

  set_moms();

  set_ins();

  allocate_val();
  if (!load_ave)
    allocate();

  if (compute_mpcac)
  {
    compute_mPCAC("");
    compute_mPCAC("sea");
  }

  if (UseEffMass)
  {
    eff_mass = read_eff_mass(path_to_ens + "eff_mass_array");
    eff_mass_time = read_eff_mass_time(path_to_ens + "eff_mass_array_time");

    if (_nm_Sea > 1)
      eff_mass_sea = read_eff_mass_sea(path_to_ens + "eff_mass_sea_array");
  }

  if (!load_ave)
  {
    switch (get_scheme())
    {
    case RI_MOM:
      ri_mom();
      break;
    case SMOM:
      cout << "SMOM!" << endl;
      smom();
      break;
    case ERR:
      cout << "Invalid scheme." << endl;
      exit(0);
      break;
    }

    compute_Zq();
    compute_Zbil();
  }
}

oper_t oper_t::average_r()
{
  cout << endl;
  cout << "----- r average -----" << endl << endl;

  oper_t out = (*this);

  out._nr = 1;
  out._nm = _nm;
  out._nmr = (out._nm) * (out._nr);

  out.allocate_val();

  out.eff_mass = eff_mass;
  out.eff_mass_sea = eff_mass_sea;

  if (!load_ave)
  {
    out.allocate();

#pragma omp parallel for collapse(5)
    for (int ilinmom = 0; ilinmom < _linmoms; ilinmom++)
      for (int iproj = 0; iproj < sigma::nproj; iproj++)
        for (int ins = 0; ins < sigma::nins; ins++)
          for (int m = 0; m < _nm; m++)
            for (int ijack = 0; ijack < njacks; ijack++)
              for (int r = 0; r < _nr; r++)
              {
                const int coeff = (iproj == sigma::SIGMA3 and r == 1) ? -1 : +1;

                (out.sigma)[ilinmom][iproj][ins][ijack][m] +=
                    coeff * sigma[ilinmom][iproj][ins][ijack][r + _nr * m] / _nr;
              }

    out.compute_Zq();

#pragma omp parallel for collapse(6)
    for (int ibilmom = 0; ibilmom < _bilmoms; ibilmom++)
      for (int mA = 0; mA < _nm; mA++)
        for (int mB = 0; mB < _nm; mB++)
          for (int ijack = 0; ijack < njacks; ijack++)
            for (int ibil = 0; ibil < 5; ibil++)
              for (int ins = 0; ins < gbil::nins; ins++)
                for (int r = 0; r < _nr; r++)
                  (out.jG)[ibilmom][ins][ibil][ijack][mA][mB] +=
                      jG[ibilmom][ins][ibil][ijack][r + _nr * mA][r + _nr * mB] / _nr;

    out.compute_Zbil();
  }

  return out;
}

oper_t oper_t::chiral_extr()
{
  cout << endl;
  cout << "----- valence chiral extrapolation -----" << endl << endl;

  oper_t out = (*this);

  out._nr = _nr;
  out._nm = 1;
  out._nmr = (out._nm) * (out._nr);

  out.allocate_val();
  out.allocate();

  out.eff_mass = eff_mass;
  out.eff_mass_sea = eff_mass_sea;

  // average of eff_mass
#warning Da cancellare quando si useranno davvero le eff_mass
  vvd_t M_eff;
  vvd_t dM_eff;
  vvd_t M2_eff;
  vvd_t M2_eff_err;

  if (UseEffMass)
  {
    vvd_t M_eff = get<0>(ave_err(eff_mass));
    vvd_t dM_eff = get<0>(ave_err(eff_mass_corr));
    for (int m = 0; m < _nm; m++)
      printf("m: %d \t M_eff: %lg \t dM_eff: %lg\n", m, M_eff[m][m], dM_eff[m][m]);

    // printing the squared mass
    vvvd_t eff_mass_sqr(vvd_t(vd_t(0.0, _nm), _nm), njacks);

    for (int ijack = 0; ijack < njacks; ijack++)
      for (int m1 = 0; m1 < _nm; m1++)
        for (int m2 = 0; m2 < _nm; m2++)
          eff_mass_sqr[ijack][m1][m2] = eff_mass[ijack][m1][m2] * eff_mass[ijack][m1][m2];

    vvd_t M2_eff = get<0>(ave_err(eff_mass_sqr));
    vvd_t M2_eff_err = get<1>(ave_err(eff_mass_sqr));

    for (int m = 0; m < _nm; m++)
      printf("m: %d \t M2_eff: %lg +- %lg\n", m, M2_eff[m][m], M2_eff_err[m][m]);
  }

  // range for fit Zq
  int x_min_q = 0;
  int x_max_q = _nm - 1;

  // range for fit bilinears
  int x_min = 0;
  int x_max = _nm * (_nm + 1) / 2 - 1;
  int x_min_P = 0;
  int x_max_P = _nm - 1;

  // number of fit parameters for sigma
  int npar_sigma = 2;
  // number of fit parameters for bilinears
  // int npar_bil[5]={3,2,3,2,2};
  int npar_bil[5] = {2, 2, 3, 2, 2};
  int npar_bil_max = *max_element(npar_bil, npar_bil + sizeof(npar_bil) / sizeof(npar_bil[0]));

  bool linear = true, constant = false, quadratic = false;

  if (strcmp(chir_ansatz_val.c_str(), "constant") == 0)
  {
    linear = false;
    constant = true;

    npar_sigma--;

    for (int ibil = 0; ibil < nbil; ibil++)
      npar_bil[ibil]--;
    npar_bil_max--;
  }
  else if (strcmp(chir_ansatz_val.c_str(), "quadratic") == 0)
  {
    linear = false;
    quadratic = true;

    npar_sigma++;

    for (int ibil = 0; ibil < nbil; ibil++)
      npar_bil[ibil]++;
    npar_bil_max++;
  }

  // extrapolate sigma

  vvd_t sigma_pars_QCD(vd_t(0.0, npar_sigma), njacks);

  // #pragma omp parallel for
  for (int ilinmom = 0; ilinmom < _linmoms; ilinmom++)
    for (int iproj = 0; iproj < sigma::nproj; iproj++)
      for (int r = 0; r < _nr; r++)
        for (int ins = 0; ins < sigma::nins; ins++) // not collapsed
        {
          vvd_t coord_sigma(vd_t(0.0, _nm), npar_sigma);

          vvd_t sigma_r(vd_t(0.0, _nm), njacks);
          vd_t sig_ave_r(0.0, _nmr), sqr_sig_ave_r(0.0, _nmr), sig_err_r(0.0, _nmr);

          for (int m = 0; m < _nm; m++)
          {
            int mr = r + _nr * m;

            if (!UseEffMass)
            {
              if (constant)
              {
                coord_sigma[0][m] = 1.0;
              }
              else if (linear)
              {
                coord_sigma[0][m] = 1.0;
                coord_sigma[1][m] = mass_val[m];
              }
              else if (quadratic)
              {
                coord_sigma[0][m] = 1.0;
                coord_sigma[1][m] = mass_val[m];
                coord_sigma[2][m] = mass_val[m] * mass_val[m];
              }
            }
            else if (UseEffMass)
            {
              if (!linear)
              {
                cout << "Only linear fit implemented when using EffMass!" << endl;
                exit(0);
              }

              coord_sigma[0][m] = 1.0;
              coord_sigma[1][m] = pow(M_eff[m][m], 2.0);
            }

            for (int ijack = 0; ijack < njacks; ijack++)
            {
              sigma_r[ijack][m] = sigma[ilinmom][iproj][ins][ijack][mr];

              sig_ave_r[m] += sigma[ilinmom][iproj][ins][ijack][mr] / njacks;
              sqr_sig_ave_r[m] += sigma[ilinmom][iproj][ins][ijack][mr] *
                                  sigma[ilinmom][iproj][ins][ijack][mr] / njacks;
            }
            sig_err_r[m] = sqrt((double)(njacks - 1)) *
                           sqrt(fabs(sqr_sig_ave_r[mr] - sig_ave_r[mr] * sig_ave_r[mr]));
          }

          vvd_t sigma_pars = polyfit(coord_sigma, npar_sigma, sig_err_r, sigma_r, x_min_q, x_max_q);

          // save fit parameters to be used to subtract dM
          if (ins == sigma::LO)
            for (int ijack = 0; ijack < njacks; ijack++)
              sigma_pars_QCD[ijack] = sigma_pars[ijack];

          for (int ijack = 0; ijack < njacks; ijack++)
          {
            (out.sigma)[ilinmom][iproj][ins][ijack][r] = sigma_pars[ijack][0];

            if (iproj == 0 and ins == 0 and r == 0 and linear)
              (out.bval)[ilinmom][ijack][0] = sigma_pars[ijack][1];
          }
        }

  out.compute_Zq();

  // extrapolate bilinears

  vvd_t gbil_pars_QCD(vd_t(0.0, npar_bil_max), njacks);

  // #pragma omp parallel for //collapse(4)
  for (int ibilmom = 0; ibilmom < _bilmoms; ibilmom++)
    for (int r1 = 0; r1 < _nr; r1++)
      for (int r2 = 0; r2 < _nr; r2++)
        for (int ibil = 0; ibil < nbil; ibil++)
          for (int ins = 0; ins < gbil::nins; ins++) // not collapsed
          {
            if (ibil != 2 and !quadratic) /* not for P and only if not quadratic */
            {
              vd_t x_bil(0.0, _nm * (_nm + 1) / 2);

              vvd_t coord_bil(vd_t(0.0, _nm * (_nm + 1) / 2), npar_bil_max);

              vvd_t jG_r1_r2(vd_t(0.0, _nm * (_nm + 1) / 2), njacks);

              vd_t G_ave_r1_r2(0.0, _nm * (_nm + 1) / 2);
              vd_t sqr_G_ave_r1_r2(0.0, _nm * (_nm + 1) / 2);
              vd_t G_err_r1_r2(0.0, _nm * (_nm + 1) / 2);

              int ieq = 0;
              for (int m1 = 0; m1 < _nm; m1++)
                for (int m2 = m1; m2 < _nm; m2++)
                {

                  int mr1 = r1 + _nr * m1;
                  int mr2 = r2 + _nr * m2;

                  if (!UseEffMass)
                  {
                    x_bil[ieq] = mass_val[m1] + mass_val[m2];

                    if (constant)
                    {
                      coord_bil[0][ieq] = 1.0;                                 // 1
                      coord_bil[1][ieq] = 1.0 / (mass_val[m1] + mass_val[m2]); // 1/(am1+am2)
                    }
                    else if (linear)
                    {
                      coord_bil[0][ieq] = 1.0;                         // 1
                      coord_bil[1][ieq] = mass_val[m1] + mass_val[m2]; // (am1+am2)
                      coord_bil[2][ieq] = 1.0 / coord_bil[1][ieq];     // 1/(am1+am2)
                    }
                    else if (quadratic)
                    {
                      coord_bil[0][ieq] = 1.0;                                   // 1
                      coord_bil[1][ieq] = mass_val[m1] + mass_val[m2];           // (am1+am2)
                      coord_bil[2][ieq] = coord_bil[1][ieq] * coord_bil[1][ieq]; // (am1+am2)^2
                      coord_bil[3][ieq] = 1.0 / coord_bil[1][ieq];               // 1/(am1+am2)
                    }
                  }
                  else if (UseEffMass)
                  {
                    x_bil[ieq] = pow((M_eff[m1][m2] + M_eff[m2][m1]) / 2.0, 2.0);

                    if (!linear)
                    {
                      cout << "Only linear fit implemented when using EffMass!" << endl;
                      exit(0);
                    }

                    coord_bil[0][ieq] = 1.0;
                    // M^2 (averaged over equivalent combinations)
                    coord_bil[1][ieq] = pow((M_eff[m1][m2] + M_eff[m2][m1]) / 2.0, 2.0);
                    // 1/M^2
                    coord_bil[2][ieq] = 1.0 / coord_bil[1][ieq];
                  }

                  for (int ijack = 0; ijack < njacks; ijack++)
                  {
                    jG_r1_r2[ijack][ieq] = jG[ibilmom][ins][ibil][ijack][mr1][mr2];

                    G_ave_r1_r2[ieq] += jG_r1_r2[ijack][ieq] / njacks;
                    sqr_G_ave_r1_r2[ieq] += jG_r1_r2[ijack][ieq] * jG_r1_r2[ijack][ieq] / njacks;
                  }
                  G_err_r1_r2[ieq] =
                      sqrt((double)(njacks - 1)) *
                      sqrt(fabs(sqr_G_ave_r1_r2[ieq] - G_ave_r1_r2[ieq] * G_ave_r1_r2[ieq]));

                  ieq++;
                }

              vvd_t jG_pars =
                  polyfit(coord_bil, npar_bil[ibil], G_err_r1_r2, jG_r1_r2, x_min, x_max);

              // save fit parameters to be used to subtract dM
              if (ins == gbil::LO)
                for (int ijack = 0; ijack < njacks; ijack++)
                {
                  gbil_pars_QCD[ijack][0] = jG_pars[ijack][0];
                  gbil_pars_QCD[ijack][1] = jG_pars[ijack][1];

                  if (npar_bil[ibil] < npar_bil_max)
                    gbil_pars_QCD[ijack][2] = 0.0;
                  else
                    gbil_pars_QCD[ijack][2] = jG_pars[ijack][2];
                }

              for (int ijack = 0; ijack < njacks; ijack++)
              {
                (out.jG)[ibilmom][ins][ibil][ijack][r1][r2] = jG_pars[ijack][0];

                if (ins == 0 and r1 == 0 and r2 == 0 and linear)
                  (out.bval)[ibilmom][ijack][ibil + 1] = jG_pars[ijack][1];
              }

              if (ibilmom % 20 == 0 and r1 == 0 and r2 == 0)
              // if(ibilmom<20 and r1==0 and r2==0)
              {
                plot_bil_chir_extr(ibilmom, ins, ibil, x_bil, G_ave_r1_r2, G_err_r1_r2, jG_pars,
                                   ""); /* (mom,ins,bil,x,y,dy,jpars) */
              }
            }
            else /* for P of if quadratic */
            {
              vd_t x_bil(0.0, _nm);

              vvd_t coord_bil(vd_t(0.0, _nm), npar_bil_max);

              vvd_t jG_r1_r2(vd_t(0.0, _nm), njacks);

              vd_t G_ave_r1_r2(0.0, _nm);
              vd_t sqr_G_ave_r1_r2(0.0, _nm);
              vd_t G_err_r1_r2(0.0, _nm);

              for (int m1 = 0; m1 < _nm; m1++)
              // for(int m1=0; m1<4; m1++)
              {
                int m2 = m1;

                int mr1 = r1 + _nr * (m1);
                int mr2 = r2 + _nr * (m2);

                if (!UseEffMass)
                {
                  x_bil[m1] = mass_val[m1] + mass_val[m2];

                  if (constant)
                  {
                    coord_bil[0][m1] = 1.0;                                 // 1
                    coord_bil[1][m1] = 1.0 / (mass_val[m1] + mass_val[m2]); // 1/(am1+am2)
                  }
                  else if (linear)
                  {
                    coord_bil[0][m1] = 1.0;                         // 1
                    coord_bil[1][m1] = mass_val[m1] + mass_val[m2]; // (am1+am2)
                    coord_bil[2][m1] = 1.0 / coord_bil[1][m1];      // 1/(am1+am2)
                  }
                  else if (quadratic)
                  {
                    coord_bil[0][m1] = 1.0;                                 // 1
                    coord_bil[1][m1] = mass_val[m1] + mass_val[m2];         // (am1+am2)
                    coord_bil[2][m1] = coord_bil[1][m1] * coord_bil[1][m1]; // (am1+am2)^2
                    coord_bil[3][m1] = 1.0 / coord_bil[1][m1];              // 1/(am1+am2)
                  }
                }
                // else if(UseEffMass)
                // {
                //     x_bil[m1] = pow((M_eff[m1][m2]+M_eff[m2][m1])/2.0,2.0);
                //
                //     if(!linear)
                //     {
                //         cout<<"Only linear fit implemented when using EffMass!"<<endl;
                //         exit(0);
                //     }
                //
                //     coord_bil[0][m1] = 1.0;
                //     // M^2 (averaged over equivalent combinations)
                //     coord_bil[1][m1] = pow((M_eff[m1][m2]+M_eff[m2][m1])/2.0,2.0);
                //     // 1/M^2
                //     coord_bil[2][m1] = 1.0/coord_bil[1][m1];
                // }

                for (int ijack = 0; ijack < njacks; ijack++)
                {
                  jG_r1_r2[ijack][m1] = jG[ibilmom][ins][ibil][ijack][mr1][mr2];

                  G_ave_r1_r2[m1] += jG_r1_r2[ijack][m1] / njacks;
                  sqr_G_ave_r1_r2[m1] += jG_r1_r2[ijack][m1] * jG_r1_r2[ijack][m1] / njacks;
                }
                G_err_r1_r2[m1] =
                    sqrt((double)(njacks - 1)) *
                    sqrt(fabs(sqr_G_ave_r1_r2[m1] - G_ave_r1_r2[m1] * G_ave_r1_r2[m1]));
              } // masses loop

              vvd_t jG_pars =
                  polyfit(coord_bil, npar_bil[ibil], G_err_r1_r2, jG_r1_r2, x_min_P, x_max_P);

              // save fit parameters to be used to subtract dM
              if (ins == gbil::LO)
                for (int ijack = 0; ijack < njacks; ijack++)
                {
                  gbil_pars_QCD[ijack][0] = jG_pars[ijack][0];
                  gbil_pars_QCD[ijack][1] = jG_pars[ijack][1];

                  if (npar_bil[ibil] < npar_bil_max)
                    gbil_pars_QCD[ijack][2] = 0.0;
                  else
                    gbil_pars_QCD[ijack][2] = jG_pars[ijack][2];
                }

              for (int ijack = 0; ijack < njacks; ijack++)
              {
                (out.jG)[ibilmom][ins][ibil][ijack][r1][r2] = jG_pars[ijack][0];

                if (ins == 0 and r1 == 0 and r2 == 0 and linear)
                  (out.bval)[ibilmom][ijack][ibil + 1] = jG_pars[ijack][1];
              }

              // if(ibilmom%20==0 and r1==0 and r2==0)
              if (ibilmom < 20 and r1 == 0 and r2 == 0)
              {
                plot_bil_chir_extr(ibilmom, ins, ibil, x_bil, G_ave_r1_r2, G_err_r1_r2, jG_pars,
                                   ""); /* (mom,ins,bil,x,y,dy,jpars) */
              }
            }
          }

  out.compute_Zbil();

  return out;
}

oper_t chiral_sea_extr(voper_t in)
{
  cout << endl;
  cout << "----- chiral sea extrapolation -----" << endl << endl;

  oper_t out = in[0]; //?

  int nmSea = in[0]._nm_Sea;
  int _linmoms = in[0]._linmoms;
  int _bilmoms = in[0]._bilmoms;

  out.allocate_val();
  out.allocate();

  out.path_to_ens = in[0].path_to_beta + in[0]._beta_label + "/";

  vd_t x(0.0, nmSea);

  // range for fit
  int x_min = 0;
  int x_max = nmSea - 1;

  for (int msea = 0; msea < nmSea; msea++)
  {
    if (UseEffMass)
      x[msea] = get<0>(ave_err(in[msea].eff_mass_sea));
    else
      x[msea] = in[msea]._mu_sea;
  }

  // number of fit parameters for sigma
  int npar_sigma = 2;
  // number of fit parameters for bilinears
  int npar_bil[5] = {2, 2, 2, 2, 2};
  int npar_bil_max = *max_element(npar_bil, npar_bil + sizeof(npar_bil) / sizeof(npar_bil[0]));

  bool linear = true, constant = false, quadratic = false;

  if (strcmp(chir_ansatz_sea.c_str(), "constant") == 0)
  {
    linear = false;
    constant = true;

    npar_sigma--;

    for (int ibil = 0; ibil < nbil; ibil++)
      npar_bil[ibil]--;
    npar_bil_max--;
  }
  else if (strcmp(chir_ansatz_sea.c_str(), "quadratic") == 0)
  {
    linear = false;
    quadratic = true;

    npar_sigma++;

    for (int ibil = 0; ibil < nbil; ibil++)
      npar_bil[ibil]++;
    npar_bil_max++;
  }

  // extrapolate Zq
#pragma omp parallel for
  for (int ilinmom = 0; ilinmom < _linmoms; ilinmom++)
  {
    vvd_t coord_Zq(vd_t(0.0, nmSea), npar_sigma); // coords at fixed r

    vvd_t y_Zq(vd_t(0.0, nmSea), njacks);
    vd_t dy_Zq(0.0, nmSea);

    for (int msea = 0; msea < nmSea; msea++)
    {
      coord_Zq[0][msea] = 1.0;

      if (!UseEffMass)
      {
        if (constant)
        {
          coord_Zq[0][msea] = 1.0;
        }
        else if (linear)
        {
          coord_Zq[0][msea] = 1.0;
          coord_Zq[1][msea] = x[msea];
        }
        else if (quadratic)
        {
          coord_Zq[0][msea] = 1.0;
          coord_Zq[1][msea] = x[msea];
          coord_Zq[2][msea] = x[msea] * x[msea];
        }
      }
      else if (UseEffMass)
      {
        if (!linear)
        {
          cout << "Only linear fit implemented when using EffMass!" << endl;
          exit(0);
        }

        coord_Zq[1][msea] = pow(x[msea], 2.0);
      }

      for (int ijack = 0; ijack < njacks; ijack++)
        y_Zq[ijack][msea] = in[msea].jZq[ilinmom][ijack][0];

      dy_Zq[msea] = (get<1>(ave_err_Zq(in[msea].jZq)))[ilinmom][0];
    }

    vvd_t Zq_pars = polyfit(coord_Zq, npar_sigma, dy_Zq, y_Zq, x_min, x_max);

    for (int ijack = 0; ijack < njacks; ijack++)
    {
      (out.jZq)[ilinmom][ijack][0] = Zq_pars[ijack][0];

      if (linear)
        (out.bsea)[ilinmom][ijack][0] = Zq_pars[ijack][1];
    }
  }

  // extrapolate bilinears
#pragma omp parallel for collapse(2)
  for (int ibilmom = 0; ibilmom < _bilmoms; ibilmom++)
    for (int ibil = 0; ibil < nbil; ibil++)
    {
      vvd_t coord_bil(vd_t(0.0, nmSea), npar_bil_max);

      vvd_t y_Z(vd_t(0.0, nmSea), njacks);
      vd_t dy_Z(0.0, nmSea);
      vd_t y_Z_ave(0.0, nmSea);

      for (int msea = 0; msea < nmSea; msea++)
      {
        coord_bil[0][msea] = 1.0;
        if (!UseEffMass)
        {
          if (constant)
          {
            coord_bil[0][msea] = 1.0; // 1
          }
          else if (linear)
          {
            coord_bil[0][msea] = 1.0;     // 1
            coord_bil[1][msea] = x[msea]; // (amsea)
          }
          else if (quadratic)
          {
            coord_bil[0][msea] = 1.0;               // 1
            coord_bil[1][msea] = x[msea];           // (amsea)
            coord_bil[2][msea] = x[msea] * x[msea]; // (amsea)^2
          }
        }
        else if (UseEffMass)
        {
          if (!linear)
          {
            cout << "Only linear fit implemented when using EffMass!" << endl;
            exit(0);
          }

          coord_bil[1][msea] = pow(x[msea], 2.0);
        }

        for (int ijack = 0; ijack < njacks; ijack++)
          y_Z[ijack][msea] = in[msea].jZ[ibilmom][ibil][ijack][0][0];

        y_Z_ave[msea] = (get<0>(ave_err_Z(in[msea].jZ)))[ibilmom][ibil][0][0];
        dy_Z[msea] = (get<1>(ave_err_Z(in[msea].jZ)))[ibilmom][ibil][0][0];
      }

      vvd_t jZ_pars = polyfit(coord_bil, npar_bil_max, dy_Z, y_Z, x_min, x_max);

      for (int ijack = 0; ijack < njacks; ijack++)
      {
        (out.jZ)[ibilmom][ibil][ijack][0][0] = jZ_pars[ijack][0];

        if (linear)
          (out.bsea)[ibilmom][ijack][ibil + 1] = jZ_pars[ijack][1];
      }

      if (ibilmom % 1 == 0)
      {
        out.plot_bil_chir_extr(ibilmom, 0, ibil, x, y_Z_ave, dy_Z, jZ_pars,
                               "sea"); /* (mom,ins,bil,x,y,dy,jpars) */
      }
    }

    // extrapolate ZV/ZA and ZP/ZS and ZA/ZV
#pragma omp parallel for collapse(2)
  for (int ibilmom = 0; ibilmom < _bilmoms; ibilmom++)
    for (int iloop = 0; iloop < 3; iloop++)
    {
      vvd_t coord_bil(vd_t(0.0, nmSea), npar_bil_max);

      vvd_t y_Z(vd_t(0.0, nmSea), njacks);
      vd_t dy_Z(0.0, nmSea);
      vd_t y_Z_ave(0.0, nmSea);

      vector<jZbil_t> jzz;

      for (int msea = 0; msea < nmSea; msea++)
      {
        coord_bil[0][msea] = 1.0;
        if (!UseEffMass)
        {
          if (constant)
          {
            coord_bil[0][msea] = 1.0; // 1
          }
          else if (linear)
          {
            coord_bil[0][msea] = 1.0;     // 1
            coord_bil[1][msea] = x[msea]; // (amsea)
          }
          else if (quadratic)
          {
            coord_bil[0][msea] = 1.0;               // 1
            coord_bil[1][msea] = x[msea];           // (amsea)
            coord_bil[2][msea] = x[msea] * x[msea]; // (amsea)^2
          }
        }
        else if (UseEffMass)
        {
          if (!linear)
          {
            cout << "Only linear fit implemented when using EffMass!" << endl;
            exit(0);
          }

          coord_bil[1][msea] = pow(x[msea], 2.0);
        }

        if (iloop == 0)
          jzz = in[msea].jZVoverZA;
        else if (iloop == 1)
          jzz = in[msea].jZPoverZS;
        else
          jzz = in[msea].jZAoverZV;

        for (int ijack = 0; ijack < njacks; ijack++)
          y_Z[ijack][msea] = jzz[ibilmom][0][ijack][0][0];

        y_Z_ave[msea] = (get<0>(ave_err_Z(jzz)))[ibilmom][0][0][0];
        dy_Z[msea] = (get<1>(ave_err_Z(jzz)))[ibilmom][0][0][0];
      }

      vvd_t jZ_pars = polyfit(coord_bil, npar_bil_max, dy_Z, y_Z, x_min, x_max);

      for (int ijack = 0; ijack < njacks; ijack++)
      {
        if (iloop == 0)
          (out.jZVoverZA)[ibilmom][0][ijack][0][0] = jZ_pars[ijack][0];
        else if (iloop == 1)
          (out.jZPoverZS)[ibilmom][0][ijack][0][0] = jZ_pars[ijack][0];
        else
          (out.jZAoverZV)[ibilmom][0][ijack][0][0] = jZ_pars[ijack][0];
      }

      //            if(ibilmom%20==0)
      //            {
      //                out.plot_bil_chir_extr(ibilmom,0,ibil,x,y_Z_ave,dy_Z,jZ_pars,"sea");   /*
      //                (mom,ins,bil,x,y,dy,jpars) */
      //            }
    }

  return out;
}

int mom_list_xyz(vector<coords_t> &mom_list, const size_t imom)
{
  return abs(mom_list[imom][1]) * abs(mom_list[imom][2]) * abs(mom_list[imom][3]);
}

/////////
oper_t oper_t::filter_moms()
{
  cout << endl;
  cout << "Filtering democratic momenta -- ";

  oper_t out = (*this);

  int count_filtered = 0;

  // number of filtered moms
  for (int imom = 0; imom < _linmoms; imom++)
    if (filt_moms[imom])
      count_filtered++;
  cout << "found: " << count_filtered << " filtered linmoms." << endl;

  (out.mom_list).resize(count_filtered);
  (out.p).resize(count_filtered);
  (out.p_tilde).resize(count_filtered);
  (out.Np).resize(count_filtered);

  out._linmoms = count_filtered;
  (out.linmoms).resize(out._linmoms);
  out._bilmoms = count_filtered;
  (out.bilmoms).resize(out._bilmoms);

  (out.p2).resize(out._linmoms);
  (out.p2_tilde).resize(out._linmoms);
  (out.p4).resize(out._linmoms);
  (out.p4_tilde).resize(out._linmoms);

  // fill the new linmoms, p2(tilde), p4(tilde)
  int ifilt = 0;
  for (int imom = 0; imom < _linmoms; imom++)
    if (filt_moms[imom])
    {
      (out.mom_list)[ifilt] = mom_list[imom];
      (out.p)[ifilt] = p[imom];
      (out.p_tilde)[ifilt] = p_tilde[imom];
      (out.Np)[ifilt] = Np[imom];

      (out.linmoms)[ifilt] = {ifilt};
      (out.bilmoms)[ifilt] = {ifilt, ifilt, ifilt};
      if (scheme == "SMOM")
      {
        cout << "Filter not implemented in SMOM." << endl;
        exit(0);
      }

      (out.p2)[ifilt] = p2[imom];
      (out.p2_tilde)[ifilt] = p2_tilde[imom];
      (out.p4)[ifilt] = p4[imom];
      (out.p4_tilde)[ifilt] = p4_tilde[imom];

      ifilt++;
    }

  print_vec(out.p2, path_print + "p2_filtmoms.txt");
  print_vec(out.p2_tilde, path_print + "p2_tilde_filtmoms.txt");

  out.allocate_val();

  out.eff_mass = eff_mass;
  out.eff_mass_sea = eff_mass_sea;

  if (!load_ave)
  {
    out.allocate();

    ifilt = 0;
    for (int imom = 0; imom < _linmoms; imom++)
      if (filt_moms[imom])
      {
        (out.sigma)[ifilt] = sigma[imom];
        ifilt++;
      }

    out.compute_Zq();

    ifilt = 0;
    for (int imom = 0; imom < _linmoms; imom++)
      if (filt_moms[imom])
      {
        (out.jG)[ifilt] = jG[imom];
        ifilt++;
      }

    out.compute_Zbil();
  }
  return out;
}

oper_t oper_t::filter_momsZ()
{
  cout << endl;
  cout << "Filtering democratic momenta on Z's -- ";

  oper_t out = (*this);

  int count_filtered = 0;

  // number of filtered moms
  for (int imom = 0; imom < _linmoms; imom++)
    if (filt_moms[imom])
      count_filtered++;
  cout << "found: " << count_filtered << " filtered linmoms." << endl;

  (out.mom_list).resize(count_filtered);
  (out.p).resize(count_filtered);
  (out.p_tilde).resize(count_filtered);
  (out.Np).resize(count_filtered);

  out._linmoms = count_filtered;
  (out.linmoms).resize(out._linmoms);
  out._bilmoms = count_filtered;
  (out.bilmoms).resize(out._bilmoms);

  (out.p2).resize(out._linmoms);
  (out.p2_tilde).resize(out._linmoms);
  (out.p4).resize(out._linmoms);
  (out.p4_tilde).resize(out._linmoms);

  // fill the new linmoms, p2(tilde), p4(tilde)
  int ifilt = 0;
  for (int imom = 0; imom < _linmoms; imom++)
    if (filt_moms[imom])
    {
      (out.mom_list)[ifilt] = mom_list[imom];
      (out.p)[ifilt] = p[imom];
      (out.p_tilde)[ifilt] = p_tilde[imom];
      (out.Np)[ifilt] = Np[imom];

      //            (out.linmoms)[ifilt]=linmoms[imom];
      //            (out.bilmoms)[ifilt]=bilmoms[imom];
      (out.linmoms)[ifilt] = {ifilt};
      (out.bilmoms)[ifilt] = {ifilt, ifilt, ifilt};
      if (scheme == "SMOM")
      {
        cout << "Filter not implemented in SMOM." << endl;
        exit(0);
      }

      (out.p2)[ifilt] = p2[imom];
      (out.p2_tilde)[ifilt] = p2_tilde[imom];
      (out.p4)[ifilt] = p4[imom];
      (out.p4_tilde)[ifilt] = p4_tilde[imom];

      ifilt++;
    }

  print_vec(out.p2, path_print + "p2_filtmoms.txt");
  print_vec(out.p2_tilde, path_print + "p2_tilde_filtmoms.txt");

  out.allocate_val();

  out.eff_mass = eff_mass;
  out.eff_mass_sea = eff_mass_sea;

  if (!load_ave)
  {
    out.allocate();

    ifilt = 0;
    for (int imom = 0; imom < _linmoms; imom++)
      if (filt_moms[imom])
      {
        (out.jZq)[ifilt] = jZq[imom];
        ifilt++;
      }

    // out.compute_Zq();

    ifilt = 0;
    for (int imom = 0; imom < _linmoms; imom++)
      if (filt_moms[imom])
      {
        (out.jZ)[ifilt] = jZ[imom];
        ifilt++;
      }

    // out.compute_Zbil();
  }
  return out;
}

oper_t oper_t::average_equiv_moms()
{
  cout << endl;
  cout << "Averaging over the equivalent momenta -- ";

  oper_t out = (*this);

  // Find equivalent linmoms
  int tag = 0, tag_aux = 0;
  double eps = 1.0e-15;

  vector<int> tag_lin_vector;
  tag_lin_vector.push_back(0);

  // Tag assignment to linmoms
  for (int imom = 0; imom < _linmoms; imom++)
  {
    int count_no = 0;

    for (int j = 0; j < imom; j++)
    {
      bool cond{2.0 * abs(p2_tilde[j] - p2_tilde[imom]) < eps * (p2_tilde[j] + p2_tilde[imom]) &&
                mom_list_xyz(mom_list, j) == mom_list_xyz(mom_list, imom) &&
                2.0 * abs(abs(p[j][0]) - abs(p[imom][0])) < eps * (abs(p[j][0]) + abs(p[imom][0]))};

      //            bool cond{2.0*abs(p2[j]-p2[imom])<eps*(p2[j]+p2[imom])};

      if (cond)
        tag_aux = tag_lin_vector[j];
      else
        count_no++;

      if (count_no == imom)
      {
        tag++;
        tag_lin_vector.push_back(tag);
      }
      else if (j == imom - 1)
        tag_lin_vector.push_back(tag_aux);
    }
  }

  // number of equivalent linmoms
  int neq_lin_moms = tag + 1;

  (out.mom_list).resize(neq_lin_moms);
  (out.p).resize(neq_lin_moms);
  (out.p_tilde).resize(neq_lin_moms);
  (out.filt_moms).resize(neq_lin_moms);
  (out.Np).resize(neq_lin_moms);

  (out.p2).resize(neq_lin_moms);
  (out.p2_tilde).resize(neq_lin_moms);
  (out.p4).resize(neq_lin_moms);
  (out.p4_tilde).resize(neq_lin_moms);

  out._linmoms = neq_lin_moms;
  (out.linmoms).resize(out._linmoms);

  cout << "found: " << out._linmoms << " equivalent linmoms ";

  // count the number of momenta for each tag
  vector<int> count_tag_lin_vector(out._linmoms);
  int count = 0;
  for (int tag = 0; tag < out._linmoms; tag++)
  {
    count = 0;
    for (int imom = 0; imom < _linmoms; imom++)
      if (tag_lin_vector[imom] == tag)
        count++;

    count_tag_lin_vector[tag] = count;
  }

  for (int tag = 0; tag < out._linmoms; tag++)
    for (int imom = 0; imom < _linmoms; imom++)
      if (tag_lin_vector[imom] == tag)
      {
        // fill the new linmoms and p2tilde
        (out.mom_list)[tag] = mom_list[imom];
        (out.p)[tag] = p[imom];
        (out.p_tilde)[tag] = p_tilde[imom];

        (out.filt_moms)[tag] = filt_moms[imom];
        (out.Np)[tag] = Np[imom];

        (out.linmoms)[tag] = {tag};
        (out.p2)[tag] = p2[imom];
        (out.p2_tilde)[tag] = p2_tilde[imom];
        (out.p4)[tag] = p4[imom];
        (out.p4_tilde)[tag] = p4_tilde[imom];
      }

  print_vec(out.p2, path_print + "p2_eqmoms.txt");
  print_vec(out.p2_tilde, path_print + "p2_tilde_eqmoms.txt");

  vector<int> tag_bil_vector = tag_lin_vector;
  int neq_bil_moms = neq_lin_moms;

  out._bilmoms = neq_bil_moms;
  cout << "and " << neq_bil_moms << " equivalent bilmoms " << endl << endl;
  (out.bilmoms).resize(out._bilmoms);

  // count the different tags
  vector<int> count_tag_bil_vector(out._bilmoms);
  count = 0;
  for (int tag = 0; tag < out._bilmoms; tag++)
  {
    count = 0;
    for (int imom = 0; imom < _bilmoms; imom++)
      if (tag_bil_vector[imom] == tag)
        count++;

    count_tag_bil_vector[tag] = count;
  }

  for (int tag = 0; tag < out._bilmoms; tag++)
    for (int ibilmom = 0; ibilmom < _bilmoms; ibilmom++)
      if (tag_bil_vector[ibilmom] == tag)
      {
        // fill the new bilmoms
        //                const int imom0=bilmoms[tag][0]; // k
        //                const int imom1=bilmoms[tag][1]; // p1
        //                const int imom2=bilmoms[tag][2]; // p2
        if (scheme == "SMOM")
        {
          cout << "Average not implemented in SMOM." << endl;
          exit(0);
        }

        //                out.bilmoms[tag] = {imom0,imom1,imom2};
        out.bilmoms[tag] = {tag, tag, tag};
      }

  out.allocate_val();

  out.eff_mass = eff_mass;
  out.eff_mass_sea = eff_mass_sea;

  if (!load_ave)
  {
    out.allocate();

    for (int tag = 0; tag < neq_lin_moms; tag++)
      for (int imom = 0; imom < _linmoms; imom++)
        if (tag_lin_vector[imom] == tag)
          for (int ijack = 0; ijack < njacks; ijack++)
            for (int mr = 0; mr < _nmr; mr++)
              for (int iproj = 0; iproj < sigma::nproj; iproj++)
                for (int ins = 0; ins < sigma::nins; ins++)
                  (out.sigma)[tag][iproj][ins][ijack][mr] +=
                      sigma[imom][iproj][ins][ijack][mr] / count_tag_lin_vector[tag];

    out.compute_Zq();

    for (int tag = 0; tag < neq_bil_moms; tag++)
      for (int imom = 0; imom < _bilmoms; imom++)
        if (tag_bil_vector[imom] == tag)
          for (int ibil = 0; ibil < 5; ibil++)
            for (int ijack = 0; ijack < njacks; ijack++)
              for (int ins = 0; ins < gbil::nins; ins++)
                for (int mr1 = 0; mr1 < _nmr; mr1++)
                  for (int mr2 = 0; mr2 < _nmr; mr2++)
                    (out.jG)[tag][ins][ibil][ijack][mr1][mr2] +=
                        jG[imom][ins][ibil][ijack][mr1][mr2] / count_tag_bil_vector[tag];

    out.compute_Zbil();
  }
  return out;
}

oper_t oper_t::average_equiv_momsZ()
{
  cout << endl;
  cout << "Averaging over the equivalent momenta on Z's -- ";

  oper_t out = (*this);

  // Find equivalent linmoms
  int tag = 0, tag_aux = 0;
  double eps = 1.0e-15;

  vector<int> tag_lin_vector;
  tag_lin_vector.push_back(0);

  // Tag assignment to linmoms
  for (int imom = 0; imom < _linmoms; imom++)
  {
    int count_no = 0;

    for (int j = 0; j < imom; j++)
    {
      bool cond{2.0 * abs(p2_tilde[j] - p2_tilde[imom]) < eps * (p2_tilde[j] + p2_tilde[imom]) &&
                mom_list_xyz(mom_list, j) == mom_list_xyz(mom_list, imom) &&
                2.0 * abs(abs(p[j][0]) - abs(p[imom][0])) < eps * (abs(p[j][0]) + abs(p[imom][0]))};

      //            bool cond{2.0*abs(p2[j]-p2[imom])<eps*(p2[j]+p2[imom])};

      if (cond)
        tag_aux = tag_lin_vector[j];
      else
        count_no++;

      if (count_no == imom)
      {
        tag++;
        tag_lin_vector.push_back(tag);
      }
      else if (j == imom - 1)
        tag_lin_vector.push_back(tag_aux);
    }
  }

  // number of equivalent linmoms
  int neq_lin_moms = tag + 1;

  (out.mom_list).resize(neq_lin_moms);
  (out.p).resize(neq_lin_moms);
  (out.p_tilde).resize(neq_lin_moms);
  (out.filt_moms).resize(neq_lin_moms);
  (out.Np).resize(neq_lin_moms);

  (out.p2).resize(neq_lin_moms);
  (out.p2_tilde).resize(neq_lin_moms);
  (out.p4).resize(neq_lin_moms);
  (out.p4_tilde).resize(neq_lin_moms);

  out._linmoms = neq_lin_moms;
  (out.linmoms).resize(out._linmoms);

  cout << "found: " << out._linmoms << " equivalent linmoms ";

  // count the number of momenta for each tag
  vector<int> count_tag_lin_vector(out._linmoms);
  int count = 0;
  for (int tag = 0; tag < out._linmoms; tag++)
  {
    count = 0;
    for (int imom = 0; imom < _linmoms; imom++)
      if (tag_lin_vector[imom] == tag)
        count++;

    count_tag_lin_vector[tag] = count;
  }

  for (int tag = 0; tag < out._linmoms; tag++)
    for (int imom = 0; imom < _linmoms; imom++)
      if (tag_lin_vector[imom] == tag)
      {
        // fill the new linmoms and p2tilde
        (out.mom_list)[tag] = mom_list[imom];
        (out.p)[tag] = p[imom];
        (out.p_tilde)[tag] = p_tilde[imom];

        (out.filt_moms)[tag] = filt_moms[imom];
        (out.Np)[tag] = Np[imom];

        (out.linmoms)[tag] = {tag};
        (out.p2)[tag] = p2[imom];
        (out.p2_tilde)[tag] = p2_tilde[imom];
        (out.p4)[tag] = p4[imom];
        (out.p4_tilde)[tag] = p4_tilde[imom];
      }

  print_vec(out.p2, path_print + "p2_eqmoms.txt");
  print_vec(out.p2_tilde, path_print + "p2_tilde_eqmoms.txt");

  vector<int> tag_bil_vector = tag_lin_vector;
  int neq_bil_moms = neq_lin_moms;

  out._bilmoms = neq_bil_moms;
  cout << "and " << neq_bil_moms << " equivalent bilmoms " << endl << endl;
  (out.bilmoms).resize(out._bilmoms);

  // count the different tags
  vector<int> count_tag_bil_vector(out._bilmoms);
  count = 0;
  for (int tag = 0; tag < out._bilmoms; tag++)
  {
    count = 0;
    for (int imom = 0; imom < _bilmoms; imom++)
      if (tag_bil_vector[imom] == tag)
        count++;

    count_tag_bil_vector[tag] = count;
  }

  for (int tag = 0; tag < out._bilmoms; tag++)
    for (int ibilmom = 0; ibilmom < _bilmoms; ibilmom++)
      if (tag_bil_vector[ibilmom] == tag)
      {
        // fill the new bilmoms
        //                const int imom0=bilmoms[tag][0]; // k
        //                const int imom1=bilmoms[tag][1]; // p1
        //                const int imom2=bilmoms[tag][2]; // p2
        if (scheme == "SMOM")
        {
          cout << "Average not implemented in SMOM." << endl;
          exit(0);
        }

        //                out.bilmoms[tag] = {imom0,imom1,imom2};
        out.bilmoms[tag] = {tag, tag, tag};
      }

  out.allocate_val();

  out.eff_mass = eff_mass;
  out.eff_mass_sea = eff_mass_sea;

  if (!load_ave)
  {
    out.allocate();

    for (int tag = 0; tag < neq_lin_moms; tag++)
      for (int imom = 0; imom < _linmoms; imom++)
        if (tag_lin_vector[imom] == tag)
          for (int ijack = 0; ijack < njacks; ijack++)
            for (int mr = 0; mr < _nmr; mr++)
              (out.jZq)[tag][ijack][mr] += jZq[imom][ijack][mr] / count_tag_lin_vector[tag];

    for (int tag = 0; tag < neq_bil_moms; tag++)
      for (int imom = 0; imom < _bilmoms; imom++)
        if (tag_bil_vector[imom] == tag)
          for (int ibil = 0; ibil < 5; ibil++)
            for (int ijack = 0; ijack < njacks; ijack++)
              for (int mr1 = 0; mr1 < _nmr; mr1++)
                for (int mr2 = 0; mr2 < _nmr; mr2++)
                  (out.jZ)[tag][ibil][ijack][mr1][mr2] +=
                      jZ[imom][ibil][ijack][mr1][mr2] / count_tag_bil_vector[tag];
  }
  return out;
}

oper_t oper_t::pick_democratic()
{
  cout << endl;
  cout << "Picking only spatial-democratic momenta -- ";

  oper_t out = (*this);

  (out.mom_list).clear();
  (out.p).clear();
  (out.p_tilde).clear();
  (out.filt_moms).clear();
  (out.Np).clear();

  (out.p2).clear();
  (out.p2_tilde).clear();
  (out.p4).clear();
  (out.p4_tilde).clear();

  (out.linmoms).clear();
  (out.bilmoms).clear();

  // Filter spatial-democratic momenta with pt > 0 and {px,py,pz} > 0
  int count_dem = 0;
  for (int imom = 0; imom < _linmoms; imom++)
  {
    bool cond{mom_list[imom][1] == mom_list[imom][2] && mom_list[imom][2] == mom_list[imom][3] &&
              mom_list[imom][0] > 0 && mom_list[imom][1] > 0};

    if (cond)
    {
      (out.mom_list).push_back(mom_list[imom]);
      (out.p).push_back(p[imom]);
      (out.p_tilde).push_back(p_tilde[imom]);
      (out.filt_moms).push_back(filt_moms[imom]);
      (out.Np).push_back(Np[imom]);

      (out.p2).push_back(p2[imom]);
      (out.p2_tilde).push_back(p2_tilde[imom]);
      (out.p4).push_back(p4[imom]);
      (out.p4_tilde).push_back(p4_tilde[imom]);

      (out.linmoms).push_back(linmoms[imom]);
      (out.bilmoms).push_back(bilmoms[imom]);

      count_dem++;
    }
  }

  out._linmoms = count_dem;
  out._bilmoms = count_dem;

  cout << "found: " << out._linmoms << " spatial-democratic momenta linmoms ";
  cout << "and " << out._bilmoms << " equivalent bilmoms " << endl << endl;

  print_vec(out.p2, path_print + "p2_dem_moms.txt");
  print_vec(out.p2_tilde, path_print + "p2_tilde_dem_moms.txt");

  if (!load_ave)
  {
    out.allocate();

    int count_dem = 0;
    for (int imom = 0; imom < _linmoms; imom++)
    {
      bool cond{mom_list[imom][1] == mom_list[imom][2] && mom_list[imom][2] == mom_list[imom][3] &&
                mom_list[imom][0] > 0 && mom_list[imom][1] > 0};

      if (cond)
      {
        for (int ijack = 0; ijack < njacks; ijack++)
          for (int mr = 0; mr < _nmr; mr++)
          {
            for (int iproj = 0; iproj < sigma::nproj; iproj++)
              for (int ins = 0; ins < sigma::nins; ins++)
                (out.sigma)[count_dem][iproj][ins][ijack][mr] = sigma[imom][iproj][ins][ijack][mr];

            for (int ibil = 0; ibil < 5; ibil++)
              for (int ins = 0; ins < gbil::nins; ins++)
                for (int mr2 = 0; mr2 < _nmr; mr2++)
                  (out.jG)[count_dem][ins][ibil][ijack][mr][mr2] =
                      jG[imom][ins][ibil][ijack][mr][mr2];
          }

        count_dem++;
      }
    }
  }

  out.compute_Zq();
  out.compute_Zbil();

  return out;
}

void oper_t::print(const string suffix)
{
  print_vec_bin(sigma, path_print + "sigmas_" + suffix);
  print_vec_bin(jG, path_print + "jG_" + suffix);
}

void oper_t::printZ(const string suffix)
{
  print_vec_bin(jZq, path_print + "jZq_" + suffix);
  print_vec_bin(jZ, path_print + "jZ_" + suffix);
}

void oper_t::load(const string suffix)
{
  (*this).allocate();

  ifstream sigma_data(path_print + "sigmas_" + suffix);
  ifstream jG_data(path_print + "jG_" + suffix);

  if (sigma_data.good() and jG_data.good())
  {
    cout << "Loading averaged quantities from files in \"" << path_print << "\"." << endl;
    read_vec_bin(sigma, path_print + "sigmas_" + suffix);
    read_vec_bin(jG, path_print + "jG_" + suffix);
  }
  else
  {
    cout << "Cannot open files of " << suffix << " quantities in \"" << path_print << "\"." << endl;
    exit(0);
  }

  if (suffix.compare("chir") == 1)
  {
    ifstream jZq_data(path_print + "jZq_" + suffix);
    ifstream jZ_data(path_print + "jZ_" + suffix);

    if (jZq_data.good() and jZ_data.good())
    {
      cout << "Loading Z from files in \"" << path_print << "\"." << endl;
      read_vec_bin(jZq, path_print + "jZq_" + suffix);
      read_vec_bin(jZ, path_print + "jZ_" + suffix);
    }
    else
    {
      cout << "Cannot open files of " << suffix << " Z in \"" << path_print << "\"." << endl;
      exit(0);
    }
  }

  (*this).compute_Zq();
  (*this).compute_Zbil();
}

oper_t oper_t::a2p2_extr(int b)
{
  cout << endl;
  cout << "----- [M1]  extrapolation to p2 = 0 -----" << endl << endl;

  oper_t out = (*this);

  out.linmoms = vector<array<int, 1>>{{0}};
  out.bilmoms = vector<array<int, 3>>{{0, 0, 0}};

  out._linmoms = 1;
  out._bilmoms = 1;

  out.allocate_val();
  out.allocate();

  out.eff_mass = (*this).eff_mass;
  out.eff_mass_sea = (*this).eff_mass_sea;

  cout << "p2 range (physical units):   " << p2min << " - " << p2max << endl;

  double ainv2 = ainv[b] * ainv[b];

  int npar = 2;
  vvd_t coord(vd_t(0.0, _linmoms), npar);
  for (int j = 0; j < _linmoms; j++)
  {
    // linear fit in physical units
    coord[0][j] = 1.0;
    coord[1][j] = p2[j] * ainv2; // Gev^2
  }

  vd_t jlincoeff(0.0, njacks), jchisq(0.0, njacks);

  // Interpolating Zq
  vvd_t y_Zq(vd_t(0.0, _linmoms), njacks);           // [njacks][moms]
  vd_t dy_Zq(0.0, _linmoms);                         // [moms]
  vvd_t dy_Zq_tmp = get<1>(ave_err_Zq((*this).jZq)); // [moms][nmr]

  for (int imom = 0; imom < _linmoms; imom++)
  {
    for (int ijack = 0; ijack < njacks; ijack++)
    {
      y_Zq[ijack][imom] = jZq[imom][ijack][0];
    }
    dy_Zq[imom] = dy_Zq_tmp[imom][0];
  }

  vvd_t jZq_pars = polyfit(coord, npar, dy_Zq, y_Zq, p2min, p2max); // [ijack][ipar]

  for (int ijack = 0; ijack < njacks; ijack++)
  {
    (out.jZq)[0][ijack][0] = jZq_pars[ijack][0];
    jlincoeff[ijack] = jZq_pars[ijack][1];
    /**/
    jchisq[ijack] = jZq_pars[ijack][npar];
  }

  cout << "  -- lincoeff[q] = " << get<0>(ave_err(jlincoeff)) << "+/-" << get<1>(ave_err(jlincoeff))
       << endl;
  cout << "    ** chisqr[q] = " << get<0>(ave_err(jchisq)) << "+/-" << get<1>(ave_err(jchisq))
       << endl
       << endl;

  // Interpolating Zbil
  vvd_t y_Zbil(vd_t(0.0, _bilmoms), njacks);           // [njacks][moms]
  vd_t dy_Zbil(0.0, _bilmoms);                         // [moms]
  vvvvd_t dy_Zbil_tmp = get<1>(ave_err_Z((*this).jZ)); // [moms][nbil][nmr][nmr]

  vector<string> str_bil = {"S", "V", "P", "A", "T"};

  for (int ibil = 0; ibil < nbil; ibil++)
  {
    for (int imom = 0; imom < _bilmoms; imom++)
    {
      for (int ijack = 0; ijack < njacks; ijack++)
      {
        y_Zbil[ijack][imom] = jZ[imom][ibil][ijack][0][0];
      }
      dy_Zbil[imom] = dy_Zbil_tmp[imom][ibil][0][0];
    }

    vvd_t jZ_pars = polyfit(coord, npar, dy_Zbil, y_Zbil, p2min, p2max); // [ijack][ipar]

    for (int ijack = 0; ijack < njacks; ijack++)
    {
      (out.jZ)[0][ibil][ijack][0][0] = jZ_pars[ijack][0];
      jlincoeff[ijack] = jZ_pars[ijack][1];
      /**/
      jchisq[ijack] = jZ_pars[ijack][npar];
    }

    cout << "  -- lincoeff[" << str_bil[ibil] << "] = " << get<0>(ave_err(jlincoeff)) << "+/-"
         << get<1>(ave_err(jlincoeff)) << endl;
    cout << "    ** chisqr[" << str_bil[ibil] << "] = " << get<0>(ave_err(jchisq)) << "+/-"
         << get<1>(ave_err(jchisq)) << endl
         << endl;
  }

  // Interpolating ZVoverZA and ZPoverZS and ZAoverZV
  vvd_t y_ZVovZA(vd_t(0.0, _bilmoms), njacks);                  // [njacks][moms]
  vd_t dy_ZVovZA(0.0, _bilmoms);                                // [moms]
  vvvvd_t dy_ZVovZA_tmp = get<1>(ave_err_Z((*this).jZVoverZA)); // [moms][nbil][nmr][nmr]
  vvd_t y_ZPovZS(vd_t(0.0, _bilmoms), njacks);                  // [njacks][moms]
  vd_t dy_ZPovZS(0.0, _bilmoms);                                // [moms]
  vvvvd_t dy_ZPovZS_tmp = get<1>(ave_err_Z((*this).jZPoverZS)); // [moms][nbil][nmr][nmr]
  vvd_t y_ZAovZV(vd_t(0.0, _bilmoms), njacks);                  // [njacks][moms]
  vd_t dy_ZAovZV(0.0, _bilmoms);                                // [moms]
  vvvvd_t dy_ZAovZV_tmp = get<1>(ave_err_Z((*this).jZAoverZV)); // [moms][nbil][nmr][nmr]

  for (int imom = 0; imom < _bilmoms; imom++)
  {
    for (int ijack = 0; ijack < njacks; ijack++)
    {
      y_ZVovZA[ijack][imom] = jZVoverZA[imom][0][ijack][0][0];
      y_ZPovZS[ijack][imom] = jZPoverZS[imom][0][ijack][0][0];
      y_ZAovZV[ijack][imom] = jZAoverZV[imom][0][ijack][0][0];
    }
    dy_ZVovZA[imom] = dy_ZVovZA_tmp[imom][0][0][0];
    dy_ZPovZS[imom] = dy_ZPovZS_tmp[imom][0][0][0];
    dy_ZAovZV[imom] = dy_ZAovZV_tmp[imom][0][0][0];
  }

  vvd_t jZVovZA_pars = polyfit(coord, npar, dy_ZVovZA, y_ZVovZA, p2min, p2max); // [ijack][ipar]
  vvd_t jZPovZS_pars = polyfit(coord, npar, dy_ZPovZS, y_ZPovZS, p2min, p2max); // [ijack][ipar]
  vvd_t jZAovZV_pars = polyfit(coord, npar, dy_ZAovZV, y_ZAovZV, p2min, p2max); // [ijack][ipar]

  for (int ijack = 0; ijack < njacks; ijack++)
  {
    (out.jZVoverZA)[0][0][ijack][0][0] = jZVovZA_pars[ijack][0];
    (out.jZPoverZS)[0][0][ijack][0][0] = jZPovZS_pars[ijack][0];
    (out.jZAoverZV)[0][0][ijack][0][0] = jZAovZV_pars[ijack][0];
  }

  return out;
}

oper_t oper_t::a2p2_extr_with_pole(int b)
{
  cout << endl;
  cout << "----- [M3] extrapolation to p2 = 0 [with pole 1/p2] -----" << endl << endl;

  oper_t out = (*this);

  out.linmoms = vector<array<int, 1>>{{0}};
  out.bilmoms = vector<array<int, 3>>{{0, 0, 0}};

  out._linmoms = 1;
  out._bilmoms = 1;

  out.allocate_val();
  out.allocate();

  out.eff_mass = (*this).eff_mass;
  out.eff_mass_sea = (*this).eff_mass_sea;

  double ainv2 = ainv[b] * ainv[b];

  double _p2min = 0; // GeV^2
  double _p2max = 0;
  if (nbeta == 1)
  {
    _p2min = p2min_M3_M4[0];
    _p2max = p2max_M3_M4[0];
  }
  else if (nbeta == 3)
  {
    if (fabs(ainv[b] - 2.1218) < 1e-10)
    {
      _p2min = p2min_M3_M4[0];
      _p2max = p2max_M3_M4[0];
    }
    else if (fabs(ainv[b] - 2.4635) < 1e-10)
    {
      _p2min = p2min_M3_M4[1];
      _p2max = p2max_M3_M4[1];
    }
    else if (fabs(ainv[b] - 2.8934) < 1e-10)
    {
      _p2min = p2min_M3_M4[2];
      _p2max = p2max_M3_M4[2];
    }
  }

  cout << "p2 range (physical units):   " << _p2min << " - " << _p2max << endl;

  int npar = 3;
  vvd_t coord(vd_t(0.0, _linmoms), npar);
  for (int j = 0; j < _linmoms; j++)
  {
    // linear fit in physical units
    coord[0][j] = 1.0;
    coord[1][j] = p2[j] * ainv2;       // GeV^2
    coord[2][j] = 1.0 / (coord[1][j]); // 1/GeV^2
  }

  vd_t jpole(0.0, njacks), jlincoeff(0.0, njacks), jchisq(0.0, njacks);

  // Interpolating Zq
  vvd_t y_Zq(vd_t(0.0, _linmoms), njacks);           // [njacks][moms]
  vd_t dy_Zq(0.0, _linmoms);                         // [moms]
  vvd_t dy_Zq_tmp = get<1>(ave_err_Zq((*this).jZq)); // [moms][nmr]

  for (int imom = 0; imom < _linmoms; imom++)
  {
    for (int ijack = 0; ijack < njacks; ijack++)
    {
      y_Zq[ijack][imom] = jZq[imom][ijack][0];
    }
    dy_Zq[imom] = dy_Zq_tmp[imom][0];
  }

  vvd_t jZq_pars = polyfit(coord, npar, dy_Zq, y_Zq, _p2min, _p2max); // [ijack][ipar]

  for (int ijack = 0; ijack < njacks; ijack++)
  {
    (out.jZq)[0][ijack][0] = jZq_pars[ijack][0];
    jlincoeff[ijack] = jZq_pars[ijack][1];
    jpole[ijack] = jZq_pars[ijack][2];
    /**/
    jchisq[ijack] = jZq_pars[ijack][npar];
  }

  cout << "  -- pole[q] = " << get<0>(ave_err(jpole)) << "+/-" << get<1>(ave_err(jpole)) << endl;
  cout << "  -- lincoeff[q] = " << get<0>(ave_err(jlincoeff)) << "+/-" << get<1>(ave_err(jlincoeff))
       << endl;
  cout << "    ** chisqr[q] = " << get<0>(ave_err(jchisq)) << "+/-" << get<1>(ave_err(jchisq))
       << endl
       << endl;

  // Interpolating Zbil
  vvd_t y_Zbil(vd_t(0.0, _bilmoms), njacks);           // [njacks][moms]
  vd_t dy_Zbil(0.0, _bilmoms);                         // [moms]
  vvvvd_t dy_Zbil_tmp = get<1>(ave_err_Z((*this).jZ)); // [moms][nbil][nmr][nmr]

  vector<double> p2min_bil = {_p2min, _p2min, _p2min, _p2min, _p2min}; /* {S,V,P,A,T} */
  vector<string> str_bil = {"S", "V", "P", "A", "T"};

  for (int ibil = 0; ibil < nbil; ibil++)
  {
    for (int imom = 0; imom < _bilmoms; imom++)
    {
      for (int ijack = 0; ijack < njacks; ijack++)
      {
        y_Zbil[ijack][imom] = jZ[imom][ibil][ijack][0][0];
      }
      dy_Zbil[imom] = dy_Zbil_tmp[imom][ibil][0][0];
    }

    vvd_t jZ_pars = polyfit(coord, npar, dy_Zbil, y_Zbil, p2min_bil[ibil], _p2max); // [ijack][ipar]

    for (int ijack = 0; ijack < njacks; ijack++)
    {
      (out.jZ)[0][ibil][ijack][0][0] = jZ_pars[ijack][0];
      jlincoeff[ijack] = jZ_pars[ijack][1];
      jpole[ijack] = jZ_pars[ijack][2];
      /**/
      jchisq[ijack] = jZ_pars[ijack][npar];
    }

    cout << "  -- pole[" << str_bil[ibil] << "] = " << get<0>(ave_err(jpole)) << "+/-"
         << get<1>(ave_err(jpole)) << endl;
    cout << "  -- lincoeff[" << str_bil[ibil] << "] = " << get<0>(ave_err(jlincoeff)) << "+/-"
         << get<1>(ave_err(jlincoeff)) << endl;
    cout << "    ** chisqr[" << str_bil[ibil] << "] = " << get<0>(ave_err(jchisq)) << "+/-"
         << get<1>(ave_err(jchisq)) << endl
         << endl;
  }

  // Interpolating ZVoverZA and ZPoverZS
  vvd_t y_ZVovZA(vd_t(0.0, _bilmoms), njacks);                  // [njacks][moms]
  vd_t dy_ZVovZA(0.0, _bilmoms);                                // [moms]
  vvvvd_t dy_ZVovZA_tmp = get<1>(ave_err_Z((*this).jZVoverZA)); // [moms][nbil][nmr][nmr]
  vvd_t y_ZPovZS(vd_t(0.0, _bilmoms), njacks);                  // [njacks][moms]
  vd_t dy_ZPovZS(0.0, _bilmoms);                                // [moms]
  vvvvd_t dy_ZPovZS_tmp = get<1>(ave_err_Z((*this).jZPoverZS)); // [moms][nbil][nmr][nmr]
  vvd_t y_ZAovZV(vd_t(0.0, _bilmoms), njacks);                  // [njacks][moms]
  vd_t dy_ZAovZV(0.0, _bilmoms);                                // [moms]
  vvvvd_t dy_ZAovZV_tmp = get<1>(ave_err_Z((*this).jZAoverZV)); // [moms][nbil][nmr][nmr]

  for (int imom = 0; imom < _bilmoms; imom++)
  {
    for (int ijack = 0; ijack < njacks; ijack++)
    {
      y_ZVovZA[ijack][imom] = jZVoverZA[imom][0][ijack][0][0];
      y_ZPovZS[ijack][imom] = jZPoverZS[imom][0][ijack][0][0];
      y_ZAovZV[ijack][imom] = jZAoverZV[imom][0][ijack][0][0];
    }
    dy_ZVovZA[imom] = dy_ZVovZA_tmp[imom][0][0][0];
    dy_ZPovZS[imom] = dy_ZPovZS_tmp[imom][0][0][0];
    dy_ZAovZV[imom] = dy_ZAovZV_tmp[imom][0][0][0];
  }

  vvd_t jZVovZA_pars = polyfit(coord, npar, dy_ZVovZA, y_ZVovZA, _p2min, _p2max); // [ijack][ipar]
  vvd_t jZPovZS_pars = polyfit(coord, npar, dy_ZPovZS, y_ZPovZS, _p2min, _p2max); // [ijack][ipar]
  vvd_t jZAovZV_pars = polyfit(coord, npar, dy_ZAovZV, y_ZAovZV, _p2min, _p2max); // [ijack][ipar]

  for (int ijack = 0; ijack < njacks; ijack++)
  {
    (out.jZVoverZA)[0][0][ijack][0][0] = jZVovZA_pars[ijack][0];
    (out.jZPoverZS)[0][0][ijack][0][0] = jZPovZS_pars[ijack][0];
    (out.jZAoverZV)[0][0][ijack][0][0] = jZAovZV_pars[ijack][0];

    // parameters for ZA/ZV
    jlincoeff[ijack] = jZAovZV_pars[ijack][1];
    jpole[ijack] = jZAovZV_pars[ijack][2];
    /**/
    jchisq[ijack] = jZAovZV_pars[ijack][npar];
  }

  cout << "  -- pole[A/V] = " << get<0>(ave_err(jpole)) << "+/-" << get<1>(ave_err(jpole)) << endl;
  cout << "  -- lincoeff[A/V] = " << get<0>(ave_err(jlincoeff)) << "+/-"
       << get<1>(ave_err(jlincoeff)) << endl;
  cout << "    ** chisqr[A/V] = " << get<0>(ave_err(jchisq)) << "+/-" << get<1>(ave_err(jchisq))
       << endl
       << endl;

  return out;
}

oper_t oper_t::a2p2_extr_with_pole_and_p4(int b)
{
  cout << endl;
  cout << "----- [M4] extrapolation to p2 = 0 [with pole 1/p2 + quadratic] -----" << endl << endl;

  oper_t out = (*this);

  out.linmoms = vector<array<int, 1>>{{0}};
  out.bilmoms = vector<array<int, 3>>{{0, 0, 0}};

  out._linmoms = 1;
  out._bilmoms = 1;

  out.allocate_val();
  out.allocate();

  out.eff_mass = (*this).eff_mass;
  out.eff_mass_sea = (*this).eff_mass_sea;

  double ainv2 = ainv[b] * ainv[b];

  double _p2min = 0; // GeV^2
  double _p2max = 0;
  if (fabs(ainv[b] - 2.1218) < 1e-10)
  {
    _p2min = p2min_M3_M4[0];
    _p2max = p2max_M3_M4[0];
  }
  else if (fabs(ainv[b] - 2.4635) < 1e-10)
  {
    _p2min = p2min_M3_M4[1];
    _p2max = p2max_M3_M4[1];
  }
  else if (fabs(ainv[b] - 2.8934) < 1e-10)
  {
    _p2min = p2min_M3_M4[2];
    _p2max = p2max_M3_M4[2];
  }

  cout << "p2 range (physical units):   " << _p2min << " - " << _p2max << endl;

  int npar = 4;
  vvd_t coord(vd_t(0.0, _linmoms), npar);
  for (int j = 0; j < _linmoms; j++)
  {
    // linear fit in physical units
    coord[0][j] = 1.0;
    coord[1][j] = p2[j] * ainv2;             // GeV^2
    coord[2][j] = 1.0 / (coord[1][j]);       // 1/GeV^2
    coord[3][j] = coord[1][j] * coord[1][j]; // GeV^4
  }

  vd_t jpole(0.0, njacks), jlincoeff(0.0, njacks), jsqrcoeff(0.0, njacks), jchisq(0.0, njacks);

  // Interpolating Zq
  vvd_t y_Zq(vd_t(0.0, _linmoms), njacks);           // [njacks][moms]
  vd_t dy_Zq(0.0, _linmoms);                         // [moms]
  vvd_t dy_Zq_tmp = get<1>(ave_err_Zq((*this).jZq)); // [moms][nmr]

  for (int imom = 0; imom < _linmoms; imom++)
  {
    for (int ijack = 0; ijack < njacks; ijack++)
    {
      y_Zq[ijack][imom] = jZq[imom][ijack][0];
    }
    dy_Zq[imom] = dy_Zq_tmp[imom][0];
  }

  vvd_t jZq_pars = polyfit(coord, npar, dy_Zq, y_Zq, _p2min, _p2max); // [ijack][ipar]

  for (int ijack = 0; ijack < njacks; ijack++)
  {
    (out.jZq)[0][ijack][0] = jZq_pars[ijack][0];
    jlincoeff[ijack] = jZq_pars[ijack][1];
    jpole[ijack] = jZq_pars[ijack][2];
    jsqrcoeff[ijack] = jZq_pars[ijack][3];
    /**/
    jchisq[ijack] = jZq_pars[ijack][npar];
  }

  cout << "  -- pole[q] = " << get<0>(ave_err(jpole)) << "+/-" << get<1>(ave_err(jpole)) << endl;
  cout << "  -- lincoeff[q] = " << get<0>(ave_err(jlincoeff)) << "+/-" << get<1>(ave_err(jlincoeff))
       << endl;
  cout << "  -- sqrcoeff[q] = " << get<0>(ave_err(jsqrcoeff)) << "+/-" << get<1>(ave_err(jsqrcoeff))
       << endl;
  cout << "    ** chisqr[q] = " << get<0>(ave_err(jchisq)) << "+/-" << get<1>(ave_err(jchisq))
       << endl
       << endl;

  // Interpolating Zbil
  vvd_t y_Zbil(vd_t(0.0, _bilmoms), njacks);           // [njacks][moms]
  vd_t dy_Zbil(0.0, _bilmoms);                         // [moms]
  vvvvd_t dy_Zbil_tmp = get<1>(ave_err_Z((*this).jZ)); // [moms][nbil][nmr][nmr]

  vector<double> p2min_bil = {_p2min, _p2min, _p2min, _p2min, _p2min}; /* {S,V,P,A,T} */
  vector<string> str_bil = {"S", "V", "P", "A", "T"};

  for (int ibil = 0; ibil < nbil; ibil++)
  {
    for (int imom = 0; imom < _bilmoms; imom++)
    {
      for (int ijack = 0; ijack < njacks; ijack++)
      {
        y_Zbil[ijack][imom] = jZ[imom][ibil][ijack][0][0];
      }
      dy_Zbil[imom] = dy_Zbil_tmp[imom][ibil][0][0];
    }

    vvd_t jZ_pars = polyfit(coord, npar, dy_Zbil, y_Zbil, p2min_bil[ibil], _p2max); // [ijack][ipar]

    for (int ijack = 0; ijack < njacks; ijack++)
    {
      (out.jZ)[0][ibil][ijack][0][0] = jZ_pars[ijack][0];
      jlincoeff[ijack] = jZ_pars[ijack][1];
      jpole[ijack] = jZ_pars[ijack][2];
      jsqrcoeff[ijack] = jZ_pars[ijack][3];
      /**/
      jchisq[ijack] = jZ_pars[ijack][npar];
    }

    cout << "  -- pole[" << str_bil[ibil] << "] = " << get<0>(ave_err(jpole)) << "+/-"
         << get<1>(ave_err(jpole)) << endl;
    cout << "  -- lincoeff[" << str_bil[ibil] << "] = " << get<0>(ave_err(jlincoeff)) << "+/-"
         << get<1>(ave_err(jlincoeff)) << endl;
    cout << "  -- sqrcoeff[" << str_bil[ibil] << "] = " << get<0>(ave_err(jsqrcoeff)) << "+/-"
         << get<1>(ave_err(jsqrcoeff)) << endl;
    cout << "    ** chisqr[" << str_bil[ibil] << "] = " << get<0>(ave_err(jchisq)) << "+/-"
         << get<1>(ave_err(jchisq)) << endl
         << endl;
  }

  // Interpolating ZVoverZA and ZPoverZS
  vvd_t y_ZVovZA(vd_t(0.0, _bilmoms), njacks);                  // [njacks][moms]
  vd_t dy_ZVovZA(0.0, _bilmoms);                                // [moms]
  vvvvd_t dy_ZVovZA_tmp = get<1>(ave_err_Z((*this).jZVoverZA)); // [moms][nbil][nmr][nmr]
  vvd_t y_ZPovZS(vd_t(0.0, _bilmoms), njacks);                  // [njacks][moms]
  vd_t dy_ZPovZS(0.0, _bilmoms);                                // [moms]
  vvvvd_t dy_ZPovZS_tmp = get<1>(ave_err_Z((*this).jZPoverZS)); // [moms][nbil][nmr][nmr]
  vvd_t y_ZAovZV(vd_t(0.0, _bilmoms), njacks);                  // [njacks][moms]
  vd_t dy_ZAovZV(0.0, _bilmoms);                                // [moms]
  vvvvd_t dy_ZAovZV_tmp = get<1>(ave_err_Z((*this).jZAoverZV)); // [moms][nbil][nmr][nmr]

  for (int imom = 0; imom < _bilmoms; imom++)
  {
    for (int ijack = 0; ijack < njacks; ijack++)
    {
      y_ZVovZA[ijack][imom] = jZVoverZA[imom][0][ijack][0][0];
      y_ZPovZS[ijack][imom] = jZPoverZS[imom][0][ijack][0][0];
      y_ZAovZV[ijack][imom] = jZAoverZV[imom][0][ijack][0][0];
    }
    dy_ZVovZA[imom] = dy_ZVovZA_tmp[imom][0][0][0];
    dy_ZPovZS[imom] = dy_ZPovZS_tmp[imom][0][0][0];
    dy_ZAovZV[imom] = dy_ZAovZV_tmp[imom][0][0][0];
  }

  vvd_t jZVovZA_pars = polyfit(coord, npar, dy_ZVovZA, y_ZVovZA, _p2min, _p2max); // [ijack][ipar]
  vvd_t jZPovZS_pars = polyfit(coord, npar, dy_ZPovZS, y_ZPovZS, _p2min, _p2max); // [ijack][ipar]
  vvd_t jZAovZV_pars = polyfit(coord, npar, dy_ZAovZV, y_ZAovZV, _p2min, _p2max); // [ijack][ipar]

  for (int ijack = 0; ijack < njacks; ijack++)
  {
    (out.jZVoverZA)[0][0][ijack][0][0] = jZVovZA_pars[ijack][0];
    (out.jZPoverZS)[0][0][ijack][0][0] = jZPovZS_pars[ijack][0];
    (out.jZAoverZV)[0][0][ijack][0][0] = jZAovZV_pars[ijack][0];
  }

  return out;
}

oper_t oper_t::a2p2_extr_with_p4(int b)
{
  cout << endl;
  cout << "----- extrapolation to p2 = 0 [with quadratic] -----" << endl << endl;

  oper_t out = (*this);

  out.linmoms = vector<array<int, 1>>{{0}};
  out.bilmoms = vector<array<int, 3>>{{0, 0, 0}};

  out._linmoms = 1;
  out._bilmoms = 1;

  out.allocate_val();
  out.allocate();

  out.eff_mass = (*this).eff_mass;
  out.eff_mass_sea = (*this).eff_mass_sea;

  double ainv2 = ainv[b] * ainv[b];

  double _p2min = 0; // GeV^2
  double _p2max = 0;
  if (fabs(ainv[b] - 2.1218) < 1e-10)
  {
    _p2min = p2min_M3_M4[0];
    _p2max = p2max_M3_M4[0];
  }
  else if (fabs(ainv[b] - 2.4635) < 1e-10)
  {
    _p2min = p2min_M3_M4[1];
    _p2max = p2max_M3_M4[1];
  }
  else if (fabs(ainv[b] - 2.8934) < 1e-10)
  {
    _p2min = p2min_M3_M4[2];
    _p2max = p2max_M3_M4[2];
  }

  cout << "p2 range (physical units):   " << _p2min << " - " << _p2max << endl;

  int npar = 3;
  vvd_t coord(vd_t(0.0, _linmoms), npar);
  for (int j = 0; j < _linmoms; j++)
  {
    // linear fit in physical units
    coord[0][j] = 1.0;
    coord[1][j] = p2[j] * ainv2;             // GeV^2
    coord[2][j] = coord[1][j] * coord[1][j]; // GeV^4
  }

  vd_t jpole(0.0, njacks), jlincoeff(0.0, njacks), jsqrcoeff(0.0, njacks), jchisq(0.0, njacks);

  // Interpolating Zq
  vvd_t y_Zq(vd_t(0.0, _linmoms), njacks);           // [njacks][moms]
  vd_t dy_Zq(0.0, _linmoms);                         // [moms]
  vvd_t dy_Zq_tmp = get<1>(ave_err_Zq((*this).jZq)); // [moms][nmr]

  for (int imom = 0; imom < _linmoms; imom++)
  {
    for (int ijack = 0; ijack < njacks; ijack++)
    {
      y_Zq[ijack][imom] = jZq[imom][ijack][0];
    }
    dy_Zq[imom] = dy_Zq_tmp[imom][0];
  }

  vvd_t jZq_pars = polyfit(coord, npar, dy_Zq, y_Zq, _p2min, _p2max); // [ijack][ipar]

  for (int ijack = 0; ijack < njacks; ijack++)
  {
    (out.jZq)[0][ijack][0] = jZq_pars[ijack][0];
    jlincoeff[ijack] = jZq_pars[ijack][1];
    jsqrcoeff[ijack] = jZq_pars[ijack][2];
    /**/
    jchisq[ijack] = jZq_pars[ijack][npar];
  }

  cout << "  -- lincoeff[q] = " << get<0>(ave_err(jlincoeff)) << "+/-" << get<1>(ave_err(jlincoeff))
       << endl;
  cout << "  -- sqrcoeff[q] = " << get<0>(ave_err(jsqrcoeff)) << "+/-" << get<1>(ave_err(jsqrcoeff))
       << endl;
  cout << "    ** chisqr[q] = " << get<0>(ave_err(jchisq)) << "+/-" << get<1>(ave_err(jchisq))
       << endl
       << endl;

  // Interpolating Zbil
  vvd_t y_Zbil(vd_t(0.0, _bilmoms), njacks);           // [njacks][moms]
  vd_t dy_Zbil(0.0, _bilmoms);                         // [moms]
  vvvvd_t dy_Zbil_tmp = get<1>(ave_err_Z((*this).jZ)); // [moms][nbil][nmr][nmr]

  vector<double> p2min_bil = {_p2min, _p2min, _p2min, _p2min, _p2min}; /* {S,V,P,A,T} */
  vector<string> str_bil = {"S", "V", "P", "A", "T"};

  for (int ibil = 0; ibil < nbil; ibil++)
  {
    for (int imom = 0; imom < _bilmoms; imom++)
    {
      for (int ijack = 0; ijack < njacks; ijack++)
      {
        y_Zbil[ijack][imom] = jZ[imom][ibil][ijack][0][0];
      }
      dy_Zbil[imom] = dy_Zbil_tmp[imom][ibil][0][0];
    }

    vvd_t jZ_pars = polyfit(coord, npar, dy_Zbil, y_Zbil, p2min_bil[ibil], _p2max); // [ijack][ipar]

    for (int ijack = 0; ijack < njacks; ijack++)
    {
      (out.jZ)[0][ibil][ijack][0][0] = jZ_pars[ijack][0];
      jlincoeff[ijack] = jZ_pars[ijack][1];
      jsqrcoeff[ijack] = jZ_pars[ijack][2];
      /**/
      jchisq[ijack] = jZ_pars[ijack][npar];
    }

    cout << "  -- lincoeff[" << str_bil[ibil] << "] = " << get<0>(ave_err(jlincoeff)) << "+/-"
         << get<1>(ave_err(jlincoeff)) << endl;
    cout << "  -- sqrcoeff[" << str_bil[ibil] << "] = " << get<0>(ave_err(jsqrcoeff)) << "+/-"
         << get<1>(ave_err(jsqrcoeff)) << endl;
    cout << "    ** chisqr[" << str_bil[ibil] << "] = " << get<0>(ave_err(jchisq)) << "+/-"
         << get<1>(ave_err(jchisq)) << endl
         << endl;
  }

  // Interpolating ZVoverZA and ZPoverZS
  vvd_t y_ZVovZA(vd_t(0.0, _bilmoms), njacks);                  // [njacks][moms]
  vd_t dy_ZVovZA(0.0, _bilmoms);                                // [moms]
  vvvvd_t dy_ZVovZA_tmp = get<1>(ave_err_Z((*this).jZVoverZA)); // [moms][nbil][nmr][nmr]
  vvd_t y_ZPovZS(vd_t(0.0, _bilmoms), njacks);                  // [njacks][moms]
  vd_t dy_ZPovZS(0.0, _bilmoms);                                // [moms]
  vvvvd_t dy_ZPovZS_tmp = get<1>(ave_err_Z((*this).jZPoverZS)); // [moms][nbil][nmr][nmr]
  vvd_t y_ZAovZV(vd_t(0.0, _bilmoms), njacks);                  // [njacks][moms]
  vd_t dy_ZAovZV(0.0, _bilmoms);                                // [moms]
  vvvvd_t dy_ZAovZV_tmp = get<1>(ave_err_Z((*this).jZAoverZV)); // [moms][nbil][nmr][nmr]

  for (int imom = 0; imom < _bilmoms; imom++)
  {
    for (int ijack = 0; ijack < njacks; ijack++)
    {
      y_ZVovZA[ijack][imom] = jZVoverZA[imom][0][ijack][0][0];
      y_ZPovZS[ijack][imom] = jZPoverZS[imom][0][ijack][0][0];
      y_ZAovZV[ijack][imom] = jZAoverZV[imom][0][ijack][0][0];
    }
    dy_ZVovZA[imom] = dy_ZVovZA_tmp[imom][0][0][0];
    dy_ZPovZS[imom] = dy_ZPovZS_tmp[imom][0][0][0];
    dy_ZAovZV[imom] = dy_ZAovZV_tmp[imom][0][0][0];
  }

  vvd_t jZVovZA_pars = polyfit(coord, npar, dy_ZVovZA, y_ZVovZA, _p2min, _p2max); // [ijack][ipar]
  vvd_t jZPovZS_pars = polyfit(coord, npar, dy_ZPovZS, y_ZPovZS, _p2min, _p2max); // [ijack][ipar]
  vvd_t jZAovZV_pars = polyfit(coord, npar, dy_ZAovZV, y_ZAovZV, _p2min, _p2max); // [ijack][ipar]

  for (int ijack = 0; ijack < njacks; ijack++)
  {
    (out.jZVoverZA)[0][0][ijack][0][0] = jZVovZA_pars[ijack][0];
    (out.jZPoverZS)[0][0][ijack][0][0] = jZPovZS_pars[ijack][0];
    (out.jZAoverZV)[0][0][ijack][0][0] = jZAovZV_pars[ijack][0];
  }

  return out;
}

oper_t oper_t::remove_hadr_cont(double _ainv)
{
  cout << endl;
  cout << "----- removing hadr. cont. ~ 1/p^2 by hand -----" << endl << endl;

  oper_t out = (*this);

  double epsq = 0.0;
  vector<double> eps = {0.0, 0.24, 0.0, 0.226, 0.0}; //{S,V,P,A,T}  [GeV^2]
  double epsVA = 0.0;
  double epsAV = 0.0;
  double epsPS = 0.497;

  for (int imom = 0; imom < _bilmoms; imom++)
    for (int ijack = 0; ijack < njacks; ijack++)
    {
      (out.jZq)[imom][ijack][0] -= epsq / (p2[imom] * _ainv * _ainv);
      for (int ibil = 0; ibil < nbil; ibil++)
        (out.jZ)[imom][ibil][ijack][0][0] -= eps[ibil] / (p2[imom] * _ainv * _ainv);
      (out.jZVoverZA)[imom][0][ijack][0][0] -= epsVA / (p2[imom] * _ainv * _ainv);
      (out.jZPoverZS)[imom][0][ijack][0][0] -= epsPS / (p2[imom] * _ainv * _ainv);
      (out.jZAoverZV)[imom][0][ijack][0][0] -= epsAV / (p2[imom] * _ainv * _ainv);
    }

  return out;
}

oper_t oper_t::Z_improvement(double _ainv)
{
  cout << endl;
  cout << "----- [strategy: a1] Z-improvement -----" << endl << endl;

  oper_t out = (*this);

  out.linmoms = vector<array<int, 1>>{{0}};
  out.bilmoms = vector<array<int, 3>>{{0, 0, 0}};

  out._linmoms = 1;
  out._bilmoms = 1;

  out.allocate_val();
  out.allocate();

  out.eff_mass = (*this).eff_mass;
  out.eff_mass_sea = (*this).eff_mass_sea;

  // cout<<"p2 range (physical units):   "<<p2min<<" - "<<p2max<<endl;

  int npar = 2;
  double a2p2minn = 15.0 / (_ainv * _ainv); // in lattice units!
  double a2p2maxx = 30.0 / (_ainv * _ainv);

  int length = 0;
  for (int imom = 0; imom < _linmoms; imom++)
    if (p2[imom] > a2p2minn && p2[imom] < a2p2maxx)
      for (int jmom = imom + 1; jmom < _linmoms; jmom++)
        if (/*p2[jmom]>a2p2minn*/ (p2[jmom] - p2[imom]) > 1.0 / (_ainv * _ainv) &&
            p2[jmom] < a2p2maxx)
          length++;

  vvd_t coord(vd_t(0.0, length), npar);

  // Interpolating Z's
  vvd_t y_Zq(vd_t(0.0, length), njacks);                 // [njacks][moms]
  vvvd_t y_Zbil(vvd_t(vd_t(0.0, length), njacks), nbil); // [njacks][moms]
  vvd_t y_ZVovZA(vd_t(0.0, length), njacks);             // [njacks][moms]
  vvd_t y_ZPovZS(vd_t(0.0, length), njacks);             // [njacks][moms]
  vvd_t y_ZAovZV(vd_t(0.0, length), njacks);             // [njacks][moms]

  int l = 0;
  for (int imom = 0; imom < _linmoms; imom++)
    if (p2[imom] > a2p2minn && p2[imom] < a2p2maxx)
      for (int jmom = imom + 1; jmom < _linmoms; jmom++)
        if (/*p2[jmom]>a2p2minn*/ (p2[jmom] - p2[imom]) > 1.0 / (_ainv * _ainv) &&
            p2[jmom] < a2p2maxx)
        {
          // linear fit in physical units
          coord[0][l] = 1.0;
          coord[1][l] = (p2[imom] + p2[jmom]);

          for (int ijack = 0; ijack < njacks; ijack++)
          {
            y_Zq[ijack][l] = (jZq[jmom][ijack][0] * p2[jmom] - jZq[imom][ijack][0] * p2[imom]) /
                             (p2[jmom] - p2[imom]);
            for (int ibil = 0; ibil < nbil; ibil++)
              y_Zbil[ibil][ijack][l] = (jZ[jmom][ibil][ijack][0][0] * p2[jmom] -
                                        jZ[imom][ibil][ijack][0][0] * p2[imom]) /
                                       (p2[jmom] - p2[imom]);
            y_ZVovZA[ijack][l] = (jZVoverZA[jmom][0][ijack][0][0] * p2[jmom] -
                                  jZVoverZA[imom][0][ijack][0][0] * p2[imom]) /
                                 (p2[jmom] - p2[imom]);
            y_ZPovZS[ijack][l] = (jZPoverZS[jmom][0][ijack][0][0] * p2[jmom] -
                                  jZPoverZS[imom][0][ijack][0][0] * p2[imom]) /
                                 (p2[jmom] - p2[imom]);
            y_ZAovZV[ijack][l] = (jZAoverZV[jmom][0][ijack][0][0] * p2[jmom] -
                                  jZAoverZV[imom][0][ijack][0][0] * p2[imom]) /
                                 (p2[jmom] - p2[imom]);
          }
          l++;
        }

  vd_t dy_Zq = get<1>(ave_err(y_Zq));
  vvd_t jZq_pars = polyfit(coord, npar, dy_Zq, y_Zq, 0.0, 100.0 /*p2min,p2max*/); // [ijack][ipar]

  vd_t dy_ZVovZA = get<1>(ave_err(y_ZVovZA));
  vd_t dy_ZPovZS = get<1>(ave_err(y_ZPovZS));
  vd_t dy_ZAovZV = get<1>(ave_err(y_ZAovZV));
  vvd_t jZVovZA_pars =
      polyfit(coord, npar, dy_ZVovZA, y_ZVovZA, 0.0, 100.0 /*p2min,p2max*/); // [ijack][ipar]
  vvd_t jZPovZS_pars =
      polyfit(coord, npar, dy_ZPovZS, y_ZPovZS, 0.0, 100.0 /*p2min,p2max*/); // [ijack][ipar]
  vvd_t jZAovZV_pars =
      polyfit(coord, npar, dy_ZAovZV, y_ZAovZV, 0.0, 100.0 /*p2min,p2max*/); // [ijack][ipar]

  for (int ijack = 0; ijack < njacks; ijack++)
  {
    (out.jZq)[0][ijack][0] = jZq_pars[ijack][0];
    (out.jZVoverZA)[0][0][ijack][0][0] = jZVovZA_pars[ijack][0];
    (out.jZPoverZS)[0][0][ijack][0][0] = jZPovZS_pars[ijack][0];
    (out.jZAoverZV)[0][0][ijack][0][0] = jZAovZV_pars[ijack][0];
  }

  for (int ibil = 0; ibil < nbil; ibil++)
  {
    vd_t dy_Zbil = get<1>(ave_err(y_Zbil[ibil]));
    vvd_t jZ_pars =
        polyfit(coord, npar, dy_Zbil, y_Zbil[ibil], 0.0, 100.0 /*p2min,p2max*/); // [ijack][ipar]

    for (int ijack = 0; ijack < njacks; ijack++)
    {
      (out.jZ)[0][ibil][ijack][0][0] = jZ_pars[ijack][0];
    }
  }

  return out;
}

voper_t combined_M3(voper_t in) // M3 method combined on all betas
{

  cout << endl;
  cout << "----- Combined M3 fit on all the curves -----" << endl << endl;

  voper_t out(nbeta);

  for (int b = 0; b < nbeta; b++)
  {
    out[b] = in[b];

    out[b].linmoms = vector<array<int, 1>>{{0}};
    out[b].bilmoms = vector<array<int, 3>>{{0, 0, 0}};

    out[b]._linmoms = 1;
    out[b]._bilmoms = 1;

    out[b].allocate_val();
    out[b].allocate();
  }

  double _p2min[] = {p2min_M3_M4[0], p2min_M3_M4[1], p2min_M3_M4[2]}; // GeV^2
  double _p2max[] = {p2max_M3_M4[0], p2max_M3_M4[1], p2max_M3_M4[2]};

  // count all the momenta respecting the above criteria of p2min&p2max
  int _linmoms_tot = 0;
  for (int b = 0; b < nbeta; b++)
  {
    double ainv2 = ainv[b] * ainv[b];
    for (int j = 0; j < in[b]._linmoms; j++)
      if (in[b].p2[j] * ainv2 >= _p2min[b] && in[b].p2[j] * ainv2 <= _p2max[b])
        _linmoms_tot++;
  }

  int _linmoms_tot_complete = 0;
  for (int b = 0; b < nbeta; b++)
  {
    _linmoms_tot_complete += in[b]._linmoms;
  }
  cout << "tot moms = " << _linmoms_tot << "/" << _linmoms_tot_complete << endl;

  int npar = nbeta + 2; // {ZA,ZB,ZC}+{linear}+{pole}

  // for(int b=0; b<nbeta;b++)
  //   _linmoms_tot += in[b]._linmoms;

  vvd_t coord(vd_t(0.0, _linmoms_tot), npar);

  vvd_t y_Zq(vd_t(0.0, _linmoms_tot), njacks);                 // [njacks][moms]
  vd_t dy_Zq(0.0, _linmoms_tot);                               // [moms]
  vvvd_t y_Zbil(vvd_t(vd_t(0.0, _linmoms_tot), njacks), nbil); // [nbil][njacks][moms]
  vvd_t dy_Zbil(vd_t(0.0, _linmoms_tot), nbil);                // [nbil][moms]

  int j_tot = 0;

  for (int b = 0; b < nbeta; b++)
  {
    double ainv2 = ainv[b] * ainv[b];
    vvd_t dy_Zq_tmp = get<1>(ave_err_Zq(in[b].jZq));
    vvvvd_t dy_Zbil_tmp = get<1>(ave_err_Z(in[b].jZ)); // [moms][nbil][nmr][nmr]

    for (int j = 0; j < in[b]._linmoms; j++)
      if (in[b].p2[j] * ainv2 >= _p2min[b] && in[b].p2[j] * ainv2 <= _p2max[b])
      {
        // coord[0][j_tot] = (b==0)? 1.0 : 0.0; //Za
        // coord[1][j_tot] = (b==1)? 1.0 : 0.0; //Zb
        // coord[2][j_tot] = (b==2)? 1.0 : 0.0; //Zc
        for (int bb = 0; bb < nbeta; bb++)
          coord[bb][j_tot] = (bb == b) ? 1.0 : 0.0;            // Constant: Za,Zb,Zc
        coord[nbeta][j_tot] = in[b].p2[j];                     // a2p2 (lattice units)
        coord[nbeta + 1][j_tot] = 1.0 / (in[b].p2[j] * ainv2); // 1/GeV^2

        //   cout<<j_tot<<" "<<coord[0][j_tot]<<" "<<coord[1][j_tot]<<" "<<coord[2][j_tot]<<"
        //   "<<coord[3][j_tot]<<" "<<ainv2*coord[3][j_tot]<<" "<<coord[4][j_tot]<<endl;

        for (int ijack = 0; ijack < njacks; ijack++)
        {
          y_Zq[ijack][j_tot] = (in[b].jZq)[j][ijack][0];

          for (int ibil = 0; ibil < nbil; ibil++)
          {
            y_Zbil[ibil][ijack][j_tot] = (in[b].jZ)[j][ibil][ijack][0][0];
          }
        }

        dy_Zq[j_tot] = dy_Zq_tmp[j][0];
        for (int ibil = 0; ibil < nbil; ibil++)
          dy_Zbil[ibil][j_tot] = dy_Zbil_tmp[j][ibil][0][0];

        j_tot++;
      }
  }

  // for(int i=0; i<dy_Zq.size(); i++)
  //     cout<<"*** "<<y_Zq[0][i]<<"-"<<y_Zq[10][i]<<" "<<dy_Zq[i]<<endl;

  vvd_t jZq_pars = polyfit(coord, npar, dy_Zq, y_Zq, 0, _linmoms_tot - 1); // [ijack][ipar]

  vd_t jpole(0.0, njacks), jlincoeff(0.0, njacks), jchisq(0.0, njacks);

  for (int b = 0; b < nbeta; b++)
    cout << "p2 range (physical units):   " << _p2min[b] << " - " << _p2max[b] << endl;
  for (int ijack = 0; ijack < njacks; ijack++)
  {
    for (int b = 0; b < nbeta; b++)
      (out[b].jZq)[0][ijack][0] = jZq_pars[ijack][b];
    jlincoeff[ijack] = jZq_pars[ijack][nbeta];
    jpole[ijack] = jZq_pars[ijack][nbeta + 1];
    /**/
    jchisq[ijack] = jZq_pars[ijack][npar];
  }

  cout << "  -- pole[q] = " << get<0>(ave_err(jpole)) << "+/-" << get<1>(ave_err(jpole)) << endl;
  cout << "  -- lincoeff[q] = " << get<0>(ave_err(jlincoeff)) << "+/-" << get<1>(ave_err(jlincoeff))
       << endl;
  cout << "    ** chisqr[q] = " << get<0>(ave_err(jchisq)) << "+/-" << get<1>(ave_err(jchisq))
       << endl
       << endl;

  vector<string> str_bil = {"S", "V", "P", "A", "T"};
  for (int ibil = 0; ibil < nbil; ibil++)
  {
    vvd_t jZ_pars =
        polyfit(coord, npar, dy_Zbil[ibil], y_Zbil[ibil], 0, _linmoms_tot - 1); // [ijack][ipar]

    // vd_t jpole(0.0,njacks), jlincoeff(0.0,njacks), jchisq(0.0,njacks);

    for (int ijack = 0; ijack < njacks; ijack++)
    {
      for (int b = 0; b < nbeta; b++)
        (out[b].jZ)[0][ibil][ijack][0][0] = jZ_pars[ijack][b];
      jlincoeff[ijack] = jZ_pars[ijack][nbeta];
      jpole[ijack] = jZ_pars[ijack][nbeta + 1];
      /**/
      jchisq[ijack] = jZ_pars[ijack][npar];
    }

    cout << "  -- pole[" << str_bil[ibil] << "] = " << get<0>(ave_err(jpole)) << "+/-"
         << get<1>(ave_err(jpole)) << endl;
    cout << "  -- lincoeff[" << str_bil[ibil] << "] = " << get<0>(ave_err(jlincoeff)) << "+/-"
         << get<1>(ave_err(jlincoeff)) << endl;
    cout << "    ** chisqr[" << str_bil[ibil] << "] = " << get<0>(ave_err(jchisq)) << "+/-"
         << get<1>(ave_err(jchisq)) << endl
         << endl;
  }

  return out;
}

voper_t combined_M4(voper_t in) // M3 method combined on all betas
{

  cout << endl;
  cout << "----- Combined M4 fit on all the curves -----" << endl << endl;

  voper_t out(nbeta);

  for (int b = 0; b < nbeta; b++)
  {
    out[b] = in[b];

    out[b].linmoms = vector<array<int, 1>>{{0}};
    out[b].bilmoms = vector<array<int, 3>>{{0, 0, 0}};

    out[b]._linmoms = 1;
    out[b]._bilmoms = 1;

    out[b].allocate_val();
    out[b].allocate();
  }

  double _p2min[] = {p2min_M3_M4[0], p2min_M3_M4[1], p2min_M3_M4[2]}; // GeV^2
  double _p2max[] = {p2max_M3_M4[0], p2max_M3_M4[1], p2max_M3_M4[2]};

  // count all the momenta respecting the above criteria of p2min&p2max
  int _linmoms_tot = 0;
  for (int b = 0; b < nbeta; b++)
  {
    double ainv2 = ainv[b] * ainv[b];
    for (int j = 0; j < in[b]._linmoms; j++)
      if (in[b].p2[j] * ainv2 >= _p2min[b] && in[b].p2[j] * ainv2 <= _p2max[b])
        _linmoms_tot++;
  }

  int _linmoms_tot_complete = 0;
  for (int b = 0; b < nbeta; b++)
  {
    _linmoms_tot_complete += in[b]._linmoms;
  }
  cout << "tot moms = " << _linmoms_tot << "/" << _linmoms_tot_complete << endl;

  int npar = nbeta + 3; // {ZA,ZB,ZC}+{linear}+{quadratic}+{pole}

  vvd_t coord(vd_t(0.0, _linmoms_tot), npar);

  vvd_t y_Zq(vd_t(0.0, _linmoms_tot), njacks);                 // [njacks][moms]
  vd_t dy_Zq(0.0, _linmoms_tot);                               // [moms]
  vvvd_t y_Zbil(vvd_t(vd_t(0.0, _linmoms_tot), njacks), nbil); // [nbil][njacks][moms]
  vvd_t dy_Zbil(vd_t(0.0, _linmoms_tot), nbil);                // [nbil][moms]

  int j_tot = 0;

  for (int b = 0; b < nbeta; b++)
  {
    double ainv2 = ainv[b] * ainv[b];
    vvd_t dy_Zq_tmp = get<1>(ave_err_Zq(in[b].jZq));
    vvvvd_t dy_Zbil_tmp = get<1>(ave_err_Z(in[b].jZ)); // [moms][nbil][nmr][nmr]

    for (int j = 0; j < in[b]._linmoms; j++)
      if (in[b].p2[j] * ainv2 >= _p2min[b] && in[b].p2[j] * ainv2 <= _p2max[b])
      {
        // coord[0][j_tot] = (b==0)? 1.0 : 0.0; //Za
        // coord[1][j_tot] = (b==1)? 1.0 : 0.0; //Zb
        // coord[2][j_tot] = (b==2)? 1.0 : 0.0; //Zc
        for (int bb = 0; bb < nbeta; bb++)
          coord[bb][j_tot] = (bb == b) ? 1.0 : 0.0;            // Constant: Za,Zb,Zc
        coord[nbeta][j_tot] = in[b].p2[j];                     // a2p2 (lattice units)
        coord[nbeta + 1][j_tot] = 1.0 / (in[b].p2[j] * ainv2); // 1/GeV^2
        coord[nbeta + 2][j_tot] =
            coord[nbeta][j_tot] * coord[nbeta][j_tot]; // (a2p2)^2 (lattice units)

        for (int ijack = 0; ijack < njacks; ijack++)
        {
          y_Zq[ijack][j_tot] = (in[b].jZq)[j][ijack][0];

          for (int ibil = 0; ibil < nbil; ibil++)
          {
            y_Zbil[ibil][ijack][j_tot] = (in[b].jZ)[j][ibil][ijack][0][0];
          }
        }

        dy_Zq[j_tot] = dy_Zq_tmp[j][0];
        for (int ibil = 0; ibil < nbil; ibil++)
          dy_Zbil[ibil][j_tot] = dy_Zbil_tmp[j][ibil][0][0];

        j_tot++;
      }
  }

  vvd_t jZq_pars = polyfit(coord, npar, dy_Zq, y_Zq, 0, _linmoms_tot - 1); // [ijack][ipar]

  vd_t jpole(0.0, njacks), jlincoeff(0.0, njacks), jsqrcoeff(0.0, njacks), jchisq(0.0, njacks);

  for (int b = 0; b < nbeta; b++)
    cout << "p2 range (physical units):   " << _p2min[b] << " - " << _p2max[b] << endl;
  for (int ijack = 0; ijack < njacks; ijack++)
  {
    for (int b = 0; b < nbeta; b++)
      (out[b].jZq)[0][ijack][0] = jZq_pars[ijack][b];
    jlincoeff[ijack] = jZq_pars[ijack][nbeta];
    jpole[ijack] = jZq_pars[ijack][nbeta + 1];
    jsqrcoeff[ijack] = jZq_pars[ijack][nbeta + 2];
    /**/
    jchisq[ijack] = jZq_pars[ijack][npar];
  }

  cout << "  -- pole[q] = " << get<0>(ave_err(jpole)) << "+/-" << get<1>(ave_err(jpole)) << endl;
  cout << "  -- lincoeff[q] = " << get<0>(ave_err(jlincoeff)) << "+/-" << get<1>(ave_err(jlincoeff))
       << endl;
  cout << "  -- sqrcoeff[q] = " << get<0>(ave_err(jsqrcoeff)) << "+/-" << get<1>(ave_err(jsqrcoeff))
       << endl;
  cout << "    ** chisqr[q] = " << get<0>(ave_err(jchisq)) << "+/-" << get<1>(ave_err(jchisq))
       << endl
       << endl;

  vector<string> str_bil = {"S", "V", "P", "A", "T"};
  for (int ibil = 0; ibil < nbil; ibil++)
  {
    vvd_t jZ_pars =
        polyfit(coord, npar, dy_Zbil[ibil], y_Zbil[ibil], 0, _linmoms_tot - 1); // [ijack][ipar]

    // vd_t jpole(0.0,njacks), jlincoeff(0.0,njacks), jchisq(0.0,njacks);

    for (int ijack = 0; ijack < njacks; ijack++)
    {
      for (int b = 0; b < nbeta; b++)
        (out[b].jZ)[0][ibil][ijack][0][0] = jZ_pars[ijack][b];
      jlincoeff[ijack] = jZ_pars[ijack][nbeta];
      jpole[ijack] = jZ_pars[ijack][nbeta + 1];
      jsqrcoeff[ijack] = jZ_pars[ijack][nbeta + 2];
      /**/
      jchisq[ijack] = jZ_pars[ijack][npar];
    }

    cout << "  -- pole[" << str_bil[ibil] << "] = " << get<0>(ave_err(jpole)) << "+/-"
         << get<1>(ave_err(jpole)) << endl;
    cout << "  -- lincoeff[" << str_bil[ibil] << "] = " << get<0>(ave_err(jlincoeff)) << "+/-"
         << get<1>(ave_err(jlincoeff)) << endl;
    cout << "  -- sqrcoeff[" << str_bil[ibil] << "] = " << get<0>(ave_err(jsqrcoeff)) << "+/-"
         << get<1>(ave_err(jsqrcoeff)) << endl;
    cout << "    ** chisqr[" << str_bil[ibil] << "] = " << get<0>(ave_err(jchisq)) << "+/-"
         << get<1>(ave_err(jchisq)) << endl
         << endl;
  }

  return out;
}

voper_t combined_M5(voper_t in) // M3 method combined on all betas
{

  cout << endl;
  cout << "----- Combined M5 fit on all the curves -----" << endl << endl;

  voper_t out(nbeta);

  for (int b = 0; b < nbeta; b++)
  {
    out[b] = in[b];

    out[b].linmoms = vector<array<int, 1>>{{0}};
    out[b].bilmoms = vector<array<int, 3>>{{0, 0, 0}};

    out[b]._linmoms = 1;
    out[b]._bilmoms = 1;

    out[b].allocate_val();
    out[b].allocate();
  }

  double _p2min[] = {p2min_M3_M4[0], p2min_M3_M4[1], p2min_M3_M4[2]}; // GeV^2
  double _p2max[] = {p2max_M3_M4[0], p2max_M3_M4[1], p2max_M3_M4[2]};

  // count all the momenta respecting the above criteria of p2min&p2max
  int _linmoms_tot = 0;
  for (int b = 0; b < nbeta; b++)
  {
    double ainv2 = ainv[b] * ainv[b];
    for (int j = 0; j < in[b]._linmoms; j++)
      if (in[b].p2[j] * ainv2 >= _p2min[b] && in[b].p2[j] * ainv2 <= _p2max[b])
        _linmoms_tot++;
  }

  int _linmoms_tot_complete = 0;
  for (int b = 0; b < nbeta; b++)
  {
    _linmoms_tot_complete += in[b]._linmoms;
  }
  cout << "tot moms = " << _linmoms_tot << "/" << _linmoms_tot_complete << endl;

  int npar = 2 * nbeta + 1; // {ZA,ZB,ZC}+{linearA,linearB,linearC}+{pole}

  vvd_t coord(vd_t(0.0, _linmoms_tot), npar);

  vvd_t y_Zq(vd_t(0.0, _linmoms_tot), njacks);                 // [njacks][moms]
  vd_t dy_Zq(0.0, _linmoms_tot);                               // [moms]
  vvvd_t y_Zbil(vvd_t(vd_t(0.0, _linmoms_tot), njacks), nbil); // [nbil][njacks][moms]
  vvd_t dy_Zbil(vd_t(0.0, _linmoms_tot), nbil);                // [nbil][moms]

  int j_tot = 0;

  for (int b = 0; b < nbeta; b++)
  {
    double ainv2 = ainv[b] * ainv[b];
    vvd_t dy_Zq_tmp = get<1>(ave_err_Zq(in[b].jZq));
    vvvvd_t dy_Zbil_tmp = get<1>(ave_err_Z(in[b].jZ)); // [moms][nbil][nmr][nmr]

    for (int j = 0; j < in[b]._linmoms; j++)
      if (in[b].p2[j] * ainv2 >= _p2min[b] && in[b].p2[j] * ainv2 <= _p2max[b])
      {
        // coord[0][j_tot] = (b==0)? 1.0 : 0.0; //Za
        // coord[1][j_tot] = (b==1)? 1.0 : 0.0; //Zb
        // coord[2][j_tot] = (b==2)? 1.0 : 0.0; //Zc
        for (int bb = 0; bb < nbeta; bb++)
        {
          coord[bb][j_tot] = (bb == b) ? 1.0 : 0.0;                    // Constant: Za,Zb,Zc
          coord[bb + nbeta][j_tot] = (bb == b) ? (in[bb].p2[j]) : 0.0; // a2p2 (lattice units)
        }
        coord[2 * nbeta][j_tot] = 1.0 / (in[b].p2[j] * ainv2); // GeV^-2

        for (int ijack = 0; ijack < njacks; ijack++)
        {
          y_Zq[ijack][j_tot] = (in[b].jZq)[j][ijack][0];

          for (int ibil = 0; ibil < nbil; ibil++)
          {
            y_Zbil[ibil][ijack][j_tot] = (in[b].jZ)[j][ibil][ijack][0][0];
          }
        }

        dy_Zq[j_tot] = dy_Zq_tmp[j][0];
        for (int ibil = 0; ibil < nbil; ibil++)
          dy_Zbil[ibil][j_tot] = dy_Zbil_tmp[j][ibil][0][0];

        j_tot++;
      }
  }

  vvd_t jZq_pars = polyfit(coord, npar, dy_Zq, y_Zq, 0, _linmoms_tot - 1); // [ijack][ipar]

  vd_t jlinA(0.0, njacks), jlinB(0.0, njacks), jlinC(0.0, njacks), jpole(0.0, njacks),
      jchisq(0.0, njacks);

  for (int b = 0; b < nbeta; b++)
    cout << "p2 range (physical units):   " << _p2min[b] << " - " << _p2max[b] << endl;
  for (int ijack = 0; ijack < njacks; ijack++)
  {
    for (int b = 0; b < nbeta; b++)
      (out[b].jZq)[0][ijack][0] = jZq_pars[ijack][b];
    jlinA[ijack] = jZq_pars[ijack][nbeta];
    jlinB[ijack] = jZq_pars[ijack][nbeta + 1];
    jlinC[ijack] = jZq_pars[ijack][nbeta + 2];
    jpole[ijack] = jZq_pars[ijack][2 * nbeta];
    /**/
    jchisq[ijack] = jZq_pars[ijack][npar];
  }

  cout << "  -- lincoeffA[q] = " << get<0>(ave_err(jlinA)) << "+/-" << get<1>(ave_err(jlinA))
       << endl;
  cout << "  -- lincoeffB[q] = " << get<0>(ave_err(jlinB)) << "+/-" << get<1>(ave_err(jlinB))
       << endl;
  cout << "  -- lincoeffC[q] = " << get<0>(ave_err(jlinC)) << "+/-" << get<1>(ave_err(jlinC))
       << endl;
  cout << "  -- pole[q] = " << get<0>(ave_err(jpole)) << "+/-" << get<1>(ave_err(jpole)) << endl;
  cout << "    ** chisqr[q] = " << get<0>(ave_err(jchisq)) << "+/-" << get<1>(ave_err(jchisq))
       << endl
       << endl;

  vector<string> str_bil = {"S", "V", "P", "A", "T"};
  for (int ibil = 0; ibil < nbil; ibil++)
  {
    vvd_t jZ_pars =
        polyfit(coord, npar, dy_Zbil[ibil], y_Zbil[ibil], 0, _linmoms_tot - 1); // [ijack][ipar]

    // vd_t jpole(0.0,njacks), jlincoeff(0.0,njacks), jchisq(0.0,njacks);

    for (int ijack = 0; ijack < njacks; ijack++)
    {
      for (int b = 0; b < nbeta; b++)
        (out[b].jZ)[0][ibil][ijack][0][0] = jZ_pars[ijack][b];
      jlinA[ijack] = jZ_pars[ijack][nbeta];
      jlinB[ijack] = jZ_pars[ijack][nbeta + 1];
      jlinC[ijack] = jZ_pars[ijack][nbeta + 2];
      jpole[ijack] = jZ_pars[ijack][2 * nbeta];
      /**/
      jchisq[ijack] = jZ_pars[ijack][npar];
    }

    cout << "  -- lincoeffA[" << str_bil[ibil] << "] = " << get<0>(ave_err(jlinA)) << "+/-"
         << get<1>(ave_err(jlinA)) << endl;
    cout << "  -- lincoeffB[" << str_bil[ibil] << "] = " << get<0>(ave_err(jlinB)) << "+/-"
         << get<1>(ave_err(jlinB)) << endl;
    cout << "  -- lincoeffC[" << str_bil[ibil] << "] = " << get<0>(ave_err(jlinC)) << "+/-"
         << get<1>(ave_err(jlinC)) << endl;
    cout << "  -- pole[" << str_bil[ibil] << "] = " << get<0>(ave_err(jpole)) << "+/-"
         << get<1>(ave_err(jpole)) << endl;
    cout << "    ** chisqr[" << str_bil[ibil] << "] = " << get<0>(ave_err(jchisq)) << "+/-"
         << get<1>(ave_err(jchisq)) << endl
         << endl;
  }

  return out;
}

voper_t combined_M5_quadratic(voper_t in)
{

  cout << endl;
  cout << "----- Combined M5 + quadratic fit on all the curves -----" << endl << endl;

  voper_t out(nbeta);

  for (int b = 0; b < nbeta; b++)
  {
    out[b] = in[b];

    out[b].linmoms = vector<array<int, 1>>{{0}};
    out[b].bilmoms = vector<array<int, 3>>{{0, 0, 0}};

    out[b]._linmoms = 1;
    out[b]._bilmoms = 1;

    out[b].allocate_val();
    out[b].allocate();
  }

  double _p2min[] = {p2min_M3_M4[0], p2min_M3_M4[1], p2min_M3_M4[2]}; // GeV^2
  double _p2max[] = {p2max_M3_M4[0], p2max_M3_M4[1], p2max_M3_M4[2]};

  // count all the momenta respecting the above criteria of p2min&p2max
  int _linmoms_tot = 0;
  for (int b = 0; b < nbeta; b++)
  {
    double ainv2 = ainv[b] * ainv[b];
    for (int j = 0; j < in[b]._linmoms; j++)
      if (in[b].p2[j] * ainv2 >= _p2min[b] && in[b].p2[j] * ainv2 <= _p2max[b])
        _linmoms_tot++;
  }

  int _linmoms_tot_complete = 0;
  for (int b = 0; b < nbeta; b++)
  {
    _linmoms_tot_complete += in[b]._linmoms;
  }
  cout << "tot moms = " << _linmoms_tot << "/" << _linmoms_tot_complete << endl;

  int npar = nbeta + 3; // {ZA,ZB,ZC}+{linear}+{quad}+{pole}

  vvd_t coord(vd_t(0.0, _linmoms_tot), npar);

  vvd_t y_Zq(vd_t(0.0, _linmoms_tot), njacks);                 // [njacks][moms]
  vd_t dy_Zq(0.0, _linmoms_tot);                               // [moms]
  vvvd_t y_Zbil(vvd_t(vd_t(0.0, _linmoms_tot), njacks), nbil); // [nbil][njacks][moms]
  vvd_t dy_Zbil(vd_t(0.0, _linmoms_tot), nbil);                // [nbil][moms]

  int j_tot = 0;

  for (int b = 0; b < nbeta; b++)
  {
    double ainv2 = ainv[b] * ainv[b];
    vvd_t dy_Zq_tmp = get<1>(ave_err_Zq(in[b].jZq));
    vvvvd_t dy_Zbil_tmp = get<1>(ave_err_Z(in[b].jZ)); // [moms][nbil][nmr][nmr]

    for (int j = 0; j < in[b]._linmoms; j++)
      if (in[b].p2[j] * ainv2 >= _p2min[b] && in[b].p2[j] * ainv2 <= _p2max[b])
      {
        // coord[0][j_tot] = (b==0)? 1.0 : 0.0; //Za
        // coord[1][j_tot] = (b==1)? 1.0 : 0.0; //Zb
        // coord[2][j_tot] = (b==2)? 1.0 : 0.0; //Zc
        for (int bb = 0; bb < nbeta; bb++)
        {
          coord[bb][j_tot] = (bb == b) ? 1.0 : 0.0; // Constant: Za,Zb,Zc
        }
        coord[nbeta][j_tot] = 1.0 / (in[b].p2[j] * ainv2);       // GeV^-2
        coord[nbeta + 1][j_tot] = in[b].p2[j];                   // a2p2 [lattice units]
        coord[nbeta + 2][j_tot] = (in[b].p2[j]) * (in[b].p2[j]); // (a2p2)^2 [lattice units]

        for (int ijack = 0; ijack < njacks; ijack++)
        {
          y_Zq[ijack][j_tot] = (in[b].jZq)[j][ijack][0];

          for (int ibil = 0; ibil < nbil; ibil++)
          {
            y_Zbil[ibil][ijack][j_tot] =
                (in[b].jZ)[j][ibil][ijack][0][0]; //-0.059/(in[b].p2[j]*ainv2);
          }
        }

        dy_Zq[j_tot] = dy_Zq_tmp[j][0];
        for (int ibil = 0; ibil < nbil; ibil++)
          dy_Zbil[ibil][j_tot] = dy_Zbil_tmp[j][ibil][0][0];

        j_tot++;
      }
  }

  vvd_t jZq_pars = polyfit(coord, npar, dy_Zq, y_Zq, 0, _linmoms_tot - 1); // [ijack][ipar]

  vd_t jpole(0.0, njacks), jd2(0.0, njacks), jd4(0.0, njacks), jchisq(0.0, njacks);

  for (int b = 0; b < nbeta; b++)
    cout << "p2 range (physical units):   " << _p2min[b] << " - " << _p2max[b] << endl;
  for (int ijack = 0; ijack < njacks; ijack++)
  {
    for (int b = 0; b < nbeta; b++)
      (out[b].jZq)[0][ijack][0] = jZq_pars[ijack][b];
    jpole[ijack] = jZq_pars[ijack][nbeta];
    jd2[ijack] = jZq_pars[ijack][nbeta + 1];
    jd4[ijack] = jZq_pars[ijack][nbeta + 2];
    /**/
    jchisq[ijack] = jZq_pars[ijack][npar];
  }

  cout << "  -- pole[q] = " << get<0>(ave_err(jpole)) << "+/-" << get<1>(ave_err(jpole)) << endl;
  cout << "  -- d2[q] = " << get<0>(ave_err(jd2)) << "+/-" << get<1>(ave_err(jd2)) << endl;
  cout << "  -- d4[q] = " << get<0>(ave_err(jd4)) << "+/-" << get<1>(ave_err(jd4)) << endl;
  cout << "    ** chisqr[q] = " << get<0>(ave_err(jchisq)) << "+/-" << get<1>(ave_err(jchisq))
       << endl
       << endl;

  vector<string> str_bil = {"S", "V", "P", "A", "T"};
  for (int ibil = 0; ibil < nbil; ibil++)
  {
    vvd_t jZ_pars =
        polyfit(coord, npar, dy_Zbil[ibil], y_Zbil[ibil], 0, _linmoms_tot - 1); // [ijack][ipar]

    for (int ijack = 0; ijack < njacks; ijack++)
    {
      for (int b = 0; b < nbeta; b++)
        (out[b].jZ)[0][ibil][ijack][0][0] = jZ_pars[ijack][b];
      jpole[ijack] = jZ_pars[ijack][nbeta];
      jd2[ijack] = jZ_pars[ijack][nbeta + 1];
      jd4[ijack] = jZ_pars[ijack][nbeta + 2];
      /**/
      jchisq[ijack] = jZ_pars[ijack][npar];
    }

    cout << "  -- pole[" << str_bil[ibil] << "] = " << get<0>(ave_err(jpole)) << "+/-"
         << get<1>(ave_err(jpole)) << endl;
    cout << "  -- d2[" << str_bil[ibil] << "] = " << get<0>(ave_err(jd2)) << "+/-"
         << get<1>(ave_err(jd2)) << endl;
    cout << "  -- d4[" << str_bil[ibil] << "] = " << get<0>(ave_err(jd4)) << "+/-"
         << get<1>(ave_err(jd4)) << endl;
    cout << "    ** chisqr[" << str_bil[ibil] << "] = " << get<0>(ave_err(jchisq)) << "+/-"
         << get<1>(ave_err(jchisq)) << endl
         << endl;
  }

  return out;
}

voper_t combined_M5_log(voper_t in) // M3 method combined on all betas
{

  cout << endl;
  cout << "----- Combined M5 fit on all the curves -----" << endl << endl;

  voper_t out(nbeta);

  for (int b = 0; b < nbeta; b++)
  {
    out[b] = in[b];

    out[b].linmoms = vector<array<int, 1>>{{0}};
    out[b].bilmoms = vector<array<int, 3>>{{0, 0, 0}};

    out[b]._linmoms = 1;
    out[b]._bilmoms = 1;

    out[b].allocate_val();
    out[b].allocate();
  }

  double _p2min[] = {p2min_M3_M4[0], p2min_M3_M4[1], p2min_M3_M4[2]}; // GeV^2
  double _p2max[] = {p2max_M3_M4[0], p2max_M3_M4[1], p2max_M3_M4[2]};

  // count all the momenta respecting the above criteria of p2min&p2max
  int _linmoms_tot = 0;
  for (int b = 0; b < nbeta; b++)
  {
    double ainv2 = ainv[b] * ainv[b];
    for (int j = 0; j < in[b]._linmoms; j++)
      if (in[b].p2[j] * ainv2 >= _p2min[b] && in[b].p2[j] * ainv2 <= _p2max[b])
        _linmoms_tot++;
  }

  int _linmoms_tot_complete = 0;
  for (int b = 0; b < nbeta; b++)
  {
    _linmoms_tot_complete += in[b]._linmoms;
  }
  cout << "tot moms = " << _linmoms_tot << "/" << _linmoms_tot_complete << endl;

  int npar = nbeta + 2; // {ZA,ZB,ZC}+{linear}+{linear/log}+{linear/log^2}//+{pole}

  vvd_t coord(vd_t(0.0, _linmoms_tot), npar);

  vvd_t y_Zq(vd_t(0.0, _linmoms_tot), njacks);                 // [njacks][moms]
  vd_t dy_Zq(0.0, _linmoms_tot);                               // [moms]
  vvvd_t y_Zbil(vvd_t(vd_t(0.0, _linmoms_tot), njacks), nbil); // [nbil][njacks][moms]
  vvd_t dy_Zbil(vd_t(0.0, _linmoms_tot), nbil);                // [nbil][moms]

  int j_tot = 0;

  for (int b = 0; b < nbeta; b++)
  {
    double ainv2 = ainv[b] * ainv[b];
    vvd_t dy_Zq_tmp = get<1>(ave_err_Zq(in[b].jZq));
    vvvvd_t dy_Zbil_tmp = get<1>(ave_err_Z(in[b].jZ)); // [moms][nbil][nmr][nmr]

    for (int j = 0; j < in[b]._linmoms; j++)
      if (in[b].p2[j] * ainv2 >= _p2min[b] && in[b].p2[j] * ainv2 <= _p2max[b])
      {
        // coord[0][j_tot] = (b==0)? 1.0 : 0.0; //Za
        // coord[1][j_tot] = (b==1)? 1.0 : 0.0; //Zb
        // coord[2][j_tot] = (b==2)? 1.0 : 0.0; //Zc
        for (int bb = 0; bb < nbeta; bb++)
        {
          coord[bb][j_tot] = (bb == b) ? 1.0 : 0.0; // Constant: Za,Zb,Zc
        }
        //   coord[nbeta][j_tot] = 1.0/(in[b].p2[j]*ainv2);       // GeV^-2
        coord[nbeta][j_tot] = in[b].p2[j];                        // a2p2 [lattice units]
        coord[nbeta + 1][j_tot] = in[b].p2[j] * log(in[b].p2[j]); // a2p2*log(a2p2) [lattice units]
        //   coord[nbeta+2][j_tot] = in[b].p2[j]*log(in[b].p2[j])*log(in[b].p2[j]); //
        //   a2p2*log^2(a2p2) [lattice units]

        for (int ijack = 0; ijack < njacks; ijack++)
        {
          y_Zq[ijack][j_tot] = (in[b].jZq)[j][ijack][0];

          for (int ibil = 0; ibil < nbil; ibil++)
          {
            y_Zbil[ibil][ijack][j_tot] =
                (in[b].jZ)[j][ibil][ijack][0][0]; //-0.059/(in[b].p2[j]*ainv2);
          }
        }

        dy_Zq[j_tot] = dy_Zq_tmp[j][0];
        for (int ibil = 0; ibil < nbil; ibil++)
          dy_Zbil[ibil][j_tot] = dy_Zbil_tmp[j][ibil][0][0];

        j_tot++;
      }
  }

  vvd_t jZq_pars = polyfit(coord, npar, dy_Zq, y_Zq, 0, _linmoms_tot - 1); // [ijack][ipar]

  vd_t jd2(0.0, njacks), jd2L(0.0, njacks), jd2LL(0.0, njacks) /*, jpole(0.0,njacks)*/,
      jchisq(0.0, njacks);

  for (int b = 0; b < nbeta; b++)
    cout << "p2 range (physical units):   " << _p2min[b] << " - " << _p2max[b] << endl;
  for (int ijack = 0; ijack < njacks; ijack++)
  {
    for (int b = 0; b < nbeta; b++)
      (out[b].jZq)[0][ijack][0] = jZq_pars[ijack][b];
    //   jpole[ijack] = jZq_pars[ijack][nbeta];
    jd2[ijack] = jZq_pars[ijack][nbeta];
    jd2L[ijack] = jZq_pars[ijack][nbeta + 1];
    //   jd2LL[ijack] = jZq_pars[ijack][nbeta+2];
    /**/
    jchisq[ijack] = jZq_pars[ijack][npar];
  }

  // cout<<"  -- pole[q] = "<<get<0>(ave_err(jpole))<<"+/-"<<get<1>(ave_err(jpole))<<endl;
  cout << "  -- d2[q] = " << get<0>(ave_err(jd2)) << "+/-" << get<1>(ave_err(jd2)) << endl;
  cout << "  -- d2L[q] = " << get<0>(ave_err(jd2L)) << "+/-" << get<1>(ave_err(jd2L)) << endl;
  // cout<<"  -- d2LL[q] = "<<get<0>(ave_err(jd2LL))<<"+/-"<<get<1>(ave_err(jd2LL))<<endl;
  cout << "    ** chisqr[q] = " << get<0>(ave_err(jchisq)) << "+/-" << get<1>(ave_err(jchisq))
       << endl
       << endl;

  vector<string> str_bil = {"S", "V", "P", "A", "T"};
  for (int ibil = 0; ibil < nbil; ibil++)
  {
    vvd_t jZ_pars =
        polyfit(coord, npar, dy_Zbil[ibil], y_Zbil[ibil], 0, _linmoms_tot - 1); // [ijack][ipar]

    // vd_t jpole(0.0,njacks), jlincoeff(0.0,njacks), jchisq(0.0,njacks);

    for (int ijack = 0; ijack < njacks; ijack++)
    {
      for (int b = 0; b < nbeta; b++)
        (out[b].jZ)[0][ibil][ijack][0][0] = jZ_pars[ijack][b];
      // jpole[ijack] = jZ_pars[ijack][nbeta];
      jd2[ijack] = jZ_pars[ijack][nbeta];
      jd2L[ijack] = jZ_pars[ijack][nbeta + 1];
      // jd2LL[ijack] = jZ_pars[ijack][nbeta+2];
      /**/
      jchisq[ijack] = jZ_pars[ijack][npar];
    }

    // cout<<"  -- pole["<<str_bil[ibil]<<"] =
    // "<<get<0>(ave_err(jpole))<<"+/-"<<get<1>(ave_err(jpole))<<endl;
    cout << "  -- d2[" << str_bil[ibil] << "] = " << get<0>(ave_err(jd2)) << "+/-"
         << get<1>(ave_err(jd2)) << endl;
    cout << "  -- d2L[" << str_bil[ibil] << "] = " << get<0>(ave_err(jd2L)) << "+/-"
         << get<1>(ave_err(jd2L)) << endl;
    // cout<<"  -- d2LL["<<str_bil[ibil]<<"] =
    // "<<get<0>(ave_err(jd2LL))<<"+/-"<<get<1>(ave_err(jd2LL))<<endl;
    cout << "    ** chisqr[" << str_bil[ibil] << "] = " << get<0>(ave_err(jchisq)) << "+/-"
         << get<1>(ave_err(jchisq)) << endl
         << endl;
  }

  return out;
}

void oper_t::plot_bil_chir_extr(int mom, int i_ins, int ibil, vd_t x, vd_t y, vd_t dy, vvd_t jpars,
                                string suffix)
{
  // this choice is relative to the twisted basis
  vector<string> bil = {"S", "V", "P", "A", "T"};
  vector<string> ins_str = {"LO", "EM"}; /* valid only if ntypes=3 */

  if (ntypes != 3 and ntypes != 1)
  {
    cout << "plot of chiral extrapolation not implemented for ntypes=" << ntypes << "!" << endl;
    exit(0);
  }

  ofstream plot_data;
  ofstream pars_data;

  plot_data.open(path_to_ens + "plots/chir_extr_Gbil_" + bil[ibil] + "_" + ins_str[i_ins] +
                 "_mom_" + to_string(mom) + (suffix != "" ? ("_" + suffix) : string("")) + ".txt");
  pars_data.open(path_to_ens + "plots/chir_extr_Gbil_" + bil[ibil] + "_" + ins_str[i_ins] +
                 "_mom_" + to_string(mom) + "_pars" + (suffix != "" ? ("_" + suffix) : string("")) +
                 ".txt");

  int npar = (int)jpars[0].size();

  vd_t pars_ave(0.0, npar);
  vd_t sqr_pars_ave(0.0, npar);
  vd_t pars_err(0.0, npar);

  for (int ipar = 0; ipar < npar; ipar++)
  {
    for (int ijack = 0; ijack < njacks; ijack++)
    {
      pars_ave[ipar] += jpars[ijack][ipar] / njacks;
      sqr_pars_ave[ipar] += jpars[ijack][ipar] * jpars[ijack][ipar] / njacks;
    }
    pars_err[ipar] = sqrt((double)(njacks - 1)) *
                     sqrt(fabs(sqr_pars_ave[ipar] - pars_ave[ipar] * pars_ave[ipar]));

    if (ipar == npar - 1)
      pars_data << endl;
    pars_data << pars_ave[ipar] << "\t" << pars_err[ipar] << endl; /* includes also the chi2 */
  }

  plot_data << "0.0\t" << pars_ave[0] << "\t" << pars_err[0] << endl;
  for (int i = 0; i < (int)x.size(); i++)
    plot_data << x[i] << "\t" << y[i] << "\t" << dy[i] << endl;
}

void oper_t::plot(const string suffix)
{
  oper_t in = (*this);

  // Zq
  Zq_tup Zq_ave_err = ave_err_Zq(in.jZq);
  vvd_t Zq_ave = get<0>(Zq_ave_err); //[imom][mr]
  vvd_t Zq_err = get<1>(Zq_ave_err); //[imom][mr]

  // Gbil
  tuple<vvvvvd_t, vvvvvd_t> Gbil_ave_err = ave_err(in.jG);
  vvvvvd_t G_ave = get<0>(Gbil_ave_err); //[imom][ins][ibil][mrA][mrB]
  vvvvvd_t G_err = get<1>(Gbil_ave_err); //[imom][ins][ibil][mrA][mrB]

  // Zbil
  Zbil_tup Zbil_ave_err = ave_err_Z(in.jZ);
  vvvvd_t Z_ave = get<0>(Zbil_ave_err); //[imom][ibil][mr1][mr2]
  vvvvd_t Z_err = get<1>(Zbil_ave_err); //[imom][ibil][mr1][mr2]

  // ZV/ZA and ZP/ZS and ZA/ZV (S,V,P,A,T)
  Zbil_tup ZVovZA_ave_err = ave_err_Z(in.jZVoverZA);
  Zbil_tup ZPovZS_ave_err = ave_err_Z(in.jZPoverZS);
  Zbil_tup ZAovZV_ave_err = ave_err_Z(in.jZAoverZV);
  vvvvd_t ZVovZA_ave = get<0>(ZVovZA_ave_err); //[imom][0][mr1][mr2]
  vvvvd_t ZPovZS_ave = get<0>(ZPovZS_ave_err);
  vvvvd_t ZAovZV_ave = get<0>(ZAovZV_ave_err);
  vvvvd_t ZVovZA_err = get<1>(ZVovZA_ave_err);
  vvvvd_t ZPovZS_err = get<1>(ZPovZS_ave_err);
  vvvvd_t ZAovZV_err = get<1>(ZAovZV_ave_err);

  // bval/bsea
  Zq_tup bval_ave_err = ave_err_Zq(in.bval);
  Zq_tup bsea_ave_err = ave_err_Zq(in.bsea);
  vvd_t bval_ave = get<0>(bval_ave_err);
  vvd_t bval_err = get<1>(bval_ave_err);
  vvd_t bsea_ave = get<0>(bsea_ave_err);
  vvd_t bsea_err = get<1>(bsea_ave_err);

  // stepfunc
  tuple<vd_t, vd_t> stepfunc_ave_err = ave_err(in.stepfunc);
  vd_t stepfunc_ave = get<0>(stepfunc_ave_err);
  vd_t stepfunc_err = get<1>(stepfunc_ave_err);

  // this choice is relative to the twisted basis
  vector<string> bil = {"S", "V", "P", "A", "T"};

  ofstream Zq_data;
  ofstream Zq_p2_data;

  vector<ofstream> Gbil_data(nbil);
  vector<ofstream> Gbil_p2_data(nbil);

  vector<ofstream> Zbil_data(nbil);
  vector<ofstream> Zbil_p2_data(nbil);

  ofstream ZVovZA_data, ZPovZS_data, ZAovZV_data;
  ofstream ZVovZA_p2_data, ZPovZS_p2_data, ZAovZV_p2_data;

  vector<ofstream> bval_data(nbil + 1), bval_p2_data(nbil + 1); // bil + Zq
  vector<ofstream> bsea_data(nbil + 1), bsea_p2_data(nbil + 1);

  vector<ofstream> stepfunc_p2_data(nbil + 1 + 3); // bil + Zq + ZA/ZV,ZP/ZS,ZV/ZA

  int m_tmp = 0;

  cout << "Plotting Zq";

  if (_nr == 1) // only after r-average
  {
    for (int m = 0; m < in._nm; m++)
    {
      m_tmp = mass_val[m] * 10000;

      if (in._nm == 1)
      {
        Zq_data.open(path_to_ens + "plots/Zq" + (suffix != "" ? ("_" + suffix) : string("")) +
                     ".txt");
        Zq_p2_data.open(path_to_ens + "plots/Zq" + (suffix != "" ? ("_" + suffix) : string("")) +
                        "_p2.txt");
      }
      else
      {
        Zq_data.open(path_to_ens + "plots/Zq" + (suffix != "" ? ("_" + suffix) : string("")) +
                     "_m" + to_string(m_tmp) + ".txt");
        Zq_p2_data.open(path_to_ens + "plots/Zq" + (suffix != "" ? ("_" + suffix) : string("")) +
                        "_p2_m" + to_string(m_tmp) + ".txt");
      }

      for (int imom = 0; imom < in._linmoms; imom++)
      {
        Zq_data << (in.p2_tilde)[imom] << "\t" << Zq_ave[imom][m] << "\t" << Zq_err[imom][m]
                << endl;
        Zq_p2_data << (in.p2)[imom] << "\t" << Zq_ave[imom][m] << "\t" << Zq_err[imom][m] << endl;
      }

      Zq_data.close();
      Zq_p2_data.close();
    }
  }
  else
  {
    Zq_data.open(path_to_ens + "plots/Zq" + (suffix != "" ? ("_" + suffix) : string("")) + ".txt");
    Zq_p2_data.open(path_to_ens + "plots/Zq" + (suffix != "" ? ("_" + suffix) : string("")) +
                    "_p2.txt");

    for (int imom = 0; imom < in._linmoms; imom++)
    {
      Zq_data << (in.p2_tilde)[imom] << "\t" << Zq_ave[imom][0] << "\t" << Zq_err[imom][0] << endl;
      Zq_p2_data << (in.p2)[imom] << "\t" << Zq_ave[imom][0] << "\t" << Zq_err[imom][0] << endl;
    }

    Zq_data.close();
    Zq_p2_data.close();
  }

  cout << ", Zbil";

  if (in._nr == 1)
  {
    for (int m = 0; m < in._nm; m++)
    {
      m_tmp = mass_val[m] * 10000;

      for (int ibil = 0; ibil < nbil; ibil++)
      {
        if (suffix == "ave")
        {
          Gbil_data[ibil].open(path_to_ens + "plots/G" + bil[ibil] +
                               (suffix != "" ? ("_" + suffix) : string("")) + "_m" +
                               to_string(m_tmp) + ".txt");
          Gbil_p2_data[ibil].open(path_to_ens + "plots/G" + bil[ibil] +
                                  (suffix != "" ? ("_" + suffix) : string("")) + "_p2_m" +
                                  to_string(m_tmp) + ".txt");
        }

        Zbil_data[ibil].open(path_to_ens + "plots/Z" + bil[ibil] +
                             (suffix != "" ? ("_" + suffix) : string("")) + "_m" +
                             to_string(m_tmp) + ".txt");
        Zbil_p2_data[ibil].open(path_to_ens + "plots/Z" + bil[ibil] +
                                (suffix != "" ? ("_" + suffix) : string("")) + "_p2_m" +
                                to_string(m_tmp) + ".txt");

        for (int imom = 0; imom < in._bilmoms; imom++)
        {
          int imomk = imom; // NB: it works only for RIMOM!

          Gbil_data[ibil] << (in.p2_tilde)[imomk] << "\t" << G_ave[imom][0][ibil][m][m] << "\t"
                          << G_err[imom][0][ibil][m][m] << endl;
          Gbil_p2_data[ibil] << (in.p2)[imomk] << "\t" << G_ave[imom][0][ibil][m][m] << "\t"
                             << G_err[imom][0][ibil][m][m] << endl;

          Zbil_data[ibil] << (in.p2_tilde)[imomk] << "\t" << Z_ave[imom][ibil][m][m] << "\t"
                          << Z_err[imom][ibil][m][m] << endl;
          Zbil_p2_data[ibil] << (in.p2)[imomk] << "\t" << Z_ave[imom][ibil][m][m] << "\t"
                             << Z_err[imom][ibil][m][m] << endl;
        }

        if (suffix == "ave")
        {
          Gbil_data[ibil].close();
          Gbil_p2_data[ibil].close();
        }

        Zbil_data[ibil].close();
        Zbil_p2_data[ibil].close();
      }
    }
  }
  else
  {
    for (int ibil = 0; ibil < nbil; ibil++)
    {
      Zbil_data[ibil].open(path_to_ens + "plots/Z" + bil[ibil] +
                           (suffix != "" ? ("_" + suffix) : string("")) + ".txt");
      Zbil_p2_data[ibil].open(path_to_ens + "plots/Z" + bil[ibil] +
                              (suffix != "" ? ("_" + suffix) : string("")) + "_p2.txt");

      for (int imom = 0; imom < in._bilmoms; imom++)
      {
        int imomk = imom; // NB: it works only for RIMOM!

        Zbil_data[ibil] << (in.p2_tilde)[imomk] << "\t" << Z_ave[imom][ibil][0][0] << "\t"
                        << Z_err[imom][ibil][0][0] << endl;
        Zbil_p2_data[ibil] << (in.p2)[imomk] << "\t" << Z_ave[imom][ibil][0][0] << "\t"
                           << Z_err[imom][ibil][0][0] << endl;
      }

      Zbil_data[ibil].close();
      Zbil_p2_data[ibil].close();
    }
  }

  cout << ", ZV/ZA, ZP/ZS and ZA/ZV";
  ZVovZA_data.open(path_to_ens + "plots/ZVovZA" + (suffix != "" ? ("_" + suffix) : string("")) +
                   ".txt");
  ZPovZS_data.open(path_to_ens + "plots/ZPovZS" + (suffix != "" ? ("_" + suffix) : string("")) +
                   ".txt");
  ZAovZV_data.open(path_to_ens + "plots/ZAovZV" + (suffix != "" ? ("_" + suffix) : string("")) +
                   ".txt");
  ZVovZA_p2_data.open(path_to_ens + "plots/ZVovZA" + (suffix != "" ? ("_" + suffix) : string("")) +
                      "_p2.txt");
  ZPovZS_p2_data.open(path_to_ens + "plots/ZPovZS" + (suffix != "" ? ("_" + suffix) : string("")) +
                      "_p2.txt");
  ZAovZV_p2_data.open(path_to_ens + "plots/ZAovZV" + (suffix != "" ? ("_" + suffix) : string("")) +
                      "_p2.txt");

  for (int imom = 0; imom < in._bilmoms; imom++)
  {
    //            int imomq = in.bilmoms[imom][0];
    //            cout<<"imomq: "<<imomq<<endl;
    //            int imomk = in.linmoms[imomq][0];
    int imomk = imom; // NB: it works only for RIMOM!

    // ZVovZA_data<<(in.p2_tilde)[imomk]<<"\t"<<ZVovZA_ave[imom][0][0]<<"\t"<<ZVovZA_err[imom][0][0]<<endl;
    // ZPovZS_data<<(in.p2_tilde)[imomk]<<"\t"<<ZPovZS_ave[imom][0][0]<<"\t"<<ZPovZS_err[imom][0][0]<<endl;
    // ZAovZV_data<<(in.p2_tilde)[imomk]<<"\t"<<ZAovZV_ave[imom][0][0]<<"\t"<<ZAovZV_err[imom][0][0]<<endl;

    // ZVovZA_p2_data<<(in.p2)[imomk]<<"\t"<<ZVovZA_ave[imom][0][0]<<"\t"<<ZVovZA_err[imom][0][0]<<endl;
    // ZPovZS_p2_data<<(in.p2)[imomk]<<"\t"<<ZPovZS_ave[imom][0][0]<<"\t"<<ZPovZS_err[imom][0][0]<<endl;
    // ZAovZV_p2_data<<(in.p2)[imomk]<<"\t"<<ZAovZV_ave[imom][0][0]<<"\t"<<ZAovZV_err[imom][0][0]<<endl;

    ZVovZA_data << (in.p2_tilde)[imomk] << "\t" << ZVovZA_ave[imom][0][0][0] << "\t"
                << ZVovZA_err[imom][0][0][0] << endl;
    ZPovZS_data << (in.p2_tilde)[imomk] << "\t" << ZPovZS_ave[imom][0][0][0] << "\t"
                << ZPovZS_err[imom][0][0][0] << endl;
    ZAovZV_data << (in.p2_tilde)[imomk] << "\t" << ZAovZV_ave[imom][0][0][0] << "\t"
                << ZAovZV_err[imom][0][0][0] << endl;

    ZVovZA_p2_data << (in.p2)[imomk] << "\t" << ZVovZA_ave[imom][0][0][0] << "\t"
                   << ZVovZA_err[imom][0][0][0] << endl;
    ZPovZS_p2_data << (in.p2)[imomk] << "\t" << ZPovZS_ave[imom][0][0][0] << "\t"
                   << ZPovZS_err[imom][0][0][0] << endl;
    ZAovZV_p2_data << (in.p2)[imomk] << "\t" << ZAovZV_ave[imom][0][0][0] << "\t"
                   << ZAovZV_err[imom][0][0][0] << endl;
  }

  if (suffix == "chir" and strcmp(chir_ansatz_val.c_str(), "linear") == 0)
  {
    vector<string> bb = {"q", "S", "V", "P", "A", "T"};

    cout << ", bval ";
    for (int ib = 0; ib < nbil + 1; ib++)
    {
      bval_data[ib].open(path_to_ens + "plots/bval_" + bb[ib] + ".txt");
      bval_p2_data[ib].open(path_to_ens + "plots/bval_" + bb[ib] + "_p2.txt");

      for (int imom = 0; imom < in._linmoms; imom++)
      {
        bval_data[ib] << (in.p2_tilde)[imom] << "\t" << bval_ave[imom][ib] << "\t"
                      << bval_err[imom][ib] << endl;
        bval_p2_data[ib] << (in.p2)[imom] << "\t" << bval_ave[imom][ib] << "\t"
                         << bval_err[imom][ib] << endl;
      }
    }
  }
  if (suffix == "sea" and strcmp(chir_ansatz_sea.c_str(), "linear") == 0)
  {
    vector<string> bb = {"q", "S", "V", "P", "A", "T"};

    cout << ", bsea ";
    for (int ib = 0; ib < nbil + 1; ib++)
    {
      bsea_data[ib].open(path_to_ens + "plots/bsea_" + bb[ib] + ".txt");
      bsea_p2_data[ib].open(path_to_ens + "plots/bsea_" + bb[ib] + "_p2.txt");

      for (int imom = 0; imom < in._linmoms; imom++)
      {
        bsea_data[ib] << (in.p2_tilde)[imom] << "\t" << bsea_ave[imom][ib] << "\t"
                      << bsea_err[imom][ib] << endl;
        bsea_p2_data[ib] << (in.p2)[imom] << "\t" << bsea_ave[imom][ib] << "\t"
                         << bsea_err[imom][ib] << endl;
      }
    }
  }
  if (suffix == "evo" or suffix == "cont")
  {
    vector<string> bb = {"q", "S", "V", "P", "A", "T", "VA", "SP", "AV"};

    cout << ", stepfunc ";
    for (int ib = 0; ib < nbil + 1 + 2; ib++)
    {
      stepfunc_p2_data[ib].open(path_to_ens + "plots/stepfunc_" + bb[ib] + "_p2_" + suffix +
                                ".txt");
      stepfunc_p2_data[ib] << stepfunc_ave[ib] << "\t" << stepfunc_err[ib] << endl;
    }
  }

  cout << " in \"" << path_to_ens << "plots\"" << endl;
}