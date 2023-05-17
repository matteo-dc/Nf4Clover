#include "aliases.hpp"
#include "allocate.hpp"
#include "global.hpp"
#include <iostream>
#include <omp.h>
#include <tuple>

// average bilinears
tuple<vvvvvd_t, vvvvvd_t> ave_err(vector<jproj_t> jG)
{
  int _bilmoms = (int)jG.size();
  int _nins = (int)jG[0].size();
  int _nbil = (int)jG[0][0].size();
  int _njacks = (int)jG[0][0][0].size();
  int _nmr = (int)jG[0][0][0][0].size();

  vvvvvd_t G_ave(vvvvd_t(vvvd_t(vvd_t(vd_t(0.0, _nmr), _nmr), _nbil), _nins), _bilmoms);
  vvvvvd_t sqr_G_ave = G_ave;
  vvvvvd_t G_err = G_ave;

  for (int imom = 0; imom < _bilmoms; imom++)
  {
#pragma omp parallel for collapse(4)
    for (int ins = 0; ins < _nins; ins++)
      for (int ibil = 0; ibil < _nbil; ibil++)
        for (int mrA = 0; mrA < _nmr; mrA++)
          for (int mrB = 0; mrB < _nmr; mrB++)
            for (int ijack = 0; ijack < _njacks; ijack++)
            {
              G_ave[imom][ins][ibil][mrA][mrB] += jG[imom][ins][ibil][ijack][mrA][mrB] / _njacks;
              sqr_G_ave[imom][ins][ibil][mrA][mrB] +=
                  (jG[imom][ins][ibil][ijack][mrA][mrB] * jG[imom][ins][ibil][ijack][mrA][mrB]) /
                  _njacks;
            }
#pragma omp parallel for collapse(4)
    for (int ins = 0; ins < _nins; ins++)
      for (int ibil = 0; ibil < _nbil; ibil++)
        for (int mrA = 0; mrA < _nmr; mrA++)
          for (int mrB = 0; mrB < _nmr; mrB++)
          {
            G_err[imom][ins][ibil][mrA][mrB] =
                sqrt((double)(njacks - 1)) *
                sqrt(fabs(sqr_G_ave[imom][ins][ibil][mrA][mrB] -
                          G_ave[imom][ins][ibil][mrA][mrB] * G_ave[imom][ins][ibil][mrA][mrB]));
          }
  }

  tuple<vvvvvd_t, vvvvvd_t> tuple_ave_err(G_ave, G_err);

  return tuple_ave_err;
}

// average Zbil
tuple<vvvvd_t, vvvvd_t> ave_err_Z(vector<jZbil_t> jZ)
{
  int _bilmoms = (int)jZ.size();
  int _nbil = (int)jZ[0].size();
  int _njacks = (int)jZ[0][0].size();
  int _nmr = (int)jZ[0][0][0].size();

  vvvvd_t Z_ave(vvvd_t(vvd_t(vd_t(0.0, _nmr), _nmr), _nbil), _bilmoms);
  vvvvd_t sqr_Z_ave = Z_ave;
  vvvvd_t Z_err = Z_ave;

  for (int imom = 0; imom < _bilmoms; imom++)
  {
#pragma omp parallel for collapse(3)
    for (int ibil = 0; ibil < _nbil; ibil++)
      for (int mrA = 0; mrA < _nmr; mrA++)
        for (int mrB = 0; mrB < _nmr; mrB++)
          for (int ijack = 0; ijack < _njacks; ijack++)
          {
            Z_ave[imom][ibil][mrA][mrB] += jZ[imom][ibil][ijack][mrA][mrB] / _njacks;
            sqr_Z_ave[imom][ibil][mrA][mrB] +=
                (jZ[imom][ibil][ijack][mrA][mrB] * jZ[imom][ibil][ijack][mrA][mrB]) / _njacks;
          }
#pragma omp parallel for collapse(3)
    for (int ibil = 0; ibil < _nbil; ibil++)
      for (int mrA = 0; mrA < _nmr; mrA++)
        for (int mrB = 0; mrB < _nmr; mrB++)
        {
          Z_err[imom][ibil][mrA][mrB] =
              sqrt((double)(njacks - 1)) *
              sqrt(fabs(sqr_Z_ave[imom][ibil][mrA][mrB] -
                        Z_ave[imom][ibil][mrA][mrB] * Z_ave[imom][ibil][mrA][mrB]));
        }
  }

  tuple<vvvvd_t, vvvvd_t> tuple_ave_err(Z_ave, Z_err);

  return tuple_ave_err;
}

// average ratios of Zbil[ibil1]/Zbil[ibil2]
tuple<vvvd_t, vvvd_t> ave_err_Z(vector<jZbil_t> jZ, const int ibil1, const int ibil2)
{
  int _bilmoms = (int)jZ.size();
  int _njacks = (int)jZ[0][0].size();
  int _nmr = (int)jZ[0][0][0].size();

  vvvd_t Z_ave(vvd_t(vd_t(0.0, _nmr), _nmr), _bilmoms);
  vvvd_t sqr_Z_ave = Z_ave;
  vvvd_t Z_err = Z_ave;

  for (int imom = 0; imom < _bilmoms; imom++)
  {
#pragma omp parallel for collapse(2)
    for (int mrA = 0; mrA < _nmr; mrA++)
      for (int mrB = 0; mrB < _nmr; mrB++)
        for (int ijack = 0; ijack < _njacks; ijack++)
        {
          double jZratio = (jZ[imom][ibil1][ijack][mrA][mrB] / jZ[imom][ibil2][ijack][mrA][mrB]);
          Z_ave[imom][mrA][mrB] += jZratio / _njacks;
          sqr_Z_ave[imom][mrA][mrB] += (jZratio * jZratio) / _njacks;
        }
#pragma omp parallel for collapse(2)
    for (int mrA = 0; mrA < _nmr; mrA++)
      for (int mrB = 0; mrB < _nmr; mrB++)
      {
        Z_err[imom][mrA][mrB] =
            sqrt((double)(njacks - 1)) *
            sqrt(fabs(sqr_Z_ave[imom][mrA][mrB] - Z_ave[imom][mrA][mrB] * Z_ave[imom][mrA][mrB]));
      }
  }

  tuple<vvvd_t, vvvd_t> tuple_ave_err(Z_ave, Z_err);

  return tuple_ave_err;
}

// average sigma
tuple<vvvvd_t, vvvvd_t> ave_err(vector<vvvvd_t> sig)
{
  int _linmoms = (int)sig.size();
  int nproj = (int)sig[0].size();
  int nins = (int)sig[0][0].size();
  int _njacks = (int)sig[0][0][0].size();
  int _nmr = (int)sig[0][0][0][0].size();

  vvvvd_t sig_ave(vvvd_t(vvd_t(vd_t(0.0, _nmr), nins), nproj), _linmoms);
  vvvvd_t sqr_sig_ave = sig_ave;
  vvvvd_t sig_err = sig_ave;

#pragma omp parallel for collapse(4)
  for (int imom = 0; imom < _linmoms; imom++)
    for (int iproj = 0; iproj < nproj; iproj++)
      for (int ins = 0; ins < nins; ins++)
        for (int mr = 0; mr < _nmr; mr++)
          for (int ijack = 0; ijack < _njacks; ijack++)
          {
            sig_ave[imom][iproj][ins][mr] += sig[imom][iproj][ins][ijack][mr] / _njacks;
            sqr_sig_ave[imom][iproj][ins][mr] +=
                sig[imom][iproj][ins][ijack][mr] * sig[imom][iproj][ins][ijack][mr] / _njacks;
          }
#pragma omp parallel for collapse(4)
  for (int imom = 0; imom < _linmoms; imom++)
    for (int iproj = 0; iproj < nproj; iproj++)
      for (int ins = 0; ins < nins; ins++)
        for (int mr = 0; mr < _nmr; mr++)
          sig_err[imom][iproj][ins][mr] =
              sqrt((double)(njacks - 1)) *
              sqrt(fabs(sqr_sig_ave[imom][iproj][ins][mr] -
                        sig_ave[imom][iproj][ins][mr] * sig_ave[imom][iproj][ins][mr]));

  tuple<vvvvd_t, vvvvd_t> tuple_ave_err(sig_ave, sig_err);

  return tuple_ave_err;
}

// average Zq
tuple<vvd_t, vvd_t> ave_err_Zq(vector<vvd_t> jZq)
{
  int _moms = (int)jZq.size();
  int _njacks = njacks;
  int _nmr = (int)jZq[0][0].size();

  vvd_t Zq_ave(vd_t(0.0, _nmr), _moms);
  vvd_t sqr_Zq_ave(vd_t(0.0, _nmr), _moms);
  vvd_t Zq_err(vd_t(0.0, _nmr), _moms);

#pragma omp parallel for collapse(2)
  for (int imom = 0; imom < _moms; imom++)
    for (int mr = 0; mr < _nmr; mr++)
      for (int ijack = 0; ijack < _njacks; ijack++)
      {
        Zq_ave[imom][mr] += jZq[imom][ijack][mr] / njacks;
        sqr_Zq_ave[imom][mr] += jZq[imom][ijack][mr] * jZq[imom][ijack][mr] / njacks;
      }
#pragma omp parallel for collapse(2)
  for (int imom = 0; imom < _moms; imom++)
    for (int mr = 0; mr < _nmr; mr++)
      Zq_err[imom][mr] = sqrt((double)(njacks - 1)) *
                         sqrt(fabs(sqr_Zq_ave[imom][mr] - Zq_ave[imom][mr] * Zq_ave[imom][mr]));

  tuple<vvd_t, vvd_t> tuple_ave_err(Zq_ave, Zq_err);

  return tuple_ave_err;
}

// average of valence effective mass and mPCAC
tuple<vvd_t, vvd_t> ave_err(vvvd_t jM)
{
  int _njacks = jM.size();
  int _nm = (int)jM[0].size();

  vvd_t M_ave(vd_t(0.0, _nm), _nm);
  vvd_t sqr_M_ave(vd_t(0.0, _nm), _nm);
  vvd_t M_err(vd_t(0.0, _nm), _nm);

  for (int mA = 0; mA < _nm; mA++)
    for (int mB = 0; mB < _nm; mB++)
      for (int ijack = 0; ijack < _njacks; ijack++)
      {
        M_ave[mA][mB] += jM[ijack][mA][mB] / _njacks;
        sqr_M_ave[mA][mB] += jM[ijack][mA][mB] * jM[ijack][mA][mB] / _njacks;
      }

  for (int mA = 0; mA < _nm; mA++)
    for (int mB = 0; mB < _nm; mB++)
      M_err[mA][mB] = sqrt((double)(njacks - 1)) *
                      sqrt(fabs(sqr_M_ave[mA][mB] - M_ave[mA][mB] * M_ave[mA][mB]));

  tuple<vvd_t, vvd_t> tuple_ave_err(M_ave, M_err);

  return tuple_ave_err;
}

// average of sea effective mass
tuple<double, double> ave_err(vd_t jM)
{
  int _njacks = jM.size();

  double M_ave = 0.0, sqr_M_ave = 0.0, M_err = 0.0;

  for (int ijack = 0; ijack < _njacks; ijack++)
  {
    M_ave += jM[ijack] / _njacks;
    sqr_M_ave += jM[ijack] * jM[ijack] / _njacks;
  }
  M_err = sqrt((double)(njacks - 1)) * sqrt(fabs(sqr_M_ave - M_ave * M_ave));

  tuple<double, double> tuple_ave_err(M_ave, M_err);

  return tuple_ave_err;
}

// average deltam
tuple<vd_t, vd_t> ave_err(vvd_t jdeltam)
{
  int _njacks = njacks;
  int _nmr = (int)jdeltam[0].size();

  vd_t deltam_ave(0.0, _nmr);
  vd_t sqr_deltam_ave(0.0, _nmr);
  vd_t deltam_err(0.0, _nmr);

  for (int mr = 0; mr < _nmr; mr++)
    for (int ijack = 0; ijack < _njacks; ijack++)
    {
      deltam_ave[mr] += jdeltam[ijack][mr] / njacks;
      sqr_deltam_ave[mr] += jdeltam[ijack][mr] * jdeltam[ijack][mr] / njacks;
    }

  for (int mr = 0; mr < _nmr; mr++)
    deltam_err[mr] = sqrt((double)(njacks - 1)) *
                     sqrt(fabs(sqr_deltam_ave[mr] - deltam_ave[mr] * deltam_ave[mr]));

  tuple<vd_t, vd_t> tuple_ave_err(deltam_ave, deltam_err);

  return tuple_ave_err;
}

// average stepscaling function
tuple<vd_t, vd_t> ave_err(vector<vd_t> jy)
{
  int _njacks = njacks;
  int _nn = (int)jy[0].size();

  vd_t y_ave(0.0, _nn);
  vd_t sqr_y_ave(0.0, _nn);
  vd_t y_err(0.0, _nn);

  for (int i = 0; i < _nn; i++)
    for (int ijack = 0; ijack < _njacks; ijack++)
    {
      y_ave[i] += jy[ijack][i] / njacks;
      sqr_y_ave[i] += jy[ijack][i] * jy[ijack][i] / njacks;
    }

  for (int i = 0; i < _nn; i++)
    y_err[i] = sqrt((double)(njacks - 1)) * sqrt(fabs(sqr_y_ave[i] - y_ave[i] * y_ave[i]));

  tuple<vd_t, vd_t> tuple_ave_err(y_ave, y_err);

  return tuple_ave_err;
}