#ifndef _AVE_ERR_HPP
#define _AVE_ERR_HPP

#include "aliases.hpp"
#include <tuple>

// average sigma
tuple<vvvvd_t, vvvvd_t> ave_err(vector<vvvvd_t> sigma);

// average of Zq
tuple<vvd_t, vvd_t> ave_err_Zq(vector<vvd_t> jZq);

// average of bilinears
tuple<vvvvvd_t, vvvvvd_t> ave_err(vector<jproj_t> jG);

// average of Zbil
tuple<vvvvd_t, vvvvd_t> ave_err_Z(vector<jZbil_t> jZ);
tuple<vvvd_t, vvvd_t> ave_err_Z(vector<jZbil_t> jZ, const int ibil1, const int ibil2);

// average of valence effective mass and mPCAC
tuple<vvd_t, vvd_t> ave_err(vvvd_t jM);

// average of sea effective mass
tuple<double, double> ave_err(vd_t jM);

// average of deltam
tuple<vd_t, vd_t> ave_err(vvd_t jdeltam);

// average stepscaling function
tuple<vd_t, vd_t> ave_err(vector<vd_t> jy);

#endif
