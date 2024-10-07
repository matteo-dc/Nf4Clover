#include "Dirac.hpp"
#include "aliases.hpp"
#include "contractions.hpp"
#include "global.hpp"
#include "jack.hpp"
#include "operations.hpp"
#include "print.hpp"
#include "prop.hpp"
#include "rotate.hpp"
#include <chrono>
#include <omp.h>

#define EXTERN_VERT
#include "vertices.hpp"

using namespace std::chrono;

namespace gbil
{
void set_ins()
{
  if (ntypes == 1)
  {
    ins_list = {LO};
    ins_tag = {"LO"};
  }

  nins = ins_list.size();
  nGamma = 16;
}
} // namespace gbil

// calculate the vertex function in a given configuration for the given equal momenta
prop_t make_vertex(const prop_t &prop1, const prop_t &prop2, const int mu)
{
  return prop1 * GAMMA[mu] * GAMMA[5] * prop2.adjoint() * GAMMA[5];
}

// compute LO vertices
void build_vert(const vvvprop_t &S1, const vvvprop_t &S2, valarray<jvert_t> &jVert)
{
#pragma omp parallel for collapse(4)
  for (int ijack = 0; ijack < njacks; ijack++)
    for (int mr_fw = 0; mr_fw < nmr; mr_fw++)
      for (int mr_bw = 0; mr_bw < nmr; mr_bw++)
        for (int igam = 0; igam < 16; igam++)
        {
          if (ntypes == 1)
          {
            // LO
            jVert[gbil::LO][ijack][mr_fw][mr_bw][igam] +=
                make_vertex(S1[ijack][qprop::LO][mr_fw], S2[ijack][qprop::LO][mr_bw], igam);
          }
        }
}

// project the amputated green function
jproj_t compute_pr_bil(vvvprop_t &jpropOUT_inv, valarray<jvert_t> &jVert, vvvprop_t &jpropIN_inv)
{
  int lambda_size = gbil::nins + 2 * jprop::nins - 2;

  vector<int> i1;
  vector<int> iv;
  vector<int> i2;
  vector<int> ip;

  if (ntypes == 1)
  {
    i1 = {jprop::LO};
    iv = {gbil::LO};
    i2 = {jprop::LO};

    ip = {gbil::LO};
  }

  jproj_t pr_bil(vvvvd_t(vvvd_t(vvd_t(vd_t(0.0, nmr), nmr), njacks), nbil), gbil::nins);

  const int ibil_of_igam[/*gbil::nGamma*/ 16] = {0, 1, 1, 1, 1, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4};

#pragma omp parallel for collapse(5)
  for (int ijack = 0; ijack < njacks; ijack++)
    for (int mr_fw = 0; mr_fw < nmr; mr_fw++)
      for (int mr_bw = 0; mr_bw < nmr; mr_bw++)
        for (int k = 0; k < lambda_size; k++)
          for (int igam = 0; igam < gbil::nGamma; igam++)
          {
            prop_t lambda_igam = jpropOUT_inv[i1[k]][ijack][mr_fw] *
                                 jVert[iv[k]][ijack][mr_fw][mr_bw][igam] * GAMMA[5] *
                                 (jpropIN_inv[i2[k]][ijack][mr_bw]).adjoint() * GAMMA[5];

            /* Usual RI' scheme projectors */
            pr_bil[ip[k]][ibil_of_igam[igam]][ijack][mr_fw][mr_bw] +=
                (lambda_igam * Proj[igam]).trace().real() / 12.0;
          }

  return pr_bil;
}

// project the amputated green function with p-dependent transverse projectors
jproj_t compute_pr_bil(vvvprop_t &jpropOUT_inv, valarray<jvert_t> &jVert, vvvprop_t &jpropIN_inv,
                       p_t p, const double p2)
{
  int lambda_size = gbil::nins + 2 * jprop::nins - 2;

  vector<int> i1;
  vector<int> iv;
  vector<int> i2;
  vector<int> ip;

  if (ntypes == 1)
  {
    i1 = {jprop::LO};
    iv = {gbil::LO};
    i2 = {jprop::LO};

    ip = {gbil::LO};
  }

  jproj_t pr_bil(vvvvd_t(vvvd_t(vvd_t(vd_t(0.0, nmr), nmr), njacks), nbil), gbil::nins);

  const int ibil_of_igam[/*gbil::nGamma*/ 16] = {0, 1, 1, 1, 1, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4};

#pragma omp parallel for collapse(5)
  for (int ijack = 0; ijack < njacks; ijack++)
    for (int mr_fw = 0; mr_fw < nmr; mr_fw++)
      for (int mr_bw = 0; mr_bw < nmr; mr_bw++)
        for (int k = 0; k < lambda_size; k++)
          for (int igam = 0; igam < gbil::nGamma; igam++)
          {
            prop_t lambda_igam = jpropOUT_inv[i1[k]][ijack][mr_fw] *
                                 jVert[iv[k]][ijack][mr_fw][mr_bw][igam] * GAMMA[5] *
                                 (jpropIN_inv[i2[k]][ijack][mr_bw]).adjoint() * GAMMA[5];

            if (ibil_of_igam[igam] == 1) /* V */
            {
              // igam = {1,2,3,4}
              int ig = igam; // ig = {1,2,3,4}

              /* Projector for RI" scheme - compatible with WI: change in vector and axial current
               */
              // for (int ig2 = 1; ig2 <= 4; ig2++)
              // {
              //   pr_bil[ip[k]][ibil_of_igam[igam]][ijack][mr_fw][mr_bw] +=
              //       (lambda_igam * Proj[ig2] * (4.0 / 3.0) *
              //        ((ig == ig2 ? 1.0 : 0.0) - p[ig % 4] * p[ig2 % 4] / p2))
              //           .trace()
              //           .real() /
              //       12.0;
              // }
              pr_bil[ip[k]][ibil_of_igam[igam]][ijack][mr_fw][mr_bw] +=
                  (lambda_igam * GAMMA[ig]).trace().real();
              for (int ig2 = 1; ig2 <= 4; ig2++)
              {
                pr_bil[ip[k]][ibil_of_igam[igam]][ijack][mr_fw][mr_bw] +=
                    -0.25 * (lambda_igam * GAMMA[ig2] * p[ig2 % 4] / p[ig % 4]).trace().real();
              }
              pr_bil[ip[k]][ibil_of_igam[igam]][ijack][mr_fw][mr_bw] /= 36.0;
            }
            else if (ibil_of_igam[igam] == 3) /* A */
            {
              // igam = {6,7,8,9}
              int ig = igam - 5; // ig = {1,2,3,4}

              /* Projector for RI" scheme - compatible with WI: change in vector and axial current
               */
              for (int ig2 = 1; ig2 <= 4; ig2++)
              {
                pr_bil[ip[k]][ibil_of_igam[igam]][ijack][mr_fw][mr_bw] +=
                    (lambda_igam * Proj[ig2 + 5] * (4.0 / 3.0) *
                     ((ig == ig2 ? 1.0 : 0.0) - p[ig % 4] * p[ig2 % 4] / p2))
                        .trace()
                        .real() /
                    12.0;
              }
            }
            else if (ibil_of_igam[igam] == 4) /* T */
            {
              // igam = {10,...,15}
              int ig = igam - 9; // ig = {1,...,6}

              size_t ind1[6] = {4, 4, 4, 2, 3, 1};
              size_t ind2[6] = {1, 2, 3, 3, 1, 2};

              /* Projector for RI" scheme */
              for (int ig2 = 1; ig2 <= 6; ig2++)
              {
                pr_bil[ip[k]][ibil_of_igam[igam]][ijack][mr_fw][mr_bw] +=
                    (lambda_igam * Proj[ig2 + 9] * (6.0 / 3.0) *
                     (

                         p2 * (ind1[ig - 1] == ind1[ig2 - 1] ? 1.0 : 0.0) *
                             (ind2[ig - 1] == ind2[ig2 - 1] ? 1.0 : 0.0) -
                         p[ind1[ig2 - 1] % 4] *
                             ((ind2[ig - 1] == ind2[ig2 - 1] ? 1.0 : 0.0) * p[ind1[ig - 1] % 4] -
                              (ind1[ig - 1] == ind2[ig2 - 1] ? 1.0 : 0.0) * p[ind2[ig - 1] % 4]) +
                         p[ind2[ig2 - 1] % 4] *
                             ((ind2[ig - 1] == ind1[ig2 - 1] ? 1.0 : 0.0) * p[ind1[ig - 1] % 4] -
                              (ind1[ig - 1] == ind1[ig2 - 1] ? 1.0 : 0.0) * p[ind2[ig - 1] % 4])) /
                     p2

                     )
                        .trace()
                        .real() /
                    12.0;
              }
            }
            else
            {
              /* Usual RI' scheme projectors */
              pr_bil[ip[k]][ibil_of_igam[igam]][ijack][mr_fw][mr_bw] +=
                  (lambda_igam * Proj[igam]).trace().real() / 12.0;
            }
          }

  return pr_bil;
}

void oper_t::compute_bil()
{
  using namespace gbil;

  ifstream jG_data(path_print + "jG");

  if (jG_data.good())
  {
    cout << "Reading bilinears from files: \"" << path_print << "jG\"" << endl << endl;
    read_vec_bin(jG, path_print + "jG");
  }
  else
  {
    cout << "Creating the vertices -- ";

    // array of input files to be read in a given conf
    FILE *input[combo];

    const vector<string> v_path = setup_read_qprop(input);

    for (int ibilmom = 0; ibilmom < _bilmoms; ibilmom++)
    {
      high_resolution_clock::time_point t0 = high_resolution_clock::now();

      cout << endl;
      cout << "\r\t bilmom = " << ibilmom + 1 << "/" << _bilmoms << endl;

      const int imom1 = bilmoms[ibilmom][1]; // p1
      const int imom2 = bilmoms[ibilmom][2]; // p2
      const bool read2 = (imom1 != imom2);

      // definition of jackknifed propagators
      /* prop1 */
      jprop_t jS1(vvprop_t(vprop_t(prop_t::Zero(), _nmr), njacks), jprop::nins);
      /* prop2 */
      jprop_t jS2(vvprop_t(vprop_t(prop_t::Zero(), _nmr), njacks), jprop::nins);

      // definition of jackknifed vertices
      valarray<jvert_t> jVert(
          jvert_t(vvvprop_t(vvprop_t(vprop_t(prop_t::Zero(), gbil::nGamma), _nmr), _nmr), njacks),
          gbil::nins);

      cout << "- Building vertices" << endl;

      double t_span1 = 0.0, t_span2 = 0.0, t_span3 = 0.0;

      for (int i_in_clust = 0; i_in_clust < clust_size; i_in_clust++)
        for (int ihit = 0; ihit < nhits; ihit++)
        {
          const int mom1 = linmoms[imom1][0];
          const int mom2 = linmoms[imom2][0];

          high_resolution_clock::time_point ta = high_resolution_clock::now();

          vvvprop_t S1 = read_qprop_mom(input, v_path, i_in_clust, ihit, mom1);
          vvvprop_t S2 = (read2) ? read_qprop_mom(input, v_path, i_in_clust, ihit, mom2) : S1;

          S1 = rotate(S1);
          S2 = (read2) ? rotate(S2) : S1;

          high_resolution_clock::time_point tb = high_resolution_clock::now();
          t_span1 += (duration_cast<duration<double>>(tb - ta)).count();

          ta = high_resolution_clock::now();

          build_prop(S1, jS1);
          if (read2)
            build_prop(S2, jS2);
          else
          {
            jS2 = jS1;
          }

          tb = high_resolution_clock::now();
          t_span2 += (duration_cast<duration<double>>(tb - ta)).count();

          ta = high_resolution_clock::now();

          build_vert(S1, S2, jVert);

          tb = high_resolution_clock::now();
          t_span3 += (duration_cast<duration<double>>(tb - ta)).count();
        }
      cout << "\t read: " << t_span1 << " s" << endl;
      cout << "\t build prop: " << t_span2 << " s" << endl;
      cout << "\t build vert: " << t_span3 << " s" << endl;

      cout << "- Jackknife of propagators and vertices" << endl;

      // jackknife averages
      /* prop1 */
      for (auto &prop1 : jS1)
        prop1 = jackknife(prop1);
      /* prop2 */
      if (read2)
        for (auto &prop2 : jS2)
          prop2 = jackknife(prop2);
      else
        jS2 = jS1;
      /* vert */
      for (int ins = 0; ins < nins; ins++)
        jVert[ins] = jackknife(jVert[ins]);

      cout << "- Inverting propagators" << endl;

      // definition of inverse propagators
      jprop_t jS1_inv(vvprop_t(vprop_t(prop_t::Zero(), _nmr), njacks), jprop::nins);
      jprop_t jS2_inv(vvprop_t(vprop_t(prop_t::Zero(), _nmr), njacks), jprop::nins);

      // invert propagators
      /* prop1 */
      jS1_inv[jprop::LO] = invert_jprop(jS1[jprop::LO]);
      if (jprop::nins != 1)
        for (int i = 1; i < jprop::nins; i++)
          jS1_inv[i] = -jS1_inv[jprop::LO] * jS1[i] * jS1_inv[jprop::LO];
      /* prop2 */
      if (read2)
      {
        jS2_inv[jprop::LO] = invert_jprop(jS2[jprop::LO]);
        if (jprop::nins != 1)
          for (int i = 1; i < jprop::nins; i++)
            jS2_inv[i] = -jS2_inv[jprop::LO] * jS2[i] * jS2_inv[jprop::LO];
      }
      else
        jS2_inv = jS1_inv;

      cout << "- Computing bilinears" << endl;

      // compute the projected green function (S,V,P,A,T)
      // jG[ibilmom] = compute_pr_bil(jS1_inv,jVert,jS2_inv);

      // using the RI" projectors for vector and axial currents
      jG[ibilmom] = compute_pr_bil(jS1_inv, jVert, jS2_inv, p_tilde[ibilmom], p2_tilde[ibilmom]);

      high_resolution_clock::time_point t1 = high_resolution_clock::now();
      duration<double> t_span = duration_cast<duration<double>>(t1 - t0);
      cout << "\t\t time: " << t_span.count() << " s" << endl;

    } // close mom loop
    cout << endl << endl;

    print_vec_bin(jG, path_print + "jG");
  }
}
