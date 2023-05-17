#include "Dirac.hpp"
#include "aliases.hpp"
#include "global.hpp"
#include "operations.hpp"
#include "sigmas.hpp"
#include <omp.h>

using namespace sigma;

void oper_t::compute_Zq()
{
  cout << "Computing Zq" << endl;

#pragma omp parallel for collapse(3)
  for (int ilinmom = 0; ilinmom < _linmoms; ilinmom++)
    for (int ijack = 0; ijack < njacks; ijack++)
      for (int mr = 0; mr < _nmr; mr++)
      {
        // LO
        jZq[ilinmom][ijack][mr] = sigma[ilinmom][SIGMA1][LO][ijack][mr];
      }
}
