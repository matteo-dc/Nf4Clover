#include "aliases.hpp"
#include "global.hpp"
#include "Dirac.hpp"
#include <stdio.h>
#include <omp.h>
#include "operations.hpp"
#include "vertices.hpp"
#include "Zbil.hpp"

void oper_t::compute_Zbil()
{
    cout<<"Computing Zbil"<<endl;

    for(int ibilmom=0;ibilmom<_bilmoms;ibilmom++)
    {
//        const int imom1=bilmoms[ibilmom][1]; // p1
//        const int imom2=bilmoms[ibilmom][2]; // p2
        const int imom1=ibilmom;
        const int imom2=ibilmom;

//compute Z's according to 'riqed.pdf', one for each momentum
#pragma omp parallel for collapse(3)
      for(int mr_fw=0;mr_fw<_nmr;mr_fw++)
        for(int mr_bw=0;mr_bw<_nmr;mr_bw++)
          for(int ijack=0;ijack<njacks;ijack++)
          {
            for(int ibil=0;ibil<nbil;ibil++)
            {
              // LO
              jZ[ibilmom][ibil][ijack][mr_fw][mr_bw] =
                sqrt(jZq[imom1][ijack][mr_fw]*jZq[imom2][ijack][mr_bw])/
                jG[ibilmom][gbil::LO][ibil][ijack][mr_fw][mr_bw];
            }
            jZVoverZA[ibilmom][0][ijack][mr_fw][mr_bw]=
              jZ[ibilmom][1][ijack][mr_fw][mr_bw]/ jZ[ibilmom][3][ijack][mr_fw][mr_bw];
            jZPoverZS[ibilmom][0][ijack][mr_fw][mr_bw]=
              jZ[ibilmom][2][ijack][mr_fw][mr_bw]/ jZ[ibilmom][0][ijack][mr_fw][mr_bw];
          }

    }// close mom loop

}
