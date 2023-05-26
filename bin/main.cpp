#ifdef HAVE_CONFIG_H
#include <config.hpp>
#endif

#include "Dirac.hpp"
#include "aliases.hpp"
#include "fit.hpp"
#include "global.hpp"
#include "operations.hpp"
#include "print.hpp"
#include "read.hpp"
#include "read_input.hpp"
#include <array>
#include <fstream>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

int main(int narg, char **arg)
{
  omp_set_nested(1);
#pragma omp parallel
#pragma omp master
  {
    system("clear");

    printf("------------------------------------------------------------------------------\n");
    printf(" ███╗   ██╗███████╗██╗  ██╗ ██████╗██╗      ██████╗ ██╗   ██╗███████╗██████╗ \n");
    printf(" ████╗  ██║██╔════╝██║  ██║██╔════╝██║     ██╔═══██╗██║   ██║██╔════╝██╔══██╗ \n");
    printf(" ██╔██╗ ██║█████╗  ███████║██║     ██║     ██║   ██║██║   ██║█████╗  ██████╔╝ \n");
    printf(" ██║╚██╗██║██╔══╝  ╚════██║██║     ██║     ██║   ██║╚██╗ ██╔╝██╔══╝  ██╔══██╗ \n");
    printf(" ██║ ╚████║██║          ██║╚██████╗███████╗╚██████╔╝ ╚████╔╝ ███████╗██║  ██║ \n");
    printf(" ╚═╝  ╚═══╝╚═╝          ╚═╝ ╚═════╝╚══════╝ ╚═════╝   ╚═══╝  ╚══════╝╚═╝  ╚═╝\n");
    printf("------------------------------------------------------------------------------\n\n\n");

    cout << "Using " << omp_get_num_threads() << " threads" << endl << endl;
  }

  char path_glb[128] = "input_glb.txt";

  if (narg > 2)
  {
    cerr << "Number of arguments not valid." << endl;
    exit(0);
  }
  if (narg == 2)
  {
    string path = arg[1];
    strcpy(path_glb, path.c_str());
  }

  cout << "Reading global input from \"" << path_glb << "\"." << endl;
  read_input_glb(path_glb);

  voper_t all(nbeta);

  for (int b = 0; b < nbeta; b++)
  {
    voper_t basic(nm_Sea[b]), rave(nm_Sea[b]), sub(nm_Sea[b]), filt(nm_Sea[b]), ave(nm_Sea[b]),
        val_chir(nm_Sea[b]);

    oper_t sea_chir, evo, cont;
    oper_t M1, M2, M1b, M2b, M3, M4;

    for (int m = 0; m < nm_Sea[b]; m++)
    {

      /*  basic  */

      basic[m].create_basic(b, m);
      if (!load_ave)
        basic[m].plot("");

      if (!only_basic)
      {

        /*  average r  */

        rave[m] = basic[m].average_r();
        if (!load_ave)
          rave[m].plot("rave");

        /* perturbative subtraction of lattice artifacts */
        if (subtraction == 0)
        {
          /*  no sub  */
          sub[m] = rave[m];
          if (!load_ave)
            sub[m].plot("sub_nosub");
        }
        else if (subtraction == 1)
        {
          /*  perturbative subtraction of O(a2)  */
          sub[m] = rave[m].subOa2(b);
          if (!load_ave)
            sub[m].plot("sub_a2");
        }
        else if (subtraction == 2)
        {
          /*  perturbative subtraction of O(ainf)  */
          sub[m] = rave[m].subOainf(b, "allmoms");
          if (!load_ave)
            sub[m].plot("sub_ainf");

          /*  perturbative subtraction of O(ainf) à la Martha  */

          // sub[m] = rave[m].subOainfMartha(b,"allmoms");
          // sub[m].plot("sub_ainf_Martha");
        }

        /*  democratic filter on momenta  */

        filt[m] = sub[m].filter_moms();
        if (!load_ave)
          filt[m].plot("filt");

        /*  average over equivalent momenta  */

        ave[m] = filt[m].average_equiv_moms();
        if (!load_ave)
          ave[m].plot("ave");

        /* store averaged ingredients */

        if (!load_ave)
          ave[m].print(load_label);
        // if(!load_ave) ave[m].print("ave_mdc");   // 04/22 [new projectors]
        // if(!load_ave) ave[m].print("ave_Enr028");    // latest working version: 10/03/22
        // if(!load_ave) ave[m].print("ave_Martha");

        /* load averaged ingredients if needed */

        if (load_ave)
        {
          ave[m].load(load_label);
          ave[m].plot("ave");
        }

        // if(load_ave) ave[m].load("ave_Enr028");  // latest working version: 10/03/22
        // if(load_ave) ave[m].load("ave_Martha");

        /*  valence chiral extr  */
        val_chir[m] = ave[m].chiral_extr();
        val_chir[m].plot("chir");

        /* store extrapolated ingredients */
        if (!load_chir)
          val_chir[m].printZ("chir");

      } // close if(!only_basic)

    } // close nm_sea

    if (!only_basic)
    {
      sea_chir = chiral_sea_extr(val_chir);
      sea_chir.plot("sea");

      evo = sea_chir.evolve(ainv[b], p2ref);
      evo.plot("evo");

      // save all the curves for combined fit
      all[b] = evo;

      // linear extrapolation
      M1 = evo.a2p2_extr(b);
      M1.plot("M1");

      // constant interpolation
      M2 = evo.interpolate_to_p2ref(b);
      M2.plot("M2");

      // linear + pole 1/p2 extrapolation
      M3 = evo.a2p2_extr_with_pole(b);
      M3.plot("M3");

      // cout<<"-----------ZA/ZV-----------"<<endl;
      // for(int ijack=0;ijack<njacks;ijack++)
      // {
      //     cout<< (M3.jZ[0][3][ijack][0][0]/M3.jZ[0][1][ijack][0][0])<<endl;
      // }
      // cout<<"--------------------------"<<endl;

      // // quadratic + pole 1/p2 extrapolation
      M4 = evo.a2p2_extr_with_pole_and_p4(b);
      M4.plot("M4");

      // quadratic
      // M4 = evo.a2p2_extr_with_p4(b);
      // M4.plot("M4");

      // /*****/
      //
      // cont = sea_chir.remove_hadr_cont(ainv[b]);
      // cont = cont.evolve(ainv[b],p2ref);
      // cont.plot("cont");
      //
      // M1b = cont.a2p2_extr(ainv[b]);
      // M1b.plot("M1b");
      //
      // M2b = cont.interpolate_to_p2ref(b);
      // M2b.plot("M2b");

      /* Method 1 : Z-improvement */
      cont = evo.Z_improvement(ainv[b]);
      cont.plot("a1");
    }

  } // close nbeta

  voper_t M3combo(nbeta);
  voper_t M4combo(nbeta);
  voper_t M5combo(nbeta);
  voper_t M5_quadcombo(nbeta);
  voper_t M5_logcombo(nbeta);

  if (!only_basic && nbeta == 3)
  {
    M3combo = combined_M3(all);
    for (int b = 0; b < nbeta; b++)
      M3combo[b].plot("M3comb");

    M4combo = combined_M4(all);
    for (int b = 0; b < nbeta; b++)
      M4combo[b].plot("M4comb");

    M5combo = combined_M5(all);
    for (int b = 0; b < nbeta; b++)
      M5combo[b].plot("M5comb");

    M5_quadcombo = combined_M5_quadratic(all);
    for (int b = 0; b < nbeta; b++)
      M5_quadcombo[b].plot("M5quadcomb");

    M5_logcombo = combined_M5_log(all);
    for (int b = 0; b < nbeta; b++)
      M5_logcombo[b].plot("M5logcomb");
  }

  cout << "DONE!" << endl;

  return 0;
}
