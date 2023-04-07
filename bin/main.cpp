#ifdef HAVE_CONFIG_H
#include <config.hpp>
#endif

#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "aliases.hpp"
#include "read.hpp"
#include "read_input.hpp"
#include "Dirac.hpp"
#include "global.hpp"
#include "operations.hpp"
#include "print.hpp"
#include "fit.hpp"

int main(int narg,char **arg)
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

        cout<<"Using "<<omp_get_num_threads()<<" threads"<<endl<<endl;
    }

    char path_glb[128]="input_glb.txt";

    if (narg>2){cerr<<"Number of arguments not valid."<<endl;
                exit(0);}
    if (narg==2){
        string path=arg[1];
        strcpy(path_glb,path.c_str());
    }

    cout<<"Reading global input from \""<<path_glb<<"\"."<<endl;
    read_input_glb(path_glb);

    int nloop;
    if(free_analysis or inte_analysis)
        nloop=1;
    if(eta_analysis)
        nloop=2;

    voper_t all(nbeta); // ONLY WORKS WITH NTHETA=1

    recompute_basic = false;

    for(int loop=0;loop<nloop;loop++)
    {
        if(nloop>1)
            cout<<" ********** LOOP "<<loop<<" ********** "<<endl;


        if(loop>0) recompute_basic=true;

        //////////

        for(int b=0; b<nbeta; b++)
        {
            vvoper_t basic(voper_t(nm_Sea[b]),ntheta);
            vvoper_t filt(voper_t(nm_Sea[b]),ntheta);
            vvoper_t ave(voper_t(nm_Sea[b]),ntheta);

            vvoper_t rave(voper_t(nm_Sea[b]),ntheta);

            vvoper_t val_chir(voper_t(nm_Sea[b]),ntheta);
            vvoper_t sub(voper_t(nm_Sea[b]),ntheta);

            voper_t sea_chir(ntheta);
            voper_t evo(ntheta);
            voper_t cont(ntheta);
            voper_t M1(ntheta);
            voper_t M2(ntheta);
            voper_t M1b(ntheta);
            voper_t M2b(ntheta);
            voper_t M3(ntheta);
            voper_t M4(ntheta);


            for(int th=0; th<ntheta; th++)
            {

                for(int m=0; m<nm_Sea[b]; m++)
                {
                    /*  basic  */

                    basic[th][m].create_basic(b,th,m);
                    if(!load_ave)  basic[th][m].plot("");

                    if(!only_basic)
                    {

                        /*  average r  */

                        rave[th][m] = basic[th][m].average_r();
                        if(!load_ave) rave[th][m].plot("rave");


                        /*  no sub  */

                        if(subtraction == 0)
                        {
                            sub[th][m] = rave[th][m];
                            if(!load_ave) sub[th][m].plot("sub_nosub");
                        }

                        /*  perturbative subtraction of O(a2)  */

                        if(subtraction == 1)
                        {
                            sub[th][m] = rave[th][m].subOa2(b);
                            if(!load_ave) sub[th][m].plot("sub_a2");
                        }

                        /*  perturbative subtraction of O(ainf)  */

                        if(subtraction == 2)
                        {
                            sub[th][m] = rave[th][m].subOainf(b,"allmoms");
                            if(!load_ave) sub[th][m].plot("sub_ainf");

                            /*  perturbative subtraction of O(ainf) à la Martha  */

                            // sub[th][m] = rave[th][m].subOainfMartha(b,"allmoms");
                            // sub[th][m].plot("sub_ainf_Martha");
                        }

                        /*  democratic filter on momenta  */

                        filt[th][m] = sub[th][m].filter_moms();
                        if(!load_ave)  filt[th][m].plot("filt");

                        /*  average over equivalent momenta  */

                        ave[th][m] = filt[th][m].average_equiv_moms();
                        if(!load_ave) ave[th][m].plot("ave");

                        /* store averaged ingredients */

                        if(!load_ave) ave[th][m].print(load_label);
                        // if(!load_ave) ave[th][m].print("ave_mdc");        // printing on other file:    04/22 [new projectors]
                        // if(!load_ave) ave[th][m].print("ave_Enr028");  // latest working version: 10/03/22
                        // if(!load_ave) ave[th][m].print("ave_Martha");

                        /* load averaged ingredients if needed */

                        if(load_ave) ave[th][m].load(load_label);
                        // if(load_ave) ave[th][m].load("ave_Enr028");  // latest working version: 10/03/22
                        // if(load_ave) ave[th][m].load("ave_Martha");


                        /////////////////////////////////////////////////////

                        /*  valence chiral extr  */
                        if(free_analysis or recompute_basic)
                        {
                            val_chir[th][m] = ave[th][m];
                            val_chir[th][m].plot("chir");

                            /* store extrapolated ingredients */
                            if(!load_chir) val_chir[th][m].printZ("chir");
                        }
                        else
                        {
                            val_chir[th][m] = ave[th][m].chiral_extr();
                            val_chir[th][m].plot("chir");

                            /* store extrapolated ingredients */
                            if(!load_chir) val_chir[th][m].printZ("chir");
                        }


                        if(eta_analysis)
                        {
                           // oper_for_eta[b][th][m][loop] = val_chir[th][m];
                            exit(0);
                        }
                        else
                        {
                          //
                        }

                    } //close if(!only_basic)

                } //close nm_sea

                if(!only_basic)
                {
                    if(eta_analysis)
                        exit(0);
                    else
                    {
                        sea_chir[th] = chiral_sea_extr(val_chir[th]);
                        sea_chir[th].plot("sea");

                        evo[th] = sea_chir[th].evolve(ainv[b],p2ref);
                        evo[th].plot("evo");

                        // save all the curves for combined fit (ONLY WORKS WITH NTHETA=1)
                        if(ntheta==1)
                        all[b] = evo[th];

                        // linear extrapolation
                        M1[th] = evo[th].a2p2_extr(b);
                        M1[th].plot("M1");

                        // constant interpolation
                        M2[th] = evo[th].interpolate_to_p2ref(b);
                        M2[th].plot("M2");

                        // linear + pole 1/p2 extrapolation
                        M3[th] = evo[th].a2p2_extr_with_pole(b);
                        M3[th].plot("M3");

                        // cout<<"-----------ZA/ZV-----------"<<endl;
                        // for(int ijack=0;ijack<njacks;ijack++)
                        // {
                        //     cout<< (M3[th].jZ[0][3][ijack][0][0]/M3[th].jZ[0][1][ijack][0][0])<<endl;
                        // }
                        // cout<<"--------------------------"<<endl;

                        // // quadratic + pole 1/p2 extrapolation
                        M4[th] = evo[th].a2p2_extr_with_pole_and_p4(b);
                        M4[th].plot("M4");

                        // quadratic
                        // M4[th] = evo[th].a2p2_extr_with_p4(b);
                        // M4[th].plot("M4");

                        // /*****/
                        //
                        // cont[th] = sea_chir[th].remove_hadr_cont(ainv[b]);
                        // cont[th] = cont[th].evolve(ainv[b],p2ref);
                        // cont[th].plot("cont");
                        //
                        // M1b[th] = cont[th].a2p2_extr(ainv[b]);
                        // M1b[th].plot("M1b");
                        //
                        // M2b[th] = cont[th].interpolate_to_p2ref(b);
                        // M2b[th].plot("M2b");

                        /* Method 1 : Z-improvement */
                        cont[th] = evo[th].Z_improvement(ainv[b]);
                        cont[th].plot("a1");

                    }
                }

            } //close ntheta (fake if Nf4Clover)

        } //close nbeta

    } //close nloop

    voper_t M3combo(nbeta);
    voper_t M4combo(nbeta);
    voper_t M5combo(nbeta);
    voper_t M5_quadcombo(nbeta);
    voper_t M5_logcombo(nbeta);

    if(!only_basic && nbeta==3)
    {
      M3combo = combined_M3(all);
      for(int b=0;b<nbeta;b++)
        M3combo[b].plot("M3comb");

      M4combo = combined_M4(all);
      for(int b=0;b<nbeta;b++)
        M4combo[b].plot("M4comb");

     M5combo = combined_M5(all);
      for(int b=0;b<nbeta;b++)
        M5combo[b].plot("M5comb");

     M5_quadcombo = combined_M5_quadratic(all);
      for(int b=0;b<nbeta;b++)
        M5_quadcombo[b].plot("M5quadcomb");        
    
     M5_logcombo = combined_M5_log(all);
      for(int b=0;b<nbeta;b++)
        M5_logcombo[b].plot("M5logcomb");
    }



    cout<<"DONE!"<<endl;

    return 0;
}
