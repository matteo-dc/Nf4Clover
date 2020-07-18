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

        printf( "---------------------------------------------\n");
        printf( " ██████╗ ██╗       ██████╗ ███████╗██████╗ \n");
        printf( " ██╔══██╗██║      ██╔═══██╗██╔════╝██╔══██╗\n");
        printf( " ██████╔╝██║█████╗██║   ██║█████╗  ██║  ██║\n");
        printf( " ██╔══██╗██║╚════╝██║▄▄ ██║██╔══╝  ██║  ██║\n");
        printf( " ██║  ██║██║      ╚██████╔╝███████╗██████╔╝\n");
        printf( " ╚═╝  ╚═╝╚═╝       ╚══▀▀═╝ ╚══════╝╚═════╝ \n");
        printf( "--------------------------------------------\n\n\n");

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
            voper_t M1(ntheta);
            voper_t M2(ntheta);

            for(int th=0; th<ntheta; th++)
            {

                for(int m=0; m<nm_Sea[b]; m++)
                {
                    /*  basic  */

                    basic[th][m].create_basic(b,th,m);
                    if(!load_ave)  basic[th][m].plot("");

                    if(!only_basic)
                    {

                        /* ///////  ENRICO  /////// */

                        /*  average r  */

                       rave[th][m] = basic[th][m].average_r();
                       if(!load_ave) rave[th][m].plot("rave");

                       /*  no sub  */

                       // sub[th][m] = rave[th][m];
                       // if(!load_ave) sub[th][m].plot("sub_nosub");

                       /*  perturbative subtraction of O(a2)  */

                       // sub[th][m] = rave[th][m].subOa2(b);
                       // if(!load_ave) sub[th][m].plot("sub_a2");

                       /*  perturbative subtraction of O(ainf)  */

                       sub[th][m] = rave[th][m].subOainf(b,"allmoms");
                       if(!load_ave) sub[th][m].plot("sub_ainf");

                       /*  democratic filter on momenta  */

                       filt[th][m] = sub[th][m].filter_moms();
                       if(!load_ave)  filt[th][m].plot("filt");

                       /*  average over equivalent momenta  */

                       ave[th][m] = filt[th][m].average_equiv_moms();
                       if(!load_ave) ave[th][m].plot("ave");

                       /* store averaged ingredients */

                        if(!load_ave) ave[th][m].print("ave_Enr028");

                        /* load averaged ingredients if needed */

                        if(load_ave) ave[th][m].load("ave_Enr028");


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
//                            oper_for_eta[b][th][m][loop] = val_chir[th][m];
                            exit(0);
                        }
                        else
                        {
//                            evo[th][m] = val_chir[th][m].evolve(ainv[b],p2ref);
//                            evo[th][m].plot("evo");
//                            sub[th][m] = evo[th][m].subOa2(b);
//                            sub[th][m].plot("sub");


//                            M1[th][b][m] = sub[th][m].a2p2_extr();
//                            M1[th][b][m].plot("M1");
//                            M2[th][b][m] = sub[th][m].interpolate_to_p2ref(b);
//                           /* M2[th][b][m] = sub[th][m].interpolate_to_ainv(b); */
//                            M2[th][b][m].plot("M2");
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

                        M1[th] = evo[th].a2p2_extr(ainv[b]);
                        M1[th].plot("M1");

                        M2[th] = evo[th].interpolate_to_p2ref(b);
                        M2[th].plot("M2");

                    }
                }

                exit(0);

            } //close ntheta

        } //close nbeta

    } //close nloop



    /*   DA MODIFICARE DA QUI IN POI  */

//    if(eta_analysis and !only_basic)
//    {
//        for(int b=0; b<nbeta; b++)
//        {
//            for(int th=0; th<ntheta; th++)
//            {
//                eta[b][th].resize(nm_Sea[b]);
//                etaM1[th][b].resize(nm_Sea[b]);
//                etaM2[th][b].resize(nm_Sea[b]);
//
//                for(int m=0; m<nm_Sea[b]; m++)
//                {
//                    eta[b][th][m] = compute_eta(oper_for_eta[b][th][m]);
////#warning uncorrelated photons
////                    eta[b][th][m] = compute_eta_uncorr(oper_for_eta[b][th][m],oper_for_eta[b][th][(m+1)%nm_Sea[b]]);
//                    eta[b][th][m].plot("eta");
//
//                    eta[b][th][m] = eta[b][th][m].evolve_mixed(ainv[b]);
//                    eta[b][th][m].plot("eta_evo");
//
//                    etaM2[th][b][m] = eta[b][th][m].interpolate_to_p2ref(b);
////                    etaM2[th][b][m] = eta[b][th][m].interpolate_to_ainv(b);
//                    etaM2[th][b][m].plot("eta_M2");
//                    etaM1[th][b][m] = eta[b][th][m].a2p2_extr();
//                    etaM1[th][b][m].plot("eta_M1");
//
//                } //close nm_sea
//            } //close theta
//        } //close beta
//
//        if(nm_Sea_max>1)
//            for(int th=0; th<ntheta; th++)
//            {
//                etaM2_sea[th] = combined_chiral_sea_extr(etaM2[th]);
//                etaM1_sea[th] = combined_chiral_sea_extr(etaM1[th]);
//
//                for(int b=0; b<nbeta; b++)
//                {
//                    etaM2_sea[th][b].plot("eta_M2_sea");
//                    etaM1_sea[th][b].plot("eta_M1_sea");
//                }
//            }
//
//        if(ntheta>1)
//        {
//            etaM2_theta = theta_average(etaM2_sea);
//            etaM1_theta = theta_average(etaM1_sea);
//
//            for(int b=0; b<nbeta; b++)
//            {
//                etaM2_theta[b].plot("eta_M2_theta");
//                etaM1_theta[b].plot("eta_M1_theta");
//            }
//        }
//    }
//    else if((free_analysis or inte_analysis) and !only_basic)
//    {
//        if(nm_Sea_max>1)
//            for(int th=0; th<ntheta; th++)
//            {
//                M2_sea[th] = combined_chiral_sea_extr(M2[th]);
//                M1_sea[th] = combined_chiral_sea_extr(M1[th]);
//
//                for(int b=0; b<nbeta; b++)
//                {
//                    M2_sea[th][b].plot("M2_sea");
//                    M1_sea[th][b].plot("M1_sea");
//                }
//            }
//
//        if(ntheta>1)
//        {
//            M2_theta = theta_average(M2_sea);
//            M1_theta = theta_average(M1_sea);
//
//            for(int b=0; b<nbeta; b++)
//            {
//                M2_theta[b].plot("M2_theta");
//                M1_theta[b].plot("M1_theta");
//            }
//        }
//    }

//    if(nbeta>1 and only_basic==0)
//    {
//        /*  a2p2->0 extrapolation  */
//
//        rave_ave_chir_sub_sea_theta_evo_extr = a2p2_extr(rave_ave_chir_sub_sea_theta_evo);
//
//        /*  continuum limit  */
//
//        cout<<"Continuum limit extrapolation:"<<endl<<endl;
////        a2p2_extr(rave_chir_sub_sea_theta_evo_ave,LO);
////        a2p2_extr(rave_chir_sub_sea_theta_evo_ave,EM);
//    }

    cout<<"DONE!"<<endl;

    return 0;
}
