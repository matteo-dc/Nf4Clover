
//**  Evolution to the scale 1/a from Nuria  **//

#include "global.hpp"
#include "aliases.hpp"
#include "operations.hpp"

#define Z3 1.2020569031595942

double alphas(int Nf,double mu2,int nloop) // mu2 is in physical units (dimensional quantity)
{

    double CF = (Nc*Nc-1.)/(2.*Nc);

    double lam0,L2,LL2,b1,b2,b3;
    double als0, als1, als2, als3;
    double beta_0,beta_1,beta_2,beta_3;

    beta_0 = (11.*Nc-2.*Nf)/3.;

    beta_1 = 34./3.*pow(Nc,2) - 10./3.*Nc*Nf-2*CF*Nf;
    beta_2 = (2857./54.)*pow(Nc,3) + pow(CF,2)*Nf -
    205./18.*CF*Nc*Nf -1415./54.*pow(Nc,2)*Nf +
    11./9.*CF*pow(Nf,2) + 79./54.*Nc*pow(Nf,2);

    beta_3 = (150653./486. - 44./9.*Z3)*pow(Nc,4) +
    (-39143./162. + 68./3.*Z3)*pow(Nc,3)*Nf +
    (7073./486. - 328./9.*Z3)*CF*pow(Nc,2)*Nf +
    (-2102./27. + 176./9.*Z3)*pow(CF,2)*Nc*Nf +
    23.*pow(CF,3)*Nf + (3965./162. + 56./9.*Z3)*pow(Nc,2)*pow(Nf,2) +
    (338./27. - 176./9.*Z3)*pow(CF,2)*pow(Nf,2) +
    (4288./243. + 112./9.*Z3)*CF*Nc*pow(Nf,2) + 53./243.*Nc*pow(Nf,3) +
    154./243.*CF*pow(Nf,3) +
    (-10./27. + 88./9.*Z3)*pow(Nc,2)*(pow(Nc,2)+36.) +
    (32./27. - 104./9.*Z3)*Nc*(pow(Nc,2)+6)*Nf +
    (-22./27. + 16./9.*Z3)*(pow(Nc,4) - 6.*pow(Nc,2) + 18.)/pow(Nc,2)*pow(Nf,2);

    b1=beta_1/beta_0/4./M_PI;
    b2=beta_2/beta_0/16./pow(M_PI,2);
    b3=beta_3/beta_0/64./pow(M_PI,3);

    lam0=LambdaQCD;

    L2   = log( mu2/(pow(lam0,2) ) );
    LL2  = log( L2 );

    als0 = 4.*M_PI/beta_0/L2;
    als1 = als0 - pow(als0,2)*b1*LL2;
    als2 = als1 + pow(als0,3)*(pow(b1,2)*(pow(LL2,2) - LL2 -1.) + b2);
    if(nloop==3)
        als3 = als2 + pow(als0,4)*(pow(b1,3)*(-pow(LL2,3)+5./2.*pow(LL2,2)+2*LL2-1./2.)-
                                   3.*b1*b2*LL2 + b3/2.);
    else
        als3 = als2;

    return als3;
}


///////////////////////////////////
// evolution from mu=p to mu0
// Z(mu0)=Z(mu) c(mu0)/c(mu)
// def: c=c(mu)/c(mu0)
// -> Z(mu0) = Z(mu) /c
//////////////////////////////////

double q_evolution_to_RIp(int Nf,double p2,double p2ref)
{
    double cmu=0.0, cmu0=0.0; // c=cmu/cmu0
    //mu_2 = p2     [from] (physical)
    //mu0_2= p2ref  [to]

    // alphas @ NNLO
    double alm, al0;
    alm=alphas(Nf,p2,3)/(4*M_PI);
    al0=alphas(Nf,p2ref,3)/(4*M_PI);

    ////////////////////////////////
    // N3LO FORMULA
    // Assuming landau gauge
    ///////////////////////////////////
    if(Nf==2){
        cmu = 1. + 2.03448 * alm + 35.9579 * pow(alm,2) + 1199.16 * pow(alm,3);
        cmu0 = 1. + 2.03448 * al0 + 35.9579 * pow(al0,2) + 1199.16 * pow(al0,3);
    }if(Nf==0){
        cmu = 1. + 2.0303 * alm + 42.1268 * pow(alm,2) + 1728.43 * pow(alm,3);
        cmu0 = 1. + 2.0303 * al0 + 42.1268 * pow(al0,2) + 1728.43 * pow(al0,3);
    }if(Nf==4){
        cmu = 1. + 2.4000 * alm + 29.6724 * pow(alm,2) + 719.141 * pow(alm,3);
        cmu0 = 1. + 2.4000 * al0 + 29.6724 * pow(al0,2) + 719.141 * pow(al0,3);
    }

    return cmu/cmu0;
}
double S_evolution_to_RIp(int Nf,double p2,double p2ref)
{
    double cmu=0.0, cmu0=0.0; // c=cmu/cmu0
    //mu_2 = p2     [from] (physical)
    //mu0_2= p2ref  [to]

    // alphas @ NNLO
    double alm, al0;
    alm=alphas(Nf,p2,3)/(4*M_PI);
    al0=alphas(Nf,p2ref,3)/(4*M_PI);

    ////////////////////////////////
    // N3LO FORMULA
    // Assuming landau gauge
    ///////////////////////////////////

    if(Nf==2){
        cmu = pow(alm,-12./29) * (1. - 8.55727 * alm - 125.423 * pow(alm,2) -
                                  3797.71 * pow(alm,3));

        cmu0 = pow(al0,-12./29) * (1. - 8.55727 * al0 - 125.423 * pow(al0,2) -
                                   3797.71 * pow(al0,3));
    }if(Nf==0){
        cmu = pow(alm,-4./11) * (1. - 8.08264 * alm - 151.012 * pow(alm,2) -
                                 5247.93 * pow(alm,3));

        cmu0 = pow(al0,-4./11) * (1. - 8.08264 * al0 - 151.012 * pow(al0,2) -
                                  5247.93 * pow(al0,3));
    }if(Nf==4){
        cmu = pow(alm,-12./25) * (1. - 9.38987 * alm - 96.2883 * pow(alm,2) -
                                  2403.82 * pow(alm,3));

        cmu0 = pow(al0,-12./25) * (1. - 9.38987 * al0 - 96.2883 * pow(al0,2) -
                                   2403.82 * pow(al0,3));
    }

    return cmu/cmu0;
}

double P_evolution_to_RIp(int Nf,double p2,double p2ref)
{
    double cmu=0.0, cmu0=0.0; // c=cmu/cmu0
    //mu_2 = p2     [from] (physical)
    //mu0_2= p2ref  [to]

    // alphas @ NNLO
    double alm, al0;
    alm=alphas(Nf,p2,3)/(4*M_PI);
    al0=alphas(Nf,p2ref,3)/(4*M_PI);

    ////////////////////////////////
    // N3LO FORMULA
    // Assuming landau gauge
    ///////////////////////////////////
    if(Nf==2){
        cmu = pow(alm,-12./29) * (1. - 8.55727 * alm - 125.423 * pow(alm,2) -
                                  3797.71 * pow(alm,3));

        cmu0 = pow(al0,-12./29) * (1. - 8.55727 * al0 - 125.423 * pow(al0,2) -
                                   3797.71 * pow(al0,3));
    }if(Nf==0){
        cmu = pow(alm,-4./11) * (1. - 8.08264 * alm - 151.012 * pow(alm,2) -
                                 5247.93 * pow(alm,3));

        cmu0 = pow(al0,-4./11) * (1. - 8.08264 * al0 - 151.012 * pow(al0,2) -
                                  5247.93 * pow(al0,3));
    }if(Nf==4){
        cmu = pow(alm,-12./25) * (1. - 9.38987 * alm - 96.2883 * pow(alm,2) -
                                  2403.82 * pow(alm,3));

        cmu0 = pow(al0,-12./25) * (1. - 9.38987 * al0 - 96.2883 * pow(al0,2) -
                                   2403.82 * pow(al0,3));
    }

    return cmu/cmu0;
}

double T_evolution_to_RIp(int Nf,double p2,double p2ref)
{
    double cmu=0.0, cmu0=0.0; // c=cmu/cmu0
    //mu_2 = p2     [from] (physical)
    //mu0_2= p2ref  [to]

    // alphas @ NNLO
    double alm, al0;
    alm=alphas(Nf,p2,3)/(4*M_PI);
    al0=alphas(Nf,p2ref,3)/(4*M_PI);

    ////////////////////////////////
    // N2LO FORMULA
    // Assuming landau gauge
    ///////////////////////////////////

    if(Nf==2){
        cmu = pow(alm,4./29) * (1. + 2.66852 * alm + 47.9701 * pow(alm,2) + 1612.4556 * pow(alm,3));

        cmu0 = pow(al0,4./29) * (1. + 2.66852 * al0 + 47.9701 * pow(al0,2) + 1612.4556 * pow(al0,3));
    }if(Nf==0){
        cmu = pow(alm,4./33) * (1. + 2.53260 * alm + 57.8740 * pow(alm,2) + 2275.2426 * pow(alm,3));

        cmu0 = pow(al0,4./33) * (1. + 2.53260 * al0 + 57.8740 * pow(al0,2) + 2275.2426 * pow(al0,3));
    }if(Nf==4){
        cmu = pow(alm,4./25) * (1. + 2.91662 * alm + 37.9471 * pow(alm,2) + 993.7488 * pow(alm,3));

        cmu0 = pow(al0,4./25) * (1. + 2.91662 * al0 + 37.9471 * pow(al0,2) + 993.7488 * pow(al0,3));
    }

    return cmu/cmu0;
}

pair<int,int> find_stepfunc_imom(const int size, const vector<double> a2p2, const double _ainv)
{
  vector<pair<double,int>> dist_list_min, dist_list_max;
  // find momentum closer to p2step_min & p2step_max
  for(int j=0;j<size;j++)
  {
    dist_list_min.push_back(make_pair( fabs(a2p2[j]*pow(_ainv,2.0) - stepfunc_min) , j ));
    dist_list_max.push_back(make_pair( fabs(a2p2[j]*pow(_ainv,2.0) - stepfunc_max) , j ));
  }
  sort(dist_list_min.begin(), dist_list_min.end());
  sort(dist_list_max.begin(), dist_list_max.end());

  return make_pair(dist_list_min[0].second,dist_list_max[0].second);
}

oper_t oper_t::evolve(const double _ainv, const double p2_evol)
{
    // p2_evol in physical units

    if(!inte_analysis)
    {
        cout<<"Function 'evolve' only works in 'inte' mode."<<endl;
        exit(0);
    }

    cout<<endl;
    cout<<"----- evolution to the scale p2 = "<<p2_evol<<" GeV^2 -----"<<endl<<endl;

    oper_t out=(*this);

    double cq=0.0, cq_step=0.0;
    vd_t cO(0.0,5), cO_step(0.0,5);

    // find p2 (physical units) for step scaling function
    cout<<"  combinations for step scaling function:"<<endl;
    pair<int,int> jminmax = find_stepfunc_imom(_linmoms,p2,_ainv);
    int step_imom_min = jminmax.first;
    int step_imom_max = jminmax.second;

    cout<<"   p2min["<<step_imom_min<<"]="<<p2[step_imom_min]*pow(_ainv,2.0)<<"  ("<<stepfunc_min<<")"<<endl;
    cout<<"   p2max["<<step_imom_max<<"]="<<p2[step_imom_max]*pow(_ainv,2.0)<<"  ("<<stepfunc_max<<")"<<endl<<endl;

    for(int ilinmom=0;ilinmom<_linmoms;ilinmom++)
    {
        cq=q_evolution_to_RIp(Nf,p2[ilinmom]*pow(_ainv,2.0),p2_evol);

        for(int ijack=0;ijack<njacks;ijack++)
          for(int mr1=0; mr1<_nmr; mr1++)
          {
              (out.jZq)[ilinmom][ijack][mr1] = jZq[ilinmom][ijack][mr1]/cq;
          }
    }

    cq_step=q_evolution_to_RIp(Nf,p2[step_imom_min]*pow(_ainv,2.0),p2[step_imom_max]*pow(_ainv,2.0));
    /**/
    cout<<"     Sigma(q) [pert.] = "<<1.0/cq_step<<endl;
    /**/
    for(int ijack=0;ijack<njacks;ijack++)
    {
      (out.stepfunc)[ijack][0] = jZq[step_imom_max][ijack][0]/jZq[step_imom_min][ijack][0];
      (out.stepfunc)[ijack][0] -= 1.0/cq_step;
    }

    for(int ibilmom=0;ibilmom<_bilmoms;ibilmom++)
    {
        // Note that ZV  ZA are RGI because they're protected by the WIs
        cO[0]=S_evolution_to_RIp(Nf,p2[ibilmom]*pow(_ainv,2.0),p2_evol); //S
        cO[1]=1.0;                                       //A
        cO[2]=P_evolution_to_RIp(Nf,p2[ibilmom]*pow(_ainv,2.0),p2_evol); //P
        cO[3]=1.0;                                       //V
        cO[4]=T_evolution_to_RIp(Nf,p2[ibilmom]*pow(_ainv,2.0),p2_evol); //T

        for(int ibil=0;ibil<5;ibil++)
            for(int ijack=0;ijack<njacks;ijack++)
                for(int mr1=0; mr1<_nmr; mr1++)
                    for(int mr2=0; mr2<_nmr; mr2++)
                    {
                        (out.jZ)[ibilmom][ibil][ijack][mr1][mr2] = jZ[ibilmom][ibil][ijack][mr1][mr2]/cO[ibil];
                        // (out.jZ_EM)[ibilmom][ibil][ijack][mr1][mr2] = jZ_EM[ibilmom][ibil][ijack][mr1][mr2]/cO[ibil];
                    }
    }

    cO_step[0]=S_evolution_to_RIp(Nf,p2[step_imom_min]*pow(_ainv,2.0),p2[step_imom_max]*pow(_ainv,2.0)); //S
    cO_step[1]=1.0;                                       //A
    cO_step[2]=P_evolution_to_RIp(Nf,p2[step_imom_min]*pow(_ainv,2.0),p2[step_imom_max]*pow(_ainv,2.0)); //P
    cO_step[3]=1.0;                                       //V
    cO_step[4]=T_evolution_to_RIp(Nf,p2[step_imom_min]*pow(_ainv,2.0),p2[step_imom_max]*pow(_ainv,2.0)); //T
    /**/
    vector<string> bil={"S","V","P","A","T"};
    for(int i=0;i<5;i++)
      cout<<"     Sigma("<<bil[i]<<") [pert.] = "<<1.0/cO_step[i]<<endl;
    cout<<endl;
    /**/
    for(int ijack=0;ijack<njacks;ijack++)
      for(int ibil=0;ibil<5;ibil++)
      {
        (out.stepfunc)[ijack][ibil+1] = jZ[step_imom_max][ibil][ijack][0][0]/jZ[step_imom_min][ibil][ijack][0][0];
        (out.stepfunc)[ijack][ibil+1] -= 1.0/cO_step[ibil];
      }

    for(int ijack=0;ijack<njacks;ijack++)
    {
      (out.stepfunc)[ijack][nbil+1] = jZVoverZA[step_imom_max][0][ijack][0][0]/jZVoverZA[step_imom_min][0][ijack][0][0];
      (out.stepfunc)[ijack][nbil+1] -= 1.0;
      (out.stepfunc)[ijack][nbil+2] = jZPoverZS[step_imom_max][0][ijack][0][0]/jZPoverZS[step_imom_min][0][ijack][0][0];
      (out.stepfunc)[ijack][nbil+2] -= 1.0;
      (out.stepfunc)[ijack][nbil+3] = jZAoverZV[step_imom_max][0][ijack][0][0]/jZAoverZV[step_imom_min][0][ijack][0][0];
      (out.stepfunc)[ijack][nbil+3] -= 1.0;
    }

    return out;


}