#include "aliases.hpp"
#include "global.hpp"
#include "operations.hpp"
#include "sigmas.hpp"
#include "vertices.hpp"
#include <omp.h>

//debug
#include <iostream>

/*****************************************************************/
/** NOTE: including only the factor 1/(16pi^2), and NOT CF*g^2  **/
/*****************************************************************/


double Zq_sub_RIprime_Pla(const double a2p2, const double a4p4, const double CSW, const double beta) //Iwasaki (RI'-MOM)
{
  return -0.00007117634677642765*a2p2 + (0.015227491298399787*a4p4)/a2p2 + 0.00960047769133904*a2p2*beta - (0.0025066440329093664*a4p4*beta)/a2p2 - 0.004455503827550094*a2p2*CSW + (0.0022450517822973216*a4p4*CSW)/a2p2 - 0.003383626323334429*a2p2*CSW*CSW + (0.0008193849167990292*a4p4*CSW*CSW)/a2p2 + 0.0008179574721126228*a2p2*log(a2p2) - (0.005329916431185477*a4p4*log(a2p2))/a2p2 - 0.0023747152416172916*a2p2*beta*log(a2p2) + 0.0015831434944115277*a2p2*CSW*log(a2p2) + 0.0015831434944115277*a2p2*CSW*CSW*log(a2p2);
}
double Zq_sub_Sigma1_Pla(const double a2p2, const double a4p4, const double CSW, const double beta, const double np) //Iwasaki (Sigma1 Variant - eq.33 of 1004.1115)
{
  return -0.00007117634677642765*a2p2 + (0.000897114646833199*a4p4)/a2p2 + 0.00960047769133904*a2p2*beta - (0.0006596431226714699*a4p4*beta)/a2p2 - 0.004455503827550094*a2p2*CSW - 0.003383626323334429*a2p2*CSW*CSW + (0.014330376651566588*a2p2)/np - (0.0018470009102378967*a2p2*beta)/np + (0.0022450517822973216*a2p2*CSW)/np + (0.0008193849167990292*a2p2*CSW*CSW)/np + 0.0008179574721126228*a2p2*log(a2p2) - 0.0023747152416172916*a2p2*beta*log(a2p2) + 0.0015831434944115277*a2p2*CSW*log(a2p2) + 0.0015831434944115277*a2p2*CSW*CSW*log(a2p2) - (0.005329916431185477*a2p2*log(a2p2))/np;
}
double Zq_sub_RIprime_Iwa(const double a2p2, const double a4p4, const double CSW, const double beta) //Iwasaki (RI'-MOM)
{
  return -0.005672860990323958*a2p2 + (0.014243712884525125*a4p4)/a2p2 + 0.00960047769133904*a2p2*beta - (0.0025066440329093664*a4p4*beta)/a2p2 - 0.003520117837775507*a2p2*CSW + (0.0011673178475060664*a4p4*CSW)/a2p2 - 0.0026499567501341303*a2p2*CSW*CSW + (0.0007616610125116473*a4p4*CSW*CSW)/a2p2 + 0.0019009331651897352*a2p2*log(a2p2) - (0.006098479826272459*a4p4*log(a2p2))/a2p2 -  0.0023747152416172916*a2p2*beta*log(a2p2) + 0.0015831434944115277*a2p2*CSW*log(a2p2) + 0.0015831434944115277*a2p2*CSW*CSW*log(a2p2);
}
double Zq_sub_Sigma1_Iwa(const double a2p2, const double a4p4, const double CSW, const double beta, const double np) //Iwasaki (Sigma1 Variant - eq.33 of 1004.1115)
{
  return -0.005672860990323958*a2p2 + (0.0006875064481731128*a4p4)/a2p2 + 0.00960047769133904*a2p2*beta - (0.0006596431226714699*a4p4*beta)/a2p2 - 0.003520117837775507*a2p2*CSW - 0.0026499567501341303*a2p2*CSW*CSW + (0.013556206436352016*a2p2)/np - (0.0018470009102378967*a2p2*beta)/np + (0.0011673178475060664*a2p2*CSW)/np + (0.0007616610125116473*a2p2*CSW*CSW)/np +  0.0019009331651897352*a2p2*log(a2p2) - 0.0023747152416172916*a2p2*beta*log(a2p2) + 0.0015831434944115277*a2p2*CSW*log(a2p2) + 0.0015831434944115277*a2p2*CSW*CSW*log(a2p2) - (0.006098479826272459*a2p2*log(a2p2))/np;
}

double G_sub_Iwa(const double a2p2, const double a4p4, const double CSW, const double beta, const int ibil) //Iwasaki
{
  /***************************/
  /** {"S","V","P","A","T"} **/
  /***************************/

  double G_sub[5]={0.0};

  /*S*/
  G_sub[0]=-0.0015980657612147536*a2p2 + (0.001590531497385448*a4p4)/a2p2 + 0.014397673260269254*a2p2*beta + (0.0007915717472057639*a4p4*beta)/a2p2 - 0.020483315101767596*a2p2*CSW + 0.001485005896759095*a2p2*CSW*CSW + 0.0031662869888230555*a2p2*log(a2p2) - 0.004749430483234583*a2p2*beta*log(a2p2) + 0.009498860966469166*a2p2*CSW*log(a2p2);

  /*V*/
  G_sub[1]=-0.0027276751592846564*a2p2 - (0.0022818374899451492*a4p4)/a2p2 + 0.004552325318635962*a2p2*beta - (0.0009894646840072048*a4p4*beta)/a2p2 + 0.005120829084735341*a2p2*CSW - 0.0017542110700737687*a2p2*CSW*CSW - 0.00013931662750821467*a2p2*log(a2p2) - 0.0011873576208086458*a2p2*beta*log(a2p2) - 0.0023747152416172916*a2p2*CSW*log(a2p2) + 0.0023747152416172916*a2p2*CSW*CSW*log(a2p2);

  /*P*/
  G_sub[2]=-0.0010650965927161962*a2p2 + (0.001590531497385448*a4p4)/a2p2 + 0.005307337913079719*a2p2*beta + (0.0007915717472057639*a4p4*beta)/a2p2 - 0.0019137787941134038*a2p2*CSW*CSW - 0.0015831434944115277*a2p2*beta*log(a2p2);

  /*A*/
  G_sub[3]=-0.002994159312760922*a2p2 - (0.0022818374899451492*a4p4)/a2p2 + 0.009097492928904988*a2p2*beta - (0.0009894646840072048*a4p4*beta)/a2p2 - 0.005120828921699463*a2p2*CSW + 0.0017542111427847394*a2p2*CSW*CSW + 0.001443826866903313*a2p2*log(a2p2) - 0.0027705011152201735*a2p2*beta*log(a2p2) + 0.0023747152416172916*a2p2*CSW*log(a2p2) - 0.0023747152416172916*a2p2*CSW*CSW*log(a2p2);

  /*T*/
  G_sub[4]=-0.003370695667882934*a2p2 - (0.0035726271523886808*a4p4)/a2p2 + 0.005815710297525325*a2p2*beta - (0.0015831434944115277*a4p4*beta)/a2p2 + 0.003413886048322409*a2p2*CSW + 0.00007146234983073258*a2p2*CSW*CSW + 0.00034195899479289027*a2p2*log(a2p2) - 0.0015831434944115277*a2p2*beta*log(a2p2) - 0.0015831434944115277*a2p2*CSW*log(a2p2);

  return G_sub[ibil];
}

double G_sub_Pla(const double a2p2, const double a4p4, const double CSW, const double beta, const int ibil) //Plaquette
{
  /***************************/
  /** {"S","V","P","A","T"} **/
  /***************************/

  double G_sub[5]={0.0};

  /*S*/
  G_sub[0]=0.0020696349336202964*a2p2 + (0.002638572490685879*a4p4)/a2p2 + 0.014397673260269254*a2p2*beta + (0.0007915717472057639*a4p4*beta)/a2p2 - 0.026285269017811638*a2p2*CSW + 0.001117135138087462*a2p2*CSW*CSW + 0.0031662869888230555*a2p2*log(a2p2) - 0.004749430483234583*a2p2*beta*log(a2p2) + 0.009498860966469166*a2p2*CSW*log(a2p2);

  /*V*/
  G_sub[1]=0.0015162493926228803*a2p2 - (0.0015831434944115277*a4p4)/a2p2 + 0.004552325318635962*a2p2*beta - (0.0009894646840072048*a4p4*beta)/a2p2 + 0.006571317085934714*a2p2*CSW - 0.002473231798463179*a2p2*CSW*CSW - 0.0011873576208086458*a2p2*log(a2p2) - 0.0011873576208086458*a2p2*beta*log(a2p2) - 0.0023747152416172916*a2p2*CSW*log(a2p2) + 0.0023747152416172916*a2p2*CSW*CSW*log(a2p2);

  /*P*/
  G_sub[2]=0.0000386171004030586*a2p2 + (0.002638572490685879*a4p4)/a2p2 + 0.005307337913079719*a2p2*beta + (0.0007915717472057639*a4p4*beta)/a2p2 - 0.0016352912863466967*a2p2*CSW*CSW - 0.0015831434944115277*a2p2*beta*log(a2p2);

  /*A*/
  G_sub[3]=0.0025317582724788627*a2p2 - (0.0015831434944115277*a4p4)/a2p2 + 0.009097492928904988*a2p2*beta - (0.0009894646840072048*a4p4*beta)/a2p2 - 0.0065713173417203145*a2p2*CSW + 0.002473231795251119*a2p2*CSW*CSW + 0.00039578587360288194*a2p2*log(a2p2) - 0.0027705011152201735*a2p2*beta*log(a2p2) + 0.0023747152416172916*a2p2*CSW*log(a2p2) - 0.0023747152416172916*a2p2*CSW*CSW*log(a2p2);

  /*T*/
  G_sub[4]=0.0023472959538066015*a2p2 - (0.0029903821561106632*a4p4)/a2p2 + 0.005815710297525325*a2p2*beta - (0.0015831434944115277*a4p4*beta)/a2p2 + 0.004380878043315581*a2p2*CSW + 0.0000863593230192642*a2p2*CSW*CSW - 0.0010554289962743518*a2p2*log(a2p2) - 0.0015831434944115277*a2p2*beta*log(a2p2) - 0.0015831434944115277*a2p2*CSW*log(a2p2);

  return G_sub[ibil];
}

double subSigma1(const int imom, const vector<double> Np, const double CSW, const vector<double> a2p2_vec, const vector<double> a4p4_vec, const bool QED)
{
  bool QCD=!QED;

  double a2p2 = a2p2_vec[imom];
  double a4p4 = a4p4_vec[imom];
  double np = Np[imom];

  double Zq_sub=0.0;

  double beta;
  if(strcmp(an_suffix.c_str(),"")==0)
  beta = 0.0; //Feynman
  else
  beta = 1.0; //Landau


  if(free_analysis) //pure QED
  {
    if(QCD)
    {
      Zq_sub = 0.0; // no subtraction in free theory
    }
    if(QED)
    {
      if(!UseSigma1)
      Zq_sub = Zq_sub_RIprime_Pla(a2p2,a4p4,CSW,beta);
      else
      Zq_sub = Zq_sub_Sigma1_Pla(a2p2,a4p4,CSW,beta,np);
    }
  }
  else if(inte_analysis) //QCD+QED  (Iwasaki)
  {
    if(QCD) //pure QCD
    {
      if(!UseSigma1)
      Zq_sub = Zq_sub_RIprime_Iwa(a2p2,a4p4,CSW,beta);
      else
      Zq_sub = Zq_sub_Sigma1_Iwa(a2p2,a4p4,CSW,beta,np);
    }
    if(QED) // QCD+QED (subtract only QED, because this is the relative EM correction)
    {
      if(!UseSigma1)
      Zq_sub = Zq_sub_RIprime_Pla(a2p2,a4p4,CSW,beta);
      else
      Zq_sub = Zq_sub_Sigma1_Pla(a2p2,a4p4,CSW,beta,np);
    }
  }

  return Zq_sub;
}

double subG(const int imom, const double CSW, const vector<double> a2p2_vec, const vector<double> a4p4_vec, const int ibil, const bool QED)
{
  bool QCD=!QED;

  double a2p2 = a2p2_vec[imom];
  double a4p4 = a4p4_vec[imom];

  double G_sub=0.0;

  double beta;
  if(strcmp(an_suffix.c_str(),"")==0)
  beta = 0.0; //Feynman
  else
  beta = 1.0; //Landau

  if(free_analysis) //pure QED
  {
    if(QCD)
    {
      G_sub = 0.0; // no subtraction in free theory
    }
    if(QED)
    {
      G_sub = G_sub_Pla(a2p2,a4p4,CSW,beta,ibil);
    }
  }
  else if(inte_analysis) //QCD+QED  (Iwasaki)
  {
    if(QCD) //pure QCD
    {
      G_sub = G_sub_Iwa(a2p2,a4p4,CSW,beta,ibil);
    }
    if(QED) // QCD+QED (subtract only QED, because this is the relative EM correction)
    {
      G_sub = G_sub_Pla(a2p2,a4p4,CSW,beta,ibil);
    }
  }

  return G_sub;

}


oper_t oper_t::subOa2(const int b)
{
  cout<<endl;
  cout<<"----- subtraction of O(a2) effects -----"<<endl<<endl;

  oper_t out=(*this);

  if(!load_ave)
  {

    double CF = 4.0/3.0;
    double g2 = 6.0/beta[b];
    if(sub_boosted)
    g2 /= plaquette;

    double CSW;
    if(csw[b]==0.0)
    CSW=0.0;
    else
    CSW=1.0; //tree-level value

    vector<double> P2(out._linmoms,0.0);
    vector<double> P4(out._linmoms,0.0);
    if(sub_ptilde)
    {P2=p2_tilde; P4=p4_tilde;}
    else
    {P2=p2; P4=p4;}

    // Zq
    for(int imom=0;imom<out._linmoms;imom++)
    {
      double sub = CF*g2*subSigma1(imom,Np,CSW,P2,P4,0);

      #pragma omp parallel for collapse(2)
      for(int ijack=0;ijack<njacks;ijack++)
      for(int mr=0;mr<out._nmr;mr++)
      {
        (out.sigma)[imom][sigma::SIGMA1][sigma::LO][ijack][mr] -= sub;
        //                (out.sigma)[imom][sigma::SIGMA1][sigma::QED][ijack][mr] -= subSigma1(imom,Np,CSW,P2,P4,1);
      }
    }

    out.compute_Zq();

    // Zbil
    for(int imom=0;imom<out._bilmoms;imom++)
    for(int ibil=0;ibil<nbil;ibil++)
    {
      double sub = CF*g2*subG(imom,CSW,P2,P4,ibil,0);

      #pragma omp parallel for collapse(3)
      for(int ijack=0;ijack<njacks;ijack++)
      for(int mr1=0;mr1<out._nmr;mr1++)
      for(int mr2=0;mr2<out._nmr;mr2++)
      {
        (out.jG)[imom][gbil::LO][ibil][ijack][mr1][mr2] -= sub;
        //                        (out.jG)[imom][gbil::QED][ibil][ijack][mr1][mr2] -= subG(imom,CSW,P2,P4,ibil,1);
      }
    }

    out.compute_Zbil();

    //    // ZV/ZA and ZP/ZS
    //    for(int imom=0;imom<out._bilmoms;imom++)
    //        for(int ijack=0;ijack<njacks;ijack++)
    //            for(int mr1=0;mr1<out._nmr;mr1++)
    //                for(int mr2=0;mr2<out._nmr;mr2++)
    //                {
    //                    //                    (out.jZVoverZA)[imom][0][ijack][mr1][mr2] = jZVoverZA[imom][0][ijack][mr1][mr2];
    //                    //                    (out.jZPoverZS)[imom][0][ijack][mr1][mr2] = jZPoverZS[imom][0][ijack][mr1][mr2];
    //                    (out.jZVoverZA)[imom][0][ijack][mr1][mr2] = (out.jZ)[imom][1][ijack][mr1][mr2]/(out.jZ)[imom][3][ijack][mr1][mr2];
    //                    (out.jZPoverZS)[imom][0][ijack][mr1][mr2] = (out.jZ)[imom][2][ijack][mr1][mr2]/(out.jZ)[imom][0][ijack][mr1][mr2];
    //                }


  }

  return out;
}

oper_t oper_t::subOainf(const int b)
{
  cout<<endl;
  cout<<"----- subtraction of O(ainf) effects -----"<<endl<<endl;

  oper_t out=(*this);

  if(!load_ave)
  {
    double CF = 4.0/3.0;
    double g2 = 6.0/beta[b];
    if(sub_boosted)
    g2 /= plaquette;

    vector<vector<double>> DeltaG_moms(out._linmoms,vector<double>(6));

    vector<ifstream> DeltaG_file(6);

    string path = out.path_to_beta+out._beta_label+"/";

    string RCs[6] = {"q","S","P","V","A","T"};
    for(int iRC=0;iRC<6;iRC++)
    {
      DeltaG_file[iRC].open(path+"subainf/DeltaG"+RCs[iRC]);
      if(!DeltaG_file[iRC].good())
      {
        cerr<<"Error opening \""<<path<<"subainf/DeltaG"<<RCs[iRC]<<"\"."<<endl;
        exit(1);
      }
      while(!DeltaG_file[iRC].eof())
      {
        for(int imom=0; imom<out._linmoms; imom++)
        DeltaG_file[iRC]>>DeltaG_moms[imom][iRC];
      }
    }

    // Zq
    for(int imom=0;imom<out._linmoms;imom++)
    {
      double sub = CF*g2*DeltaG_moms[imom][0];

      #pragma omp parallel for collapse(2)
      for(int ijack=0;ijack<njacks;ijack++)
      for(int mr=0;mr<out._nmr;mr++)
      {
        (out.sigma)[imom][sigma::SIGMA1][sigma::LO][ijack][mr] -= sub;
      }
    }
    out.compute_Zq();

    // Zbil
    for(int imom=0;imom<out._bilmoms;imom++)
    for(int ibil=0;ibil<nbil;ibil++)
    {
      double sub = CF*g2*DeltaG_moms[imom][ibil+1];

      #pragma omp parallel for collapse(3)
      for(int ijack=0;ijack<njacks;ijack++)
      for(int mr1=0;mr1<out._nmr;mr1++)
      for(int mr2=0;mr2<out._nmr;mr2++)
      {
        (out.jG)[imom][gbil::LO][ibil][ijack][mr1][mr2] -= sub;
      }
    }

    out.compute_Zbil();

  }

  return out;
}
