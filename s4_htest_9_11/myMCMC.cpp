/*
 * main.cpp
 *
 *  Created on: Oct 15, 2018
 *      Author: zhangxu
 */

#include <iostream>
#include <bits/stdc++.h>
//#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <random>
//#include <boost/random.hpp>
//#include <boost/random/normal_distribution.hpp>

#include "Constants.h"
#include "genDRAM_GibbsWithMH.h"

using namespace std;
using namespace boost::numeric::odeint;

//starting values of parameters;
const double Optimize_orig [57]= {Kf_GLC_Transport,
    Kr_GLC_Transport,
    Kf_GLC_SerPool1,
    Kf_GLC_PYR,
    Kf_GLC_PRPP ,
    Kf_SerPool3_Synthesis                   ,
    Kf_SerPool3_Degradation                 ,
    Kf_GlyPool3_Synthesis                   ,
    Kf_GlyPool3_Degradation                 ,
    Kf_SerPool1_GlyPool1                      ,
    Kr_SerPool1_GlyPool1                      ,
    Kf_Ser_Transport                        ,
    Kr_Ser_Transport                        ,
    Kf_Gly_Transport                        ,
    Kr_Gly_Transport                        ,
    Kf_SerPool2_SerPoolMitochon             ,
    Kr_SerPool2_SerPoolMitochon             ,
    Kf_GlyPool2_GlyPoolMitochon             ,
    Kr_GlyPool2_GlyPoolMitochon             ,
    Kf_SerPoolMitochon_GlyPoolMitochon      ,
    Kr_SerPoolMitochon_GlyPoolMitochon      ,
    Kf_GlyPool1_CO2                         ,
    Kr_GlyPool1_CO2                          ,
    Kf_MethyleneTHF_FormylTHF_Cytoplasm     ,
    Kr_MethyleneTHF_FormylTHF_Cytoplasm     ,
    Kf_FormylTHF_Formate_Cytoplasm          ,
    Kr_FormylTHF_Formate_Cytoplasm          ,
    Kf_GlyPoolMitochon_CO2                  ,
    Kr_GlyPoolMitochon_CO2                  ,
    Kf_MethyleneTHF_FormylTHF_Mitochon      ,
    Kr_MethyleneTHF_FormylTHF_Mitochon      ,
    Kf_FormylTHF_Formate_Mitochon           ,
    Kr_FormylTHF_Formate_Mitochon           ,
    Kf_PRPP1_GAR                            ,
    Kf_GAR_FGAR                             ,
    Kf_FGAR_AMP                             ,
    Kf_IMP_AMP                              ,
    Kf_AMP_Degradation                      ,
    Kr_MethyleneTHF_Transport               ,
    Kf_MethyleneTHF_Transport               ,
    Kr_FormylTHF_Transport                  ,
    Kf_FormylTHF_Transport                  ,
    Kr_THF_Transport                        ,
    Kf_THF_Transport                        ,
    Kr_Formate_Transport                    ,
    Kf_Formate_Transport                    ,
    Kf_PRPP2_GAR                            ,
    Kf_PRPP3_GAR                            ,
    Kf_SerPool2_GlyPool2                    ,
    Kr_SerPool2_GlyPool2                    ,
    Kf_GlyPool2_CO2                         ,
    Kr_GlyPool2_CO2                         ,
    CellVolume                                ,
    MediumVolume                            ,
    TissueVolume                            ,
    CytoplasmVolume                         ,
    MitochonVolume                          };

double odeParameters[57] = {0.000461657750477753,0.427503086639990,0.000138933974778344,0.912957765988904,0.000302264234311588,0.000264216988034605,0.528433976069209,0.000178629396143585,0.0105686795213842 ,   389.997313257504 ,   3.65988124056536 ,   0.000422747180855367 ,   0.0496192767065513   , 7.44289150598270e-05 ,   0.0178629396143585  ,  0.0109460124414652 ,   0.0140352064043982 ,   0.00966781987350446  ,  0.0490562874450916  ,  830.998723620212,    17.0631986313531 ,   54.1720082557871,    1.25536770067575 ,   1.67337425755250 ,   2.03451173487994,    0.399124270787164 ,   1000000.00000000 ,   434.000150742915 ,   1.64330339153493  ,  2.22492836752176 ,   0.600557183706194 ,   2.29538774044506 ,   44400.0000000432  ,  89.9889179513641  ,  44000.0000016264 ,   43999.9999979280  ,  0.0266554754467594 ,   0.00151731008995339,    0.189003218774087  ,  0.00751136610783774  ,  11.6715973413360  ,  0.495498097191622 ,   2.51038222698186  ,  0.0677700081258078  ,  4.57755449879482,    0.579744450285708 ,   42.5463885164964 ,   21.0649543138687 ,   85.0037419701686  ,  120.001190410154 ,   65.0980133729533  ,  1.78629396143585 ,   CellVolume  ,  MediumVolume  ,  TissueVolume ,   CytoplasmVolume ,   MitochonVolume};

extern state_type Xd_flag, Xd_flag2;
state_type Xd_flag, Xd_flag2;

const double t0 (0*60), tf (24*60), df (60), MoleculeNumberInOneNanoMole (6.02214129e14);

double mmYobs [3][95];
//={ {3.923781e+00,1.230924e+01,-3.865501e+00,3.946130e+00,7.522404e+00,-5.531175e-01,-2.857281e+00,-5.269439e+00,9.851570e+00,1.949492e+00,-2.350955e+00,-5.069652e+00,-5.535152e+00,-8.901399e+00,-4.642440e+00,3.278378e-01,4.081335e+00,6.944867e-01,1.908361e+00,2.873403e+00,-4.026351e+00,-6.389187e+00,-8.680123e+00,3.739213e+00,-1.469665e+00,2.064665e+00,-5.835460e+00,-9.284215e+00,-1.198053e+01,4.972928e+00,-1.814118e+00,2.506683e+00,-4.353117e+00,-7.220447e+00,-9.799003e+00,2.010781e+00,-3.122531e+00,-3.336353e+00,-5.852487e+00,-3.966557e+00,-6.853213e+00,-6.913261e+00,-9.824434e+00,-1.519986e+00,-5.905339e+00,-1.037060e+01,-6.191903e+00,-1.093310e+01,-8.072083e+00,-1.155661e+01,-5.437412e+00,-1.265746e+01,-4.281763e+00,-8.591180e+00,-1.274291e+01,-1.125580e+01,-1.196219e+00,-5.969395e+00,-4.697448e+00,-8.991246e+00,-1.134333e+01,-1.382845e+01,-1.573229e+01,-1.263374e+01,-1.630341e+01,-1.201151e+00,-4.343192e+00,-5.511824e+00,-8.884382e+00,-4.725345e+00,-7.176104e+00,8.523642e-01,-8.082951e+00,-7.822068e+00,-1.196336e+01,-1.028748e+01,-1.424792e+01,-1.381014e+01,-1.768831e+01,-9.117977e+00,-1.213283e+01,-1.229295e+01,-1.657801e+01,2.187902e+00,-6.864018e-01,-6.772364e-01,-3.619414e+00,-2.559915e+00,-5.532167e+00,-5.855491e+00,-9.063635e+00,-1.754138e+00,-3.969053e+00,-4.152393e+00,-7.403287e+00},
//    {4.248849e+00,1.260609e+01,-4.097609e+00,4.114594e+00,7.027729e+00,-6.217608e-01,-3.293396e+00,-6.031101e+00,8.390275e+00,1.919824e+00,-2.254387e+00,-5.132814e+00,-6.108247e+00,-8.828580e+00,-4.888986e+00,3.088754e-01,4.375514e+00,7.209071e-01,1.753225e+00,2.972377e+00,-4.419424e+00,-6.896750e+00,-8.686345e+00,3.265580e+00,-1.389189e+00,2.162320e+00,-6.126750e+00,-8.606197e+00,-1.210372e+01,4.990319e+00,-2.019362e+00,2.561960e+00,-4.923762e+00,-6.919141e+00,-9.415554e+00,1.963666e+00,-3.120290e+00,-3.317184e+00,-6.612356e+00,-4.381811e+00,-7.019332e+00,-7.663950e+00,-9.564914e+00,-1.263678e+00,-5.687745e+00,-1.135157e+01,-6.724443e+00,-1.053347e+01,-7.975994e+00,-1.219163e+01,-5.372532e+00,-1.227331e+01,-3.718268e+00,-8.550189e+00,-1.287364e+01,-1.384120e+01,-1.218322e+00,-5.573658e+00,-4.707918e+00,-8.952874e+00,-1.254950e+01,-1.403861e+01,-1.562221e+01,-1.324442e+01,-1.535809e+01,-1.179747e+00,-3.957198e+00,-5.584960e+00,-8.399715e+00,-4.609610e+00,-7.079823e+00,7.745343e-01,-7.771338e+00,-8.319837e+00,-1.102579e+01,-1.095932e+01,-1.293240e+01,-1.345523e+01,-1.562549e+01,-9.189701e+00,-1.260788e+01,-1.240774e+01,-1.512089e+01,2.351002e+00,-7.236805e-01,-8.360008e-01,-3.761491e+00,-2.812828e+00,-5.740621e+00,-5.182288e+00,-8.744091e+00,-1.509097e+00,-4.513028e+00,-4.655951e+00,-7.792016e+00},
//    {4.028333e+00,1.197165e+01,-3.902373e+00,4.095526e+00,7.344075e+00,-6.074994e-01,-2.969496e+00,-5.820330e+00,8.888721e+00,1.883694e+00,-2.331839e+00,-4.807703e+00,-5.963726e+00,-8.917468e+00,-4.756403e+00,3.457931e-01,4.202314e+00,7.054040e-01,2.016411e+00,3.058522e+00,-4.395742e+00,-6.224194e+00,-9.325178e+00,3.433540e+00,-1.468458e+00,2.162838e+00,-5.651421e+00,-8.808234e+00,-1.148919e+01,4.951381e+00,-1.807969e+00,2.672061e+00,-4.449637e+00,-6.457612e+00,-9.486866e+00,1.798178e+00,-3.230217e+00,-3.179553e+00,-6.391751e+00,-4.241342e+00,-6.557138e+00,-7.269771e+00,-1.056876e+01,-1.340226e+00,-5.243689e+00,-1.076663e+01,-7.204488e+00,-9.731367e+00,-7.568712e+00,-1.152422e+01,-5.520150e+00,-1.206502e+01,-4.528864e+00,-8.995670e+00,-1.336184e+01,-1.283621e+01,-1.064374e+00,-6.184712e+00,-4.869156e+00,-9.266543e+00,-1.263767e+01,-1.481511e+01,-1.780640e+01,-1.256325e+01,-1.497397e+01,-1.149707e+00,-4.190381e+00,-5.723500e+00,-8.168540e+00,-4.875983e+00,-7.435451e+00,8.195110e-01,-8.383229e+00,-8.277852e+00,-1.047448e+01,-9.911909e+00,-1.329591e+01,-1.317672e+01,-1.661792e+01,-1.009631e+01,-1.157019e+01,-1.155316e+01,-1.578044e+01,2.277610e+00,-7.509357e-01,-7.933384e-01,-3.788405e+00,-2.527732e+00,-6.082214e+00,-5.788140e+00,-8.220406e+00,-1.692411e+00,-4.914480e+00,-4.416052e+00,-7.782693e+00}
//};
double mmYobs2 [3][95];
double Yode [95], Yode2 [95], startp [147], startp2 [147], muprior [147],muprior2 [147], sigmaprior [52];
//double beta_tmp[95], alpha_tmp[95];
int observed [95];
vector<int> ifFit;
vector<double> sd_tune, sd_tune2;

const int nruns (20000);
const  int nburnin (8000);
const int nthinning (1);

const int N_tmp(147);
const int n_Y(95);
//extern const int fitsize(5);

int main(int argc, char **argv)
{
    int i,j; //for loop
    //initial data values
    state_type Xd;
    for(i=0; i<512; ++i) Xd[i] = 0;
    Xd[2] = 0.3 * 6.02214129e17;      // 300 n moles unlabeled cGLC, cGLC_13C0
    Xd[84] = 0.01 * 6.02214129e17; // cSerPool1_13C000D000 = 0.01 * 6.02214129e17; 0.01 u moles unlabeled cSerPool1
    Xd[438] = 0.0001 * 6.02214129e17; //cPRPP_13C0 = 0.0001 * 6.02214129e17; 0.1 n moles unlabeled cPRPP
    Xd[148] = 0.1 * 6.02214129e17; //cGlyPool1_13C00D00 = 0.1 * 6.02214129e17; 0.1 u moles unlabeled cGlyPool1
    Xd[420] = 0.00005 * 6.02214129e17; //cTHF = 0.00005 * 6.02214129e17;      % 0.05 n moles unlabeled cTHF
    Xd[404] = 0.00005 * 6.02214129e17; //cMethyleneTHF_13C0D00 = 0.00005 * 6.02214129e17; 0.05 n moles unlabeled cMethyleneTHF (no carbon label, no deterium label)
    Xd[164] = 0.01 * 6.02214129e17; //cSerPool2_13C000D000 = 0.01 * 6.02214129e17;      % 0.01 u moles unlabeled cSerPool2 (no carbon label, no deterium label)
    Xd[228] = 0.1 * 6.02214129e17; //cGlyPool2_13C00D00 = 0.1 * 6.02214129e17;      % 0.1 u moles unlabeled cGlyPool2 (no carbon label, no deterium label)
    Xd[244] = 0.01 * 6.02214129e17; //cSerPool3_13C000D000 = 0.01 * 6.02214129e17;      % 0.01 u moles unlabeled cSerPool3 (no carbon label, no deterium label)
    Xd[308] = 0.1 * 6.02214129e17; //cGlyPool3_13C00D00 = 0.1 * 6.02214129e17;      % 0.1 u moles unlabeled cGlyPool3 (no carbon label, no deterium label)
    Xd[324] = 0.002 * 6.02214129e17; //mSerPool_13C000D000 = 0.002 * 6.02214129e17;      % 0.002 u moles unlabeled mSerPool (no carbon label, no deterium label)
    Xd[388] = 0.02 * 6.02214129e17; //mGlyPool_13C00D00 = 0.02 * 6.02214129e17;      % 0.02 u moles unlabeled mGlyPool (no carbon label, no deterium label)
    Xd[437] = 0.00005 * 6.02214129e17; //mTHF = 0.00005 * 6.02214129e17;      % 0.05 n moles unlabeled cTHF
    Xd[421] = 0.00005 * 6.02214129e17; //mMethyleneTHF_13C0D00 = 0.00005 * 6.02214129e17;      % 0.05 n moles unlabeled mMethyleneTHF (no carbon label, no deterium label)
    Xd[412] = 0.00005 * 6.02214129e17; //cFormylTHF_13C0D0 = 0.00005 * 6.02214129e17;      % 0.05 n moles unlabeled cFormylTHF (no carbon label, no deterium label)
    Xd[416] = 0.00005 * 6.02214129e17; //cFormate_13C0D0 = 0.00005 * 6.02214129e17;      % 0.05 n moles unlabeled cFormate (no carbon label, no deterium label)
    Xd[433] = 0.00005 * 6.02214129e17; //mFormate_13C0D0 = 0.00005 * 6.02214129e17;      % 0.05 n moles unlabeled mFormate (no carbon label, no deterium label)
    Xd[429] = 0.00005 * 6.02214129e17; //mFormylTHF_13C0D0 = 0.00005 * 6.02214129e17;      % 0.05 n moles unlabeled mFormylTHF (no carbon label, no deterium label)
    Xd[440] = 0.0001 * 6.02214129e17; //cGAR_13C000 = 0.0001 * 6.02214129e17;      % 0.1 n moles unlabeled cGAR (no carbon label, no deterium label)
    Xd[448] = 0.0001 * 6.02214129e17; //cFGAR_13C0000D0 = 0.0001 * 6.02214129e17;      % 0.1 n moles unlabeled cFGAR (no carbon label, no deterium label)
    Xd[480] = 0.02 * 6.02214129e17; //cAMP_13C00000D00 = 0.02 * 6.02214129e17;      % 20 n moles unlabeled cAMP (no carbon label, no deterium label)
    Xd[1] = 236.310 * 6.02214129e17; //eGLC_13C1 = 236.310 * 6.02214129e17;      % 236310 n moles labeled eGLC
    Xd[4] = 2.987748945 * 6.02214129e17; //eSer_13C000D000 = 2.987748945 * 6.02214129e17;      % 2987.748945 n moles unlabeled eSer (no carbon label, no deterium label)
    Xd[68] = 4.496768784 * 6.02214129e17; //eGly_13C00D00 = 4.496768784 * 6.02214129 * 10^17;      % 4496.768784 n moles unlabeled eGly (no carbon label, no deterium label)
    
    
    for(i=0; i<NumberOfIsotopomers; ++i) Xd_flag[i] = Xd[i];//cout<<Xd_flag[i]<<'\t';}
    for(i=0; i<NumberOfIsotopomers; ++i) Xd_flag2[i] = Xd[i];
    Xd_flag2[1] = 196.292 * 6.02214129e17;      // 196292 n moles labeled eGLC;
    Xd_flag2[4] = 3.310533725 * 6.02214129e17;      //3310.533725 n moles unlabeled eSer (no carbon label, no deterium label)
    Xd_flag2[68] = 3.955700675 * 6.02214129e17;      //3955.700675 n moles unlabeled eGly (no carbon label, no deterium label)
    
    //cancer;
    integrate(PurineSynthesis, Xd, t0, tf, df);
    state_type Yd; // the observations at 24h;
    for(i = 0; i < Yd.size() ; ++i ){
        Yd[i] = Xd[i]/MoleculeNumberInOneNanoMole; // unit is nano mole
    }
    
    //exclude unobserved data with 0
    //vector<int> observed;
    j = 0;
    for (i = 0;i <512; ++i){
        if(Yd[i] != 0) {observed[j] = i; ++j;}
    }
    
    for(i=0; i<95; ++i) Yode[i] = log(Yd[observed[i]]);// number of actually observed Y is 95.
    double Y_var[95];//double Y_var [sizeof(Yode)];
    for(i=0; i<n_Y; ++i) Y_var[i] = pow(0.01*Yode[i],2);
    //epsilon ~ MVN(0, SIGMA)
    
    
    //non cancer;
    for(i=0; i<512; ++i) Xd[i] = Xd_flag2[i];
    for(i=0; i<(N_tmp-n_Y); i++) odeParameters[i] = OptimizeParameters2[i];
    integrate(PurineSynthesis, Xd, t0, tf, df);
    
    // the observations at 24h;
    for(i = 0; i < Yd.size() ; ++i ){
        Yd[i] = Xd[i]/MoleculeNumberInOneNanoMole; // unit is nano mole
    }
    
    for(i=0; i<95; ++i) Yode2[i] = log(Yd[observed[i]]);// number of actually observed Y is 95.
    double Y_var2[95];//double Y_var [sizeof(Yode)];
    for(i=0; i<n_Y; ++i) Y_var2[i] = pow(0.01*Yode2[i],2);
    //epsilon ~ MVN(0, SIGMA)
    
    
    
    for(i=0; i<57-5; ++i)
    {
        startp[i] = log(OptimizeParameters[i]);//startp.push_back(log(OptimizeParameters[i]));
        muprior[i] = log(OptimizeParameters[i]);
        startp2[i] = log(OptimizeParameters2[i]);
        muprior2[i] = log(OptimizeParameters2[i]);
    }
    for(i=57-5; i<147; ++i)
    {
        startp[i] = Y_var[i-52];//startp.insert(startp.end(),Y_var.begin(),Y_var.end());
        muprior[i] = Y_var[i-52];
        startp2[i] = Y_var2[i-52];
        muprior2[i] = Y_var2[i-52];
    }
    
    
    double sigma_par [52] = {log(0.01)-log(0.005),0-log(0.5),log(0.1)-log(0.05),0-log(0.5),log(0.01)-log(0.005),
        log(5e-4)-log(2.5e-4),log(1)-log(0.5),log(2e-4)-log(1e-4),log(0.1)-log(0.05),log(500)-log(225),
        log(10)-log(5),log(0.01)-log(0.005),log(0.1)-log(0.05),log(0.001)-log(0.0005),log(0.2)-log(0.1),
        log(0.02)-log(0.01),log(0.03)-log(0.015),log(0.01)-log(0.005),log(0.05)-log(0.025),log(1000)-log(450),
        log(20)-log(10),log(100)-log(45),log(2)-log(1),log(3)-log(1.5),log(10)-log(5),
        log(3)-log(1.5),log(1e7)-log(5e6),log(1e3)-log(500),log(4)-log(2),log(10)-log(5),
        0-log(0.5),log(10)-log(5),log(1e5)-log(5e4),log(100)-log(50),log(2e5)-log(0.9e5),
        log(2e5)-log(0.9e5),log(0.25)-log(0.15),log(0.01)-log(0.005),log(0.5)-log(0.25),log(0.3)-log(0.15),
        log(50)-log(25),log(5)-log(2.5),log(8)-log(4),log(1)-log(0.5),log(10)-log(5),
        log(10)-log(5),log(100)-log(50),log(50)-log(25),log(500)-log(225),log(500)-log(250),
        log(100)-log(50),log(3)-log(1.5)};
    for(i = 0; i<52; ++i) sigmaprior[i] = pow(sigma_par[i]/2,2);  //same for cancer and non-cancer;
    //    for(i = 52; i<147; ++i) sigmaprior[i] = 0.01*Y_var[i-52]*Y_var[i-52]; // sigma^2
    
    //    //assume the prior of var(epsilon) is inverse Gamma distribution
    //           for(i=0; i<n_Y; ++i){
    //               alpha_tmp[i] = pow(muprior[N_tmp-n_Y+i]/sigmaprior[N_tmp-n_Y+i],2) + 2;
    //               beta_tmp[i] = muprior[N_tmp-n_Y+i] * (alpha_tmp[i] - 1);
    //           }
    
    int indexFit [147];
    for(i=0; i<147; ++i) indexFit[i] = 0;
    for(i=0; i<52; ++i)
    {
        indexFit[i] = 1;
        startp[i] = log(Optimize_orig[i]);
        startp2[i] = log(Optimize_orig[i]);
    }
    
    std::mt19937 rng;
    rng.seed(4);
    std::normal_distribution<double> distribution0(0,1);
    
    for(i = 0; i<3; ++i){
        for(j = 0; j<95; ++j)
        {
            mmYobs[i][j] = Yode[j] + distribution0(rng) * pow(Y_var[j],0.5);
            mmYobs2[i][j] = Yode2[j] + distribution0(rng) * pow(Y_var2[j],0.5);
        }
    }  //mmYobs = [Yode;Yode;Yode] + randn(3,length(Yode))*R; % with 3 obs;
    
    for(i = 0; i <147; ++i)
    if(indexFit[i] != 0) ifFit.push_back(i);
    
    for(i=0; i<ifFit.size();++i)
    {
        sd_tune.push_back(2.38*2.38/ifFit.size());
        sd_tune2.push_back(2.38*2.38/ifFit.size());
    }
    
    
    genDRAM_GibbsWithMH();
    
    return 0;
} // end of main function;



