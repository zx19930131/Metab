//
//  genDRAM_GibbsWithMH.h
//
//
//  Created by X.Z on 10/4/18.
//

#ifndef genDRAM_GibbsWithMH_h
#define genDRAM_GibbsWithMH_h

#include <iostream>
#include <bits/stdc++.h>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
//#include <Eigen/core>
#include <random>
#include <numeric>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/inverse_gamma.hpp>
#include <cstdlib>
#include <fstream>
#include <string>
#include <chrono>
#include "PurineSynthesis.h"
//#include "Constants.h"
//#include "some_extern.h"


extern vector<int> ifFit;
extern vector<double> sd_tune, sd_tune2;
extern    double mmYobs [3][95], mmYobs2 [3][95];
extern    double Yode [95],Yode2[95] ;
extern    double startp [147], startp2 [147] ;
extern    double muprior [147], muprior2 [147];
extern    double sigmaprior [52];
//extern double beta_tmp[95], alpha_tmp[95];
extern    int observed [95];

extern const double t0;
extern const double tf;
extern const double df;

extern const int nruns;
extern  const  int nburnin;
extern const int nthinning;

extern const int N_tmp;
extern const int n_Y;

void lpostLik(double& lpost, double *ptmp, int ind);

void methodGibbs(int& jruns);

template <typename Derived>
void reloCov(int index, Eigen::MatrixBase<Derived>& CAdapt, Eigen::MatrixBase<Derived>& CAdapt2);

void genDRAM_GibbsWithMH();


#endif /* genDRAM_GibbsWithMH_h */
