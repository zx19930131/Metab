//
//  genDRAM_GibbsWithMH.cpp
//
//
//  Created by X.Z on 10/4/18.
//

#include "genDRAM_GibbsWithMH.h"
#include <iomanip>
// 'startp', give some vector as starting point,including sigma^2 and constat parameters;
// 'muprior', 'sigmaprior', assume \theta_i ~ MVN(mu_i, sigma_i^2) and mutually indep;
//
// 'indexFit', indicator for this parameter to be fitted or not with 1 or 0
//
// only have observations at 24h, but with multiple independent replicates for this time
// point
//
// mmYobs: a matrix, one row for one observation at 24h
//
//

extern const int cov_update_ind (30);
extern const int adaptint (5);
extern const double qcovadj (1e-6);

Eigen::MatrixXd ps_hist_total, ps_hist_total2;

Eigen::MatrixXd chaincov, chaincov2;//= Eigen::MatrixXd::Zero(ifFit.size(),ifFit.size());
double delta12;

void genDRAM_GibbsWithMH()
{
    //    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    
    //double ps [nruns] [N_tmp];
    int i,j; // for loop
    extern const int cov_update_ind;
    extern const int adaptint;
    extern const double qcovadj;
    
    double lb_old_tmp [] = {0,0.1,0,0.1,0,0,0,0,0,50,
        0,0,0,0,0,0,0,0,0,100,
        0.5,10,0,0,0.1,0,0.2e6,100,0.1,0.2,
        0.02,0.2,1e4,10,2e4,2e4,0,0,0,0,
        1,0.05,0.1,0.01,0.5,0.1,5,1,50,1,
        1,0};
    double ub_old_tmp [] = {0.01,1,0.1,1,0.01,5e-4,1,2e-4,1e-1,500,
        10,0.01,0.1,0.001,0.2,0.02,0.03,0.01,0.05,1000,
        20,100,2,3,10,3,1e7,1000,4,10,
        1,10,1e5,100,2e5,2e5,0.25,0.01,0.5,0.3,
        50,5,8,1,10,10,100,50,500,500,
        100,3};
    
    double p_curr [ifFit.size()], p_curr2 [ifFit.size()];
    double lb[ifFit.size()];
    double ub[ifFit.size()];
    for(i = 0; i < ifFit.size(); ++i)
    {
        p_curr[i] = startp[ifFit[i]];  // exlude the non-fitted parameters
        p_curr2[i] = startp2[ifFit[i]];
        lb[i] = lb_old_tmp[ifFit[i]];
        ub[i] = ub_old_tmp[ifFit[i]];
    }
    
    //reparameterize beta11 beta12 into beta11 and delta1;
    delta12 = startp[5] - startp2[5];
    
    double lpost = 0.0, lpost2 = 0.0;
    lpostLik(lpost,p_curr,1);
    double L_curr = -2*lpost;
    double L_trial;
    double L_trial2;
    
    lpostLik(lpost2, p_curr2,2);
    double L_curr2 = -2*lpost2;
    //    double L_trial_non;
    //    double L_trial2_non;
    
    double sp_trial;
    double qmu_c_curr, qmu_c_curr2;
    double qcov_c_curr, qcov_c_curr2;
    double tp_trial [ifFit.size()];
    bool bound_flag;
    double qcov2;
    double new_sp_trial;
    
    const double min_accept = 0.4;  // 0.4    0.7
    const double max_accept = 0.7;
    
    ps_hist_total = Eigen::MatrixXd::Zero(nruns*nthinning+nburnin,ifFit.size());
    ps_hist_total2 = Eigen::MatrixXd::Zero(nruns*nthinning+nburnin,ifFit.size());
    const int Cmax = 1e8;
    const double Cmin = 1e-8;
    const double Cmod = 1.1;
    //Cfactor = 0.02;
    double Cfactor [ifFit.size()], Cfactor2 [ifFit.size()];
    for(i=0; i<ifFit.size(); ++i)
    {
        Cfactor[i] = 2.38*2.38/ifFit.size()/ifFit.size();
        Cfactor2[i] = 2.38*2.38/ifFit.size()/ifFit.size();
        
    }
    
    // mcmc;
    bool dr_flag = true;
    const int drscale = 100;
    int save_q;
    double save_r;
    const int file_save = 500;// for result data saving;
    
    std::mt19937 generator;
    generator.seed(44);
    std::normal_distribution<double> distribution(0.0,1.0);
    std::uniform_real_distribution<double> distribution2(0.0,1.0);
    double randa;
    double mymin;
    double a12 = 0;
    bool qa = false;
    double a32, l2, q1, a13;
    
    
    std::cout << "MCMC sampling..." << endl;
    
    //int jthin = 1;
    int jcount = 1;
    int jrungo = -nburnin+1;
    
    // covariance update uses these to store previous values
    Eigen::MatrixXd CAdapt (ifFit.size(),ifFit.size()), CAdapt2 (ifFit.size(),ifFit.size());
    double cov11;
    Eigen::MatrixXd cov12, cov22, cov_tmp;
    chaincov = Eigen::MatrixXd::Zero(ifFit.size(),ifFit.size());
    chaincov2 = Eigen::MatrixXd::Zero(ifFit.size(),ifFit.size());
    
    const int naccepts = 30;
    int accepts [ifFit.size()][naccepts], accepts2 [ifFit.size()][naccepts]; //accept or not, only count for the most recent 30 ones
    int i_accepts [ifFit.size()];
    for(i=0; i<ifFit.size(); ++i) i_accepts[i] = 1; //= repmat(1,sum(ifFit),1);
    double accept_rate, accept_rate2; //acceptance monitor;
    
    for(int jruns = 0; jruns<(nruns*nthinning)+nburnin; ++jruns)
    {
        //cout << "j" << jruns << endl;
        
        methodGibbs(jruns); //!!
        
        //        double accept_rate, accept_rate2; //acceptance monitor;
        
        for(int jGibb = 0; jGibb<ifFit.size(); ++jGibb)
        {
            qmu_c_curr = p_curr[jGibb];
            if(ifFit[jGibb] == 5){qmu_c_curr2 = delta12;}
            else {qmu_c_curr2 = p_curr2[jGibb];}
            if(jruns<adaptint*nthinning+nburnin + cov_update_ind)
            {
                qcov_c_curr = Cfactor[jGibb];
                qcov_c_curr2 = Cfactor2[jGibb];
            }
            else{
                if(ifFit.size() != 1)
                {
                    reloCov(jGibb, CAdapt, CAdapt2);//!!
                    
                    cov11 = CAdapt(CAdapt.rows()-1,CAdapt.rows()-1);
                    cov12 = CAdapt.bottomLeftCorner(1,CAdapt.rows()-1);
                    cov22 = CAdapt.topLeftCorner(CAdapt.rows()-1,CAdapt.rows()-1);//(1:(end-1),1:(end-1));
                    cov_tmp = cov12 * (cov22.colPivHouseholderQr().solve(cov12.transpose()));
                    qcov_c_curr = sd_tune[jGibb]*(cov11 - cov_tmp(0,0));
                    
                    cov11 = CAdapt2(CAdapt2.rows()-1,CAdapt2.rows()-1);
                    cov12 = CAdapt2.bottomLeftCorner(1,CAdapt2.rows()-1);
                    cov22 = CAdapt2.topLeftCorner(CAdapt2.rows()-1,CAdapt2.rows()-1);//(1:(end-1),1:(end-1));
                    cov_tmp = cov12 * (cov22.colPivHouseholderQr().solve(cov12.transpose()));
                    qcov_c_curr2 = sd_tune2[jGibb]*(cov11 - cov_tmp(0,0));
                }
                else if(ifFit.size() == 1) {
                    qcov_c_curr = sd_tune[jGibb]*chaincov(0,0);
                    qcov_c_curr2 = sd_tune2[jGibb]*chaincov2(0,0);
                }
            }
            
            
            
            //###############################Cancer#########################################//
            sp_trial = qmu_c_curr + distribution(generator)*pow(qcov_c_curr,0.5);
            
            if(exp(sp_trial)<lb[jGibb] && ifFit[jGibb]<N_tmp-n_Y ) sp_trial = log(exp(sp_trial) + 2*(lb[jGibb] - exp(sp_trial)));
            else if(exp(sp_trial)>ub[jGibb] && ifFit[jGibb]<N_tmp-n_Y) sp_trial = log(exp(sp_trial) + 2*(ub[jGibb] - exp(sp_trial)));
            
            a12 = 0;
            qa = false;
            
            bound_flag = (exp(sp_trial)>=lb[jGibb])? true:false;
            bound_flag = (exp(sp_trial)<=ub[jGibb])? true:false;
            
            for(i=0; i<ifFit.size(); ++i) tp_trial[i] = p_curr[i];
            
            if(bound_flag)
            {
                tp_trial[jGibb] = sp_trial;
                lpostLik(lpost, tp_trial,1);
                L_trial = -2*lpost;
                
                a12 = exp(-0.5*(L_trial - L_curr));
                randa = distribution2(generator);
                mymin = (a12<=1)?a12:1;
                qa = randa <= mymin; //% accept?
            }
            
            // update
            if(qa)   // if accept this trial
            {
                p_curr[jGibb] = sp_trial;
                L_curr = L_trial;
            }
            else if(!qa && dr_flag && jrungo>4000) //do delayed rejection
            {
                qcov2 = qcov_c_curr/drscale;
                new_sp_trial = qmu_c_curr + distribution(generator)*pow(qcov2,0.5);
                
                if(exp(new_sp_trial)<lb[jGibb] && ifFit[jGibb]<N_tmp-n_Y ) new_sp_trial = log(exp(new_sp_trial) + 2*(lb[jGibb] - exp(new_sp_trial)));
                else if(exp(new_sp_trial)>ub[jGibb] && ifFit[jGibb]<N_tmp-n_Y) new_sp_trial = log(exp(new_sp_trial) + 2*(ub[jGibb] - exp(new_sp_trial)));
                
                bound_flag = (exp(new_sp_trial)>=lb[jGibb])? true:false;
                bound_flag = (exp(new_sp_trial)<=ub[jGibb])? true:false;
                
                if(bound_flag)
                {
                    tp_trial[jGibb] = new_sp_trial;
                    lpostLik(lpost, tp_trial,1);
                    L_trial2 = -2*lpost;
                    a32 = (exp(-0.5*(L_trial - L_trial2))<=1)? exp(-0.5*(L_trial - L_trial2)):1;
                    l2 = exp(-0.5*(L_trial2-L_curr));
                    q1 = exp(-0.5*(pow(sp_trial - new_sp_trial,2)/qcov_c_curr + pow(sp_trial - qmu_c_curr,2)/qcov_c_curr));
                    a13 = l2*q1*(1-a32)/(1-((a12<=1)?a12:1));
                    randa = distribution2(generator);
                    mymin = (a13<=1)?a13:1;
                    qa = randa <= mymin; //% accept?
                }
                if(qa)   // if accept this trial
                {
                    p_curr[jGibb] = new_sp_trial;
                    L_curr = L_trial2;
                }
            }
            //###############################Cancer#########################################//
            accepts[jGibb][i_accepts[jGibb]-1] = qa? 1:0;
            accept_rate = std::accumulate(&accepts[jGibb][0],&accepts[jGibb][naccepts-1],0.0)/(naccepts-1);
            
            
            //###############################Non-Cancer#####################################//
            if(ifFit[jGibb] == 5){
                sp_trial = p_curr[jGibb] - qmu_c_curr2 - distribution(generator)*pow(qcov_c_curr2,0.5);//b12 = b11 - delta12;
            }else{sp_trial = qmu_c_curr2 + distribution(generator)*pow(qcov_c_curr2,0.5);}
            
            if(exp(sp_trial)<lb[jGibb] && ifFit[jGibb]<N_tmp-n_Y) sp_trial = log(exp(sp_trial) + 2*(lb[jGibb] - exp(sp_trial)));
            else if(exp(sp_trial)>ub[jGibb] && ifFit[jGibb]<N_tmp-n_Y) sp_trial = log(exp(sp_trial) + 2*(ub[jGibb] - exp(sp_trial)));
            
            a12 = 0;
            qa = false;
            
            bound_flag = (exp(sp_trial)>=lb[jGibb])? true:false;
            bound_flag = (exp(sp_trial)<=ub[jGibb])? true:false;
            
            for(i=0; i<ifFit.size(); ++i) tp_trial[i] = p_curr2[i];
            
            if(bound_flag)
            {
                tp_trial[jGibb] = sp_trial;
                lpostLik(lpost2, tp_trial,2);
                L_trial = -2*lpost2;
                
                a12 = exp(-0.5*(L_trial - L_curr2));
                randa = distribution2(generator);
                mymin = (a12<=1)?a12:1;
                qa = randa <= mymin; //% accept?
            }
            
            // update
            if(qa)   // if accept this trial
            {
                p_curr2[jGibb] = sp_trial;
                L_curr2 = L_trial;
            }
            else if(!qa && dr_flag && jrungo>4000) //do delayed rejection
            {
                qcov2 = qcov_c_curr2/drscale;
                //                new_sp_trial = qmu_c_curr2 + distribution(generator)*pow(qcov2,0.5);
                if(ifFit[jGibb] == 5){
                    new_sp_trial = p_curr[jGibb] - qmu_c_curr2 - distribution(generator)*pow(qcov2,0.5);//b12 = b11 - delta12;
                }else{new_sp_trial = qmu_c_curr2 + distribution(generator)*pow(qcov2,0.5);}
                
                
                if(exp(new_sp_trial)<lb[jGibb] && ifFit[jGibb]<N_tmp-n_Y ) new_sp_trial = log(exp(new_sp_trial) + 2*(lb[jGibb] - exp(new_sp_trial)));
                else if(exp(new_sp_trial)>ub[jGibb] && ifFit[jGibb]<N_tmp-n_Y) new_sp_trial = log(exp(new_sp_trial) + 2*(ub[jGibb] - exp(new_sp_trial)));
                
                bound_flag = (exp(new_sp_trial)>=lb[jGibb])? true:false;
                bound_flag = (exp(new_sp_trial)<=ub[jGibb])? true:false;
                
                if(bound_flag)
                {
                    tp_trial[jGibb] = new_sp_trial;
                    lpostLik(lpost2, tp_trial,2);
                    L_trial2 = -2*lpost2;
                    a32 = (exp(-0.5*(L_trial - L_trial2))<=1)? exp(-0.5*(L_trial - L_trial2)):1;
                    l2 = exp(-0.5*(L_trial2-L_curr2));
                    q1 = exp(-0.5*(pow(sp_trial - new_sp_trial,2)/qcov_c_curr2 + pow(sp_trial - qmu_c_curr2,2)/qcov_c_curr2));
                    a13 = l2*q1*(1-a32)/(1-((a12<=1)?a12:1));
                    randa = distribution2(generator);
                    mymin = (a13<=1)?a13:1;
                    qa = randa <= mymin; //% accept?
                }
                if(qa)   // if accept this trial
                {
                    p_curr2[jGibb] = new_sp_trial;
                    L_curr2 = L_trial2;
                }
            }
            
            if(ifFit[jGibb] == 5) {delta12 = p_curr[jGibb] - p_curr2[jGibb];}
            //###############################Non-Cancer#####################################//
            accepts2[jGibb][i_accepts[jGibb]-1] = qa? 1:0;
            accept_rate2 = std::accumulate(&accepts2[jGibb][0],&accepts2[jGibb][naccepts-1],0.0)/(naccepts-1);
            
            
            
            
            i_accepts[jGibb] = i_accepts[jGibb] + 1;
            if(i_accepts[jGibb]>naccepts) i_accepts[jGibb] = 1;    //only look at the accept rate of previous 30 ones
            
            // adjust scaling during burn-in
            if(jrungo<=0){
                if(accept_rate > max_accept && Cfactor[jGibb]*Cmod<Cmax)
                Cfactor[jGibb] = Cfactor[jGibb]*Cmod;
                else if(accept_rate < min_accept && Cfactor[jGibb]/Cmod>Cmin)
                Cfactor[jGibb] = Cfactor[jGibb]/Cmod;
                
                if(accept_rate2 > max_accept && Cfactor2[jGibb]*Cmod<Cmax)
                Cfactor2[jGibb] = Cfactor2[jGibb]*Cmod;
                else if(accept_rate2 < min_accept && Cfactor2[jGibb]/Cmod>Cmin)
                Cfactor2[jGibb] = Cfactor2[jGibb]/Cmod;
            }
            else{ // adjust scaling after burn-in
                if(accept_rate > max_accept && sd_tune[jGibb]*Cmod<Cmax)
                sd_tune[jGibb] = sd_tune[jGibb]*Cmod;
                else if(accept_rate < min_accept && sd_tune[jGibb]/Cmod>Cmin)
                sd_tune[jGibb] = sd_tune[jGibb]/Cmod;
                
                if(accept_rate2 > max_accept && sd_tune2[jGibb]*Cmod<Cmax)
                sd_tune2[jGibb] = sd_tune2[jGibb]*Cmod;
                else if(accept_rate2 < min_accept && sd_tune2[jGibb]/Cmod>Cmin)
                sd_tune2[jGibb] = sd_tune2[jGibb]/Cmod;
            }
            
        }// end for jGibb
        
        
        for(i=0; i<ifFit.size(); ++i)
        {
            ps_hist_total(jruns,i) = p_curr[i];
        }
        
        for(i=0; i<ifFit.size(); ++i)
        {
            ps_hist_total2(jruns,i) = p_curr2[i];
        }
        ps_hist_total2(jruns,5) = delta12;////////////////////////////////////////////////
        
        
        //save samples
        if(jrungo>0)
        {
            jcount = jcount + 1;
            //}
        }
        
        std::div_t divre = div(jruns+1,file_save);
        save_q = divre.quot;
        save_r = divre.rem;
        if(save_r == 0 & save_q > 0)
        {
            string file_name = "mysample_curr_" + std::to_string(save_q) + ".txt";
            ofstream myfile;
            myfile.open (file_name);
            for(i=0; i<ifFit.size(); ++i)
            {
                myfile << "par" << ifFit[i]+1 << '\t';
            }
            myfile << '\n';
            
            for(i=0; i<file_save; ++i)
            {
                for(j=0; j<ifFit.size(); ++j) myfile << std::scientific << std::setprecision(8) << exp(ps_hist_total(file_save*(save_q-1)+i,j)) << '\t';// exp(ps[file_save*(save_q-1)+i][ifFit[j]]);
                myfile << '\n';
            }
            
            myfile.close();
            
            file_name = "mysample2_curr_" + std::to_string(save_q) + ".txt";
            ofstream myfile2;
            myfile2.open (file_name);
            for(i=0; i<ifFit.size(); ++i)
            {
                myfile2 << "par" << ifFit[i]+1 << '\t';
            }
            myfile2 << '\n';
            
            for(i=0; i<file_save; ++i)
            {
                for(j=0; j<ifFit.size(); ++j)
                {
                    if(ifFit[j]==5) {
                        myfile2 << std::scientific << std::setprecision(8) << ps_hist_total2(file_save*(save_q-1)+i,j) << '\t';
                    }
                    else {
                        myfile2 << std::scientific << std::setprecision(8) << exp(ps_hist_total2(file_save*(save_q-1)+i,j)) << '\t';
                    }
                }
                myfile2 << '\n';
            }
            
            myfile2.close();
        }
        
        jrungo = jrungo + 1;
    }
    
    cout << "MCMC sampling done. " << '\n';
    return;
}

void lpostLik(double& lpost, double *ptmp, int ind)
{
    //extern int observed [95];
    double inputp_tmp[N_tmp];
    int i,j; //for loop
    state_type Xd_tmp;
    extern double odeParameters[57];
    double Yode_tmp[n_Y];
    double llikp_tmp (0);
    
    if(ind == 1)
    {
        for(i=0; i<N_tmp; ++i) inputp_tmp[i] = startp[i];
        for(i = 0; i < ifFit.size(); ++i)inputp_tmp[ifFit[i]] = ptmp[i];
        
        extern state_type Xd_flag;
        for(i=0; i<512; ++i) Xd_tmp[i] = Xd_flag[i];
        for(i=0; i<(N_tmp-n_Y); i++) odeParameters[i] = exp(inputp_tmp[i]);
        integrate(PurineSynthesis, Xd_tmp, t0, tf, df);
        
        //exclude unobserved data with 0
        for(i=0; i<n_Y; ++i) Yode_tmp[i]= log(Xd_tmp[observed[i]]/MoleculeNumberInOneNanoMole);// cout<<Yode_tmp[i]<<'\t';}
        
        llikp_tmp = 0;
        
        for(i=0; i<3; ++i)
        {
            for(j=0; j<n_Y; ++j){
                llikp_tmp += - 0.5*pow(mmYobs[i][j]-Yode_tmp[j],2)/inputp_tmp[N_tmp-n_Y+j];//- 0.5*log(2.0 * M_PI)
            }
        }
        
        for(i=0; i<ifFit.size();++i)
        {
            llikp_tmp +=  - 0.5*pow(inputp_tmp[ifFit[i]]-muprior[ifFit[i]],2)/sigmaprior[ifFit[i]];
        }
    }  //################################next for non-cancer##################################;
    else if(ind == 2)
    {
        for(i=0; i<N_tmp; ++i) inputp_tmp[i] = startp2[i];
        for(i = 0; i < ifFit.size(); ++i)inputp_tmp[ifFit[i]] = ptmp[i];
        
        extern state_type Xd_flag2;
        for(i=0; i<512; ++i) Xd_tmp[i] = Xd_flag2[i];
        for(i=0; i<(N_tmp-n_Y); i++) odeParameters[i] = exp(inputp_tmp[i]);
        integrate(PurineSynthesis, Xd_tmp, t0, tf, df);
        
        //exclude unobserved data with 0
        for(i=0; i<n_Y; ++i) Yode_tmp[i]= log(Xd_tmp[observed[i]]/MoleculeNumberInOneNanoMole);// cout<<Yode_tmp[i]<<'\t';}
        
        llikp_tmp = 0;
        
        for(i=0; i<3; ++i)
        {
            for(j=0; j<n_Y; ++j){
                llikp_tmp += - 0.5*pow(mmYobs2[i][j]-Yode_tmp[j],2)/inputp_tmp[N_tmp-n_Y+j];//- 0.5*log(2.0 * M_PI)
            }
        }
        
        llikp_tmp +=  - 0.5*pow(delta12-muprior[5]+muprior2[5],2)/sigmaprior[5];
        for(i=0; i<ifFit.size();++i)
        {
            if(ifFit[i]!=5)llikp_tmp +=  - 0.5*pow(inputp_tmp[ifFit[i]]-muprior2[ifFit[i]],2)/sigmaprior[ifFit[i]];
        }
    }
    
    lpost = llikp_tmp;
}


void methodGibbs(int& jruns)
{
    std::div_t divre = div(jruns,cov_update_ind);
    double cov_r = divre.rem;
    
    if((jruns == adaptint*nthinning+nburnin + cov_update_ind) || (cov_r == 0 && jruns > adaptint*nthinning+nburnin + cov_update_ind))
    {
        Eigen::MatrixXd ps_hist_tmp(jruns-nburnin,ifFit.size());
        for(int i=0; i<(jruns-nburnin); ++i)
        for(int j=0; j< ifFit.size(); ++j) ps_hist_tmp(i,j) = ps_hist_total(nburnin+i,j);
        Eigen::MatrixXd centered = ps_hist_tmp.rowwise() - ps_hist_tmp.colwise().mean();
        Eigen::MatrixXd cov_data = (centered.adjoint() * centered) / double(ps_hist_tmp.rows() - 1);
        //.adjoint() = .transpose() here;
        chaincov = cov_data + Eigen::MatrixXd::Identity(ifFit.size(), ifFit.size())*qcovadj;
        
        for(int i=0; i<(jruns-nburnin); ++i)
        for(int j=0; j< ifFit.size(); ++j) ps_hist_tmp(i,j) = ps_hist_total2(nburnin+i,j);
        centered = ps_hist_tmp.rowwise() - ps_hist_tmp.colwise().mean();
        cov_data = (centered.adjoint() * centered) / double(ps_hist_tmp.rows() - 1);
        chaincov2 = cov_data + Eigen::MatrixXd::Identity(ifFit.size(), ifFit.size())*qcovadj;
    }
}

template <typename Derived>
void reloCov(int index, Eigen::MatrixBase<Derived>& CAdapt, Eigen::MatrixBase<Derived>& CAdapt2)
{
    //Eigen::MatrixXd reloM;
    int i,j;
    
    if(ifFit.size()>1 && index != ifFit.size() - 1)
    {
        Eigen::MatrixXd tmpM (chaincov), tmpM2 (chaincov), tmpM_n (chaincov2), tmpM2_n (chaincov2);
        
        //change the rows;
        for(i=0; i<ifFit.size(); ++i)
        for(j=0; j<ifFit.size(); ++j)
        {
            if(i >= index && i <ifFit.size()-1)
            {
                tmpM(i,j) = chaincov(i+1,j);
                tmpM_n(i,j) = chaincov2(i+1,j);
            }
            else if(i == ifFit.size()-1)
            {
                tmpM(i,j) = chaincov(index,j);
                tmpM_n(i,j) = chaincov2(index,j);
            }
        }
        tmpM2 = tmpM;
        tmpM2_n = tmpM_n;
        //change the columns;
        for(j=0; j<ifFit.size(); ++j)
        for(i=0; i<ifFit.size(); ++i)
        {
            if(j >= index && j <ifFit.size()-1)
            {
                tmpM(i,j) = tmpM2(i,j+1); tmpM_n(i,j) = tmpM2_n(i,j+1);
            }
            else if(j == ifFit.size()-1)
            {
                tmpM(i,j) = tmpM2(i,index);tmpM_n(i,j) = tmpM2_n(i,index);
            }
        }
        CAdapt = tmpM;
        CAdapt2 = tmpM_n;
    }
    else if(ifFit.size() == 1 | index == ifFit.size() - 1) {CAdapt = chaincov; CAdapt2 = chaincov2;}
    
}



























