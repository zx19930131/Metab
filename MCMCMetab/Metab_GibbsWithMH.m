% 'startp', give some vector as starting point,including sigma^2 and constat parameters;
% 'muprior', 'sigmaprior', assume \theta_i ~ MVN(mu_i, sigma_i^2) and mutually indep;
%
% 'indexFit', indicator for this parameter to be fitted or not with 1 or 0
%
%
%
%



function [ ps, ps_trial, chi2s, chi2s_trial, acceptance ] = Metab_GibbsWithMH(Yobs, Yode, startp, muprior,sigmaprior, indexFit, nruns, nburnin, nthinning,Xd,observed,MoleculeNumberInOneNanoMole,ode_fun)

%initial values of the chain
ifFit = indexFit==1;
pReset = startp;
p_curr = startp(ifFit); % exlude the non-fitted parameters
mupriorFit = muprior(ifFit);
sigmapriorFit = sigmaprior(ifFit);

ps = nan(nruns, length(startp)); % save the samples
ps_trial = nan(nruns, length(startp));
chi2s = nan(1,nruns); % -2log(post)
chi2s_trial = nan(1,nruns);
acceptance = nan(1,nruns);
nwindow = sum(ifFit)*50;

%p_trial = p_curr;
pp_trial = p_curr; % nan(1,length(p_curr); % for within Gibb usage
                         
N_tmp = length(startp);
n_Y = length(Yode);
logpostlikeli = @lpostLik;
lpost_curr = feval(logpostlikeli,p_curr);
L_curr = -2*lpost_curr;
                         
min_accept = 0.4;
max_accept = 0.7;
                         
lb = [0,0.1,0,0.1,0,0,0.01,0,0.1e-3,50,1,0,0,0,0,0,0,0,0,50,1,10,0,0,0.1,...
      0,0.2e6,100,0.1,0.2,0.02,0.2,1e4,10,5e4,5e4,0,0,0,0,0.5,0.05,0.1,0.01,0.5,0.1,...
      5,1,50,0.5,1,0,repmat(0,1,n_Y)];
% lb = [30,0.1,0.1,0.1,0.1,0.002] ;
lb = lb(ifFit);
ub = [0.01,1,0.1,1,0.01,1e-4,0.1,1e-4,5e-3,500,10,0.01,0.1,0.001,0.01,0.01,0.03,0.01,0.05,500,10,50,2,3,10,...
      3,1e7,1000,4,10,1,10,1e5,100,5e5,5e5,0.25,0.01,0.1,0.3,5,5,8,1,10,10,...
      100,50,500,10,50,0.2,repmat(1e20,1,n_Y)];
% ub = [100,5,5,5,5,0.4];
ub = ub(ifFit);
                         
method = @condi_Gibbs; % @mcmc_adaptive;
%adjust_scaling = true;
ps_hist = nan(nwindow, sum(ifFit)); % a matrix to store the history of the parameter values
ps_hist_index = 1;
CAdapt = eye(sum(ifFit));
                         
ps_hist_total = nan(nruns*nthinning+nburnin,sum(ifFit));
                         
Cmax = 1e8;
Cmin = 1e-8;
Cmod = 1.1;
Cfactor = 0.02;
%Cfactor = (2.38/sqrt(sum(indexFit)))^2 / sum(indexFit);
%cScale = 1.0;
                         
%res_curr = [];
%Sres_curr = [];
                         
jrungo = -nburnin+1;
                         
                         % additional functionality
                         %do_chain_resets = false;
                         %do_reflect_bounds = true;
                         
% mcmc
%arWaitbar(0);
naccepts = 30;
accepts = nan(1,naccepts); %accept or not, only count for the most recent 30 ones
i_accepts = 1;
fprintf('MCMC sampling...')
%tic;
%count_chain_reset = 0;
jthin = 1;
jcount = 1;
for jruns = 1:((nruns*nthinning)+nburnin)
                         qnonnanacc = ~isnan(accepts);
                         accept_rate = sum(accepts(qnonnanacc))/length(accepts(qnonnanacc));
%                          if(jrungo>0) %after burn-in
%                          arWaitbar(jruns, (nruns*nthinning)+nburnin, sprintf('MCMC run (acceptance rate %4.1f%%)', ...
%                                                                              accept_rate*100));
%                          else
%                          arWaitbar(jruns, (nruns*nthinning)+nburnin, sprintf('MCMC burn-in (acceptance rate %4.1f%%)', ...
%                                                                              accept_rate*100));
%                          end
                         
                         condi_accept = 0;
      for jGibb = 1:sum(ifFit)
                         
                         %[mu_curr, covar_curr] = feval(method, p_curr);
                         [qmu_c_curr, qcov_c_curr] = feval(method,p_curr);
                         sp_trial = qmu_c_curr + randn * qcov_c_curr; % one-dimension normal here
                         
                         L_trial = 0;
                         
                         % reflect from bejond bounds
                         if(sum(sp_trial<lb(jGibb)) + sum(sp_trial>ub(jGibb)) > 0)
                         qlb = sp_trial<lb(jGibb);
                         sp_trial(qlb) = sp_trial(qlb) + 2*(lb(jGibb) - sp_trial(qlb));
                         qub = sp_trial>ub(jGibb);
                         sp_trial(qub) = sp_trial(qub) + 2*(ub(jGibb) - sp_trial(qub));
                         end
                         
                         pp_trial(jGibb) = sp_trial;
                         
                       if(sum(sp_trial<lb(jGibb)) + sum(sp_trial>ub(jGibb)) == 0) % check bounds
                         try
                         tp_trial = p_curr;
                         tp_trial(jGibb) = sp_trial;
                         lpost_trial = feval(logpostlikeli,tp_trial);
                         L_trial = -2*lpost_trial;
                         [qmu_c_trial, qcov_c_trial] = feval(method, tp_trial);
                         Q_trial = normpdf(sp_trial, qmu_c_curr, qcov_c_curr);
                         Q_curr =normpdf(p_curr(jGibb), qmu_c_trial, qcov_c_trial);
                         a = exp(-0.5*(L_trial - L_curr)) * (Q_curr / Q_trial);
                         randa = rand;
                         qa = randa <= min([1 a]); % accept?
                         catch err_id
                         fprintf('#%i error: %s\n', jrungo, err_id.message);
                         qa = false;
                         end
                       else
                         qa = false;
                       end
                         
% adjust scaling during burn-in
if(jrungo<=0) % && adjust_scaling)
    if(accept_rate > max_accept && Cfactor*Cmod<Cmax)
        Cfactor = Cfactor*Cmod;
    elseif(accept_rate < min_accept && Cfactor/Cmod>Cmin)
        Cfactor = Cfactor/Cmod;
    end
end

% accepts(i_accepts) = qa;
% i_accepts = i_accepts + 1;
%if(i_accepts>length(accepts))   %only look at the accept rate of previous 30 ones
%   i_accepts = 1;
%end


% update
if(qa)  % if accept this trial
condi_accept = condi_accept + 1;
p_curr(jGibb) = sp_trial;
L_curr = L_trial;
end

end % end for jGibb

accepts(i_accepts) = condi_accept>0;
i_accepts = i_accepts + 1;
if(i_accepts>length(accepts))   %only look at the accept rate of previous 30 ones
    i_accepts = 1;
end


if(condi_accept > 0)
%p_curr = p_trial;
ps_hist(ps_hist_index,:) = p_curr;
ps_hist_index = ps_hist_index + 1;
   if(ps_hist_index > nwindow)
      ps_hist_index = 1;
   end
end      % HAVE QUESTION HERE


ps_hist_total(jruns,:)=p_curr;


% save samples
if(jrungo>0)
    jthin = jthin + 1;
    if(jthin > nthinning)
        ps(jcount,:) = pReset;
        ps(jcount,ifFit) = p_curr;
        ps_trial(jcount,:) = pReset;
        ps_trial(jcount,ifFit) = pp_trial; %if any elemetnts in pp_trial and p_curr is not equal, then the previous step rejected.
        chi2s(jcount) = L_curr;
        chi2s_trial(jcount) = L_trial; % sort of meaningless now
        acceptance(jcount) = accept_rate;
        jthin = 1;
        jcount = jcount + 1;
   end
end


jrungo = jrungo + 1;
end

%arWaitbar(-1);
fprintf('done \n');%, secToHMS(toc));
%mcmc_toc = toc;
%ar.p = pReset;


function [lpost] = lpostLik(ptmp)

inputp_tmp = startp;
inputp_tmp(ifFit) = ptmp;
%N_tmp = length(inputp_tmp);
%n_tmp = n_Y;
%%%%%%%%%%%%%
Parameters_tmp = [inputp_tmp(1:(N_tmp-n_Y)),0.5E-6,7.8E3,26,26*0.64,26*0.16];% OptimizeParameters((end-4):end)];
[t, Ydd_tmp] = ode15s(ode_fun, [0 : 60 : 24*60], Xd, [], Parameters_tmp);
%[t, Ydd_tmp] = ode45(ode_fun, [0 : 60 : 24*60], Xd, [], Parameters_tmp);
Yode_tmp = Ydd_tmp(end,observed)./ MoleculeNumberInOneNanoMole;
%%%%%%%%%%%%%
llikp_tmp = log(mvnpdf(Yobs, Yode_tmp, diag(inputp_tmp((end-n_Y+1):end)))); %put the var(y) in the last

%assume the prior of var(epsilon) is Gamma distribution
beta_tmp = sigmaprior((end-n_Y+1):end)./muprior((end-n_Y+1):end);
alpha_tmp = muprior((end-n_Y+1):end)./beta_tmp;

lprior_tmp = log(mvnpdf(inputp_tmp(1:(N_tmp-n_Y)),muprior(1:(N_tmp-n_Y)),diag(sigmaprior(1:(N_tmp-n_Y)))))...
+log(prod(gampdf(inputp_tmp((end-n_Y+1):end),alpha_tmp,beta_tmp)));
%prior of the parameters

lpost = llikp_tmp + lprior_tmp;
end

%function [mu, covar] = mcmc_adaptive(ptmp, ~, ~)
%        qnonnan = ~isnan(ps_hist(:,1));
%        if(jrungo<=0 | sum(qnonnan)<=1)
%            mu = ptmp;
%            covar = eye(length(ptmp)) * Cfactor;
%        elseif(jrungo>0 & sum(qnonnan)>1)
%
%            CAdapt = cov(ps_hist(qnonnan,:));
%           mu = ptmp;
%            covar = CAdapt;
%        end
%end

function [qmu_condi, qcov_condi] = condi_Gibbs(ptmp) % (... | ptmp)
      qnonnan = ~isnan(ps_hist(:,1));
      if(jrungo<=0 | ps_hist_index == 1)
           qmu_condi = ptmp(jGibb);
           qcov_condi = Cfactor;

      elseif(jrungo>0 & ps_hist_index > 1)
             if(jruns<100*nthinning+nburnin)
                 ps_hist_tmp = ps_hist(qnonnan,:);
             else
                 ps_hist_tmp = ps_hist_total(1:(jruns-1),:);
             end
             ps_hist_flag = ps_hist_tmp;

             if(sum(ifFit)==1)
                 qmu_condi = ptmp(jGibb);
                 qcov_condi = cov(ps_hist_tmp);

             elseif(sum(ifFit)>1)
                 ps_hist_tmp(:,1) = ps_hist_flag(:,jGibb);
                 ps_hist_tmp(:,jGibb) = ps_hist_flag(:,1);
%                  mu_tmp = ptmp;
%                  mu_tmp(1) = ptmp(jGibb);
%                  mu_tmp(jGibb) = ptmp(1);
%                  ptmpF = ptmp;
%                  ptmpF(1) = ptmp(jGibb);
%                  ptmpF(jGibb) = ptmp(1);

                 CAdapt = cov(ps_hist_tmp);
                 cov11 = CAdapt(1,1);
                 cov12 = CAdapt(1,2:end);
                 cov22 = CAdapt(2:end,2:end);
                 qmu_condi = ptmp(jGibb);
                 qcov_condi = cov11 - cov12 * inv(cov22) * cov12.';
             end
      end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

