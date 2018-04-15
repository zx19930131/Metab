function [ ps, ps_trial, chi2s, chi2s_trial, acceptance ] = parallelMCMCMetab(temperature,aux_initials, Yobs, Yode, startp, muprior,sigmaprior, indexFit, nruns, nburnin, nthinning)

temps = temperature;
N_chain = length(temps);
startpp = [startp;aux_initials]; % one row, one chain

%initial values of the chain
ifFit = indexFit==1;
pReset = startpp; % matrix
p_curr = startpp(:,ifFit); % exlude the non-fitted parameters

jindexoffset = 0;
ps = nan(nruns, length(startp)); % save the samples only for the first chain
ps_trial = nan(nruns, length(startp));
chi2s = nan(1,nruns); % -2log(post)
chi2s_trial = nan(1,nruns);
acceptance = nan(1,nruns);
nwindow = sum(indexFit == 1)*50;

n_Y = length(Yode);
postlikelihood = @postLik;
post_curr = arrayfun(@(n) postlikelihood(p_curr(n,:)), 1:size(p_curr,1))
%post_curr = feval(postlikelihood,p_curr);
L_curr = -2*log(post_curr);

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

lbm = repmat(lb,N_chain,1);
ubm = repmat(ub,N_chain,1);

method = @mcmc_adaptive;
adjust_scaling = true;
ps_hist = nan(nwindow, sum(indexFit == 1),N_chain); % one matrix, one chain,a matrix to store the history of the parameter values
ps_hist_index = 1;
CAdapt = eye(sum(indexFit == 1));

Cmax = 1e8;
Cmin = 1e-8;
Cmod = 1.1;
Cfactor = (2.38/sqrt(sum(indexFit)))^2 / sum(indexFit);
%cScale = 1.0;

%res_curr = [];
%Sres_curr = [];

jrungo = -nburnin+1;
    
% additional functionality
do_chain_resets = false;
do_reflect_bounds = true;

% mcmc
arWaitbar(0);
naccepts = 30;
accepts = nan(1,naccepts,N_chain); %accept or not, only count for the most recent 30 ones
i_accepts = 1;
fprintf('MCMC sampling...')
tic;
count_chain_reset = 0;
jthin = 1;
jcount = 1;
for jruns = 1:((nruns*nthinning)+nburnin)
    qnonnanacc = ~isnan(accepts(1,:,1));
    accept_rate = sum(accepts(1,qnonnanacc,1))/length(accepts(1,qnonnanacc,1)); % only for the first chain
    if(jrungo>0) %after burn-in
        arWaitbar(jruns, (nruns*nthinning)+nburnin, sprintf('MCMC run (acceptance rate %4.1f%%)', ...
            accept_rate*100));
    else
        arWaitbar(jruns, (nruns*nthinning)+nburnin, sprintf('MCMC burn-in (acceptance rate %4.1f%%)', ...
            accept_rate(1)*100));
    end
   
    p_trial = nan(N_chain,size(p_curr,2));
    for index = 1:N_chain
    [mu_curr, covar_curr] = feval(method, p_curr(index,:),index); 
    p_trial(index,:) = mvnrnd(mu_curr, covar_curr);
    end
    
    L_trial = repmat(0,1,N_chain);
    
    % reflect from bejond bounds
    if(do_reflect_bounds == true)
        if(sum(sum(p_trial<lbm,2) + sum(p_trial>ubm,2)) > 0)
            qlb = p_trial<lbm;
            p_trial(qlb) = p_trial(qlb) + 2*(lbm(qlb) - p_trial(qlb));
            qub = p_trial>ubm;
            p_trial(qub) = p_trial(qub) + 2*(ubm(qub) - p_trial(qub));
        end
%         for index = 1:N_chain
%          if(sum(p_trial(index,:)<lb) + sum(p_trial(index,:)>ub) > 0)
%              qlb = p_trial(index,:)<lb;
%              p_trial(index,qlb) = p_trial(index,qlb) + 2*(lb(qlb) - p_trial(index,qlb));
%              qub = p_trial(index,:)>ub;
%              p_trial(index,qub) = p_trial(index,qub) + 2*(ub(qub) - p_trial(index,qub));
%          end
%         end
    end
 
    
    if(sum(sum(p_trial<lbm,2) + sum(p_trial>ubm,2)) == 0) % check bounds
        % fprintf('#%i bounds ok\n', jrungo);
        try
            %res_trial = [];
            %Sres_trial = [];
            post_trial = feval(postlikelihood,p_trial); 
            L_trial = -2*log(post_trial);
            %L_trial = -2*log(mvnpdf(Yobs, Yode, diag(p_trial((end-n_Y+1):end)))*mvnpdf(p_trial,muprior,diag(sigmaprior)));
            [mu_trial, covar_trial] = feval(method, p_trial);%, res_trial, Sres_trial);
            Q_trial = mvnpdf(p_trial, mu_curr, covar_curr);
            Q_curr = mvnpdf(p_curr, mu_trial, covar_trial);
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
    if(jrungo<=0 && adjust_scaling)
        if(accept_rate > max_accept && Cfactor*Cmod<Cmax)
            Cfactor = Cfactor*Cmod;
        elseif(accept_rate < min_accept && Cfactor/Cmod>Cmin)
            Cfactor = Cfactor/Cmod;
        end
    end
    
    accepts(i_accepts) = qa;
    i_accepts = i_accepts + 1;
    if(i_accepts>length(accepts))   %only look at the accept rate of previous 30 ones
        i_accepts = 1;
    end
    
    % reset chain ?
%     if(do_chain_resets == true)
%         if(jtrials==100 && jruns+jindexoffset-1>0)
%             p_curr = ps(randi(jruns+jindexoffset-1,1),:); %set of possible initial parameter values
%             likp_curr = mvnpdf(Yobs, Yode, diag(p_curr((end-n_Y+1):end)));
%             prior_curr = mvnpdf(p_curr,muprior,diag(sigmaprior));  %prior of the parameters
%             post_curr = likp_curr * prior_curr;
%             res_curr = [];
%             Sres_curr = [];
%             L_curr = -2*log(post_curr);
%             jtrials = 1;
%             count_chain_reset = count_chain_reset + 1;
%         end
%     end
    
    % update
    if(qa)  % if accept this trial
        p_curr = p_trial;
        L_curr = L_trial;

        ps_hist(ps_hist_index,:) = p_curr;
        ps_hist_index = ps_hist_index + 1;
        if(ps_hist_index > nwindow)
            ps_hist_index = 1;
        end   
    end
    
    
    % save samples
    if(jrungo>0)
            jthin = jthin + 1;
            if(jthin > nthinning)
                ps(jcount+jindexoffset,:) = pReset;
                ps(jcount+jindexoffset,ifFit) = p_curr;
                ps_trial(jcount+jindexoffset,:) = pReset;
                ps_trial(jcount+jindexoffset,ifFit) = p_trial; %if p_trial here is not equal to p_curr, then the previous step rejected.
                chi2s(jcount+jindexoffset) = L_curr;
                chi2s_trial(jcount+jindexoffset) = L_trial;
                acceptance(jcount+jindexoffset) = accept_rate;
                jthin = 1;
                jcount = jcount + 1;
            end
        end
    
    
    jrungo = jrungo + 1;
end

arWaitbar(-1);
fprintf('done (%s, %i chain resets) \n', secToHMS(toc), count_chain_reset);
mcmc_toc = toc;
%ar.p = pReset;


function [post] = postLik(ptmp)
    
         inputp_tmp = startp;
         inputp_tmp(ifFit) = ptmp;
         likp_tmp = mvnpdf(Yobs, Yode, diag(inputp_tmp((end-n_Y+1):end))); %put the var(y) in the last
         prior_tmp = mvnpdf(inputp_tmp,muprior,diag(sigmaprior));  %prior of the parameters
         post = likp_tmp * prior_tmp;
end

function [mu, covar] = mcmc_adaptive(ptmp, index)
        qnonnan = ~isnan(ps_hist(:,1,index));
        if(jrungo<=0 | sum(qnonnan)<=1)
            mu = ptmp;
            covar = eye(length(ptmp)) * Cfactor;
        elseif(jrungo>0 & sum(qnonnan)>1)
            
            CAdapt = cov(ps_hist(qnonnan,:,index));
            mu = ptmp;
            covar = CAdapt;
        end
    end


end

