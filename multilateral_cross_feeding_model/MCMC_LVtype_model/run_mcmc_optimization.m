% This script runs parameter sensitivity analysis of Mee et al. 2014 model
% Last updated by Chen Liao, April 29, 2020

addpath('../../auxiliary_functions/mcmc');
addpath('../');

burnin = true;
if (~burnin)
    load('burnin.mat');
    previous_results = res;
    burnin = false;
end

%% read experimental data
tblPairwiseCoculture = readtable('../data/experiment_data.xls','Sheet','pairwise coculture');
obsFC = tblPairwiseCoculture{:,[2,4]};

%% read initial parameters
coefs_Mee2014 = readtable('../data/parameters.xls','Sheet','coefficients of Mee 2014 model');
k_Mee2014       = 1e9; % carrying capactiy
beta_Mee2014    = 1;   % Michaelis constant

params = cell(size(coefs_Mee2014,1)*2+2, 1);
params{1} = {'k', 1e9, 0};
params{2} = {'beta', 1, 0};
for i=1:height(coefs_Mee2014)
    params{2*(i-1)+3} = {strcat('C12_',coefs_Mee2014.Strain1{i},'_',coefs_Mee2014.Strain2{i}), coefs_Mee2014.C12(i), 0};
    params{2*(i-1)+4} = {strcat('C21_',coefs_Mee2014.Strain1{i},'_',coefs_Mee2014.Strain2{i}), coefs_Mee2014.C21(i), 0};
end

model.ssfun  = @residual_function;
model.sigma2 = 0.01^2;
model.N = size(obsFC,1)*size(obsFC,2);

if (burnin)
    options.nsimu = 10000;
else
    options.nsimu = 100000;
end
options.updatesigma = 1;
options.verbosity   = 1;

if (burnin)
    [res,chain,s2chain] = mcmcrun(model,obsFC,params,options);
else
    [res,chain,s2chain] = mcmcrun(model,obsFC,params,options,previous_results);
end
mcmcplot(chain,[],res,'chainpanel');

if (burnin)
    save('burnin.mat');
else
    save('Fig5B_MCMC_run_100000_steps_smape.mat');
end