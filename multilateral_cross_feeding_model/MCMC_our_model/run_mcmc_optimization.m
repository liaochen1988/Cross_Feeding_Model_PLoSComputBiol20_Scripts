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

%% data file
data_file = '../data/experiment_data.xls';

%% read initial parameters
tblPara = readtable('../data/parameters.xls','Sheet','byproduct fraction');
auxotroph = {'C';'F';'G';'H';'I';'K';'L';'M';'P';'R';'S';'T';'W';'Y'};
Byp_frac = tblPara{:, 2:end};

params = {
    {'V_dk_k',        5.00e-11,   0}
    {'V_dm_m',        1.76e-10,   0}
    {'V_dr_r',        3.20e-11,   0}
    {'V_dt_t',        2.79e-10,   0}
    {'Km_dh_h',       2.60e-2,    0}
    {'Km_di_i',       2.20e-1,    0}
    {'Km_dr_r',       5.00e-2,    0}
    {'Km_ds_s',       7.50e-1,    0}
    {'Km_dt_t',       5.40e-1,    0}
    {'eta_dc',        1.50e-1,    0}
    {'eta_df',        3.00e-1,    0}
    {'eta_dg',        2.00e-1,    0}
    {'eta_dh',        2.00e-1,    0}
    {'eta_di',        4.00e-1,    0}
    {'eta_dk',        0,    0}
    {'eta_dl',        2.00e-1,    0}
    {'eta_dm',        0,    0}
    {'eta_dp',        1.00e-1,    0}
    {'eta_dr',        0,    0}
    {'eta_ds',        7.50e-1,    0}
    {'eta_dt',        2.00e-1,    0}
    {'eta_dw',        5.00e-2,    0}
    {'eta_dy',        1.00e-1,    0}
    };

read_from_data = 1;
if (~read_from_data)
    for j=1:size(Byp_frac,2)
        for i=1:size(Byp_frac,1)
            params{end+1} = {strcat('varphi_',auxotroph{i},'_',auxotroph{j}), Byp_frac(i,j), 0, 1};
        end
    end
end

model.ssfun  = @residual_function;
model.sigma2 = 0.01^2;
model.N = 742;

if (burnin)
    options.nsimu = 10000;
else
    options.nsimu = 20000;
end
options.updatesigma = 1;
options.verbosity   = 1;

if (burnin)
    [res,chain,s2chain] = mcmcrun(model,data_file,params,options);
else
    [res,chain,s2chain] = mcmcrun(model,data_file,params,options,previous_results);
end
mcmcplot(chain,[],res,'chainpanel');

if (burnin)
    save('burnin.mat');
else
    save('Example3_MCMC_20000_smape.mat');
end
