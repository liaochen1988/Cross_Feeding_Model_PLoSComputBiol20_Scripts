% This script estimates posterior distribution of free parameters in the
% amino-acid-mediated cross-feeding system
%
% Last updated by Chen Liao, June 2, 2020

addpath('../');
addpath('../../auxiliary_functions');
addpath('../../auxiliary_functions/mcmc');

burnin = true;
if (~burnin)
    load('burnin.mat');
    previous_results = res;
    burnin = false;
end

%% read experimental data
Exp1Data_DeltaK     = readtable('../data/Experiment1_LysAux_monoculture_growth.xlsx');
Exp2Data_DeltaL     = readtable('../data/Experiment2_LeuAux_monoculture_growth.xlsx');
Exp3Data_DeltaKL    = readtable('../data/Experiment3_LysLeuAux_coculture_growth.xlsx');
ExpDataAll = struct('Exp1Data_DeltaK',Exp1Data_DeltaK,...
                    'Exp2Data_DeltaL',Exp2Data_DeltaL,...
                    'Exp3Data_DeltaKL',Exp3Data_DeltaKL);

% parameters
params = {
    {'Km_g',           1.75,       0}
    {'Vmax_dk_k',      8.35e-14,   0}
    {'Vmax_dl_k',      1e-16,   0}
    {'Vmax_dk_l',      1e-16,   0}
    {'Vmax_dl_l',      1.22e-13,   0}
    {'phi_dl_k',       0.01,        0,     1}
    {'phi_dk_l',       0.01,        0,     1}
    {'varphi_dl_k',    1.39e-1,    0,      1}
    {'varphi_dk_l',    3.20e-2,    0,      1}
    {'gamma_k',        5.72e12,    0}
    {'gamma_l',        2.53e12,    0}
    {'eta_dk',         0.102,      0}
    {'eta_dl',         4.00e-4,    0}
    };

model.ssfun  = @residual_function;
model.sigma2 = 0.01^2;
model.N = 129;

if (burnin)
    options.nsimu = 10000;
else
    options.nsimu = 100000;
end
options.updatesigma = 1;
options.verbosity   = 1;

if (burnin)
    [res,chain,s2chain] = mcmcrun(model,ExpDataAll,params,options);
else
    [res,chain,s2chain] = mcmcrun(model,ExpDataAll,params,options,previous_results);
end
mcmcplot(chain,[],res,'chainpanel');

if (burnin)
    save('burnin.mat');
else
    save('Example2_complete_MCMC_100000_smape.mat','res','chain');
end