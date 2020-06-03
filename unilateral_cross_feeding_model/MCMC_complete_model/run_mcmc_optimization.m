% This script runs parameter sensitivity analysis of acetate-mediated
% cross-feeding model (complete version)
%
% Last updated by Chen Liao, June 2, 2020

addpath('../../auxiliary_functions/mcmc/');
addpath('../');

burnin = true;
if (~burnin)
    load('burnin.mat');
    previous_results = res;
    burnin = false;
end

%% read experimental data
Exp1Data_CV101   = readtable('../data/Experiment1_CV101_monoculture_growth.xlsx');
Exp1Data_CV103   = readtable('../data/Experiment1_CV103_monoculture_growth.xlsx');
Exp2Data_CV101   = readtable('../data/Experiment2_CV101_acetate_inhibition.xlsx');
Exp2Data_CV103   = readtable('../data/Experiment2_CV103_acetate_inhibition.xlsx');
Exp3Data_All     = readtable('../data/Experiment3_CV101CV103_coculture.xlsx');
ExpDataAll = struct('Exp1Data_CV101',Exp1Data_CV101,...
                    'Exp1Data_CV103',Exp1Data_CV103,...
                    'Exp2Data_CV101',Exp2Data_CV101,...
                    'Exp2Data_CV103',Exp2Data_CV103,...
                    'Exp3Data_All',  Exp3Data_All);

% parameters
params = {
    {'V_1g',        5.52e-13,   0}
    {'V_3g',        6.28e-13,   0}
    {'V_1a',        2.84e-13,   0}
    {'C_1g',        0.87,       0}
    {'phi_a',       0.75,       0,      1}
    {'varphi_a',    0.45,       0,      1}
    {'gamma_g',     1e13,        0}
    {'gamma_a',     3e11,       0}
    {'I_1a',        39.45,       0}
    {'I_3a',        112.7,      0}
    };

model.ssfun  = @residual_function;
model.sigma2 = 0.01^2;
model.N = 72;

if (burnin)
    options.nsimu = 10000;
else
    options.nsimu = 20000;
end
options.updatesigma = 1;
options.verbosity   = 1;

if (burnin)
    [res,chain,s2chain] = mcmcrun(model,ExpDataAll,params,options);
else
    [res,chain,s2chain] = mcmcrun(model,ExpDataAll,params,options,previous_results);
end
mcmcplot(chain,[],res,'chainpanel');

part_i = 1;
if (burnin)
    save('burnin.mat');
else
    save(sprintf('../output/Example1_complete_MCMC_20000_smape_part%d.mat', part_i), 'chain', 'res');
end