% By Chen Liao, Memorial Sloan Kettering Cancer Center, June, 3 2020
%
% This file draws scatter plot of observed data vs simulation for our model
% as well as Lotka-Volterra-type model (parameters adopted from Mee et al.
% PNAS 2014 or optimized through MCMC)

auxotroph = {'C';'F';'G';'H';'I';'K';'L';'M';'P';'R';'S';'T';'W';'Y'};

%% Load observed data
tblPairwiseCoculture = readtable('data/experiment_data.xls','Sheet','pairwise coculture');
obsFC = tblPairwiseCoculture{:,[2,4]};

figure();

%% Simulate pairwise coculture using our manually curated parameters
% Before coculture experiments, all cell concentrations were adjusted to 10^7 cells per mL using M9 media.
% Coculture growth was performed by equal-volume inoculation of each strain at a seeding density of 107 cells per mL.
% All two-member and three-member cocultures were grown in 200 uL of M9-glucose media

tblPara = readtable('data/parameters.xls','Sheet','kinetic rates and yields');
Vmax_g   = tblPara.Vmax_g;     % Maximum glucose uptake rates, umol/h
Km_g     = tblPara.Km_g;       % Michaelis constants for glucose uptake, uM
Vmax_aa  = tblPara.Vmax_aa;    % Maximum amino acid uptake rates, umol/h
Km_aa    = tblPara.Km_aa;      % Michaelis constants for amino acid uptake, uM
Yield_g  = tblPara.Yield_g;    % Biomass yields of E. coli auxotrophs growing on glucose, 1/umol
Yield_aa = tblPara.Yield_aa;   % Biomass yields of E. coli auxotrophs growing on auxotrophic amino acids, 1/umol
Rcarbon  = tblPara.Rcarbon;    % Number of amino acid molecules produced per molecule of glucose is consumed
Cdr      = tblPara.Cdr;       % Cell death rate, 1/h

tblPara = readtable('data/parameters.xls','Sheet','byproduct fraction');
Byp_frac = tblPara{:, 2:end};

initial_glucose = 0.002/180.156*1e6/1e-3; % uM, glucose concentration 0.2% wt/vol (0.2 g in 100g water)
initial_cell_density = 1e7*1e3;           % 1e7 per mL (convert mL to L)

simFC_this_paper = zeros(size(obsFC)); % simulated fold change of pairwise cocultures
for i=1:height(tblPairwiseCoculture)
    first_strain_index = find(strcmp(auxotroph, tblPairwiseCoculture.Strain1(i)));
    second_strain_index = find(strcmp(auxotroph, tblPairwiseCoculture.Strain2(i)));

    % setup initial condition
    y0 = zeros(1+3*14,1);
    y0(1) = initial_glucose;
    y0([1+first_strain_index,15+first_strain_index,1+second_strain_index,15+second_strain_index]) = initial_cell_density;
    
    % integrate
    [~,y] = ode15s(@multi_aa_auxotroph_cf_model, ... Model
        [0,84],  ... Time interval of integration (unit: h)
        y0,      ... Initial condition
        odeset('RelTol',1e-6,'AbsTol',1e-6*ones(1,1+3*14),'NonNegative',ones(1,1+3*14)), ... Matlab options
        Vmax_g, Km_g, Vmax_aa, Km_aa, Yield_g, Yield_aa, Rcarbon, Byp_frac, zeros(size(Cdr)) ... Model parameters
        );
    simFC_this_paper(i,:) = [y(end,1+first_strain_index), y(end,1+second_strain_index)]./initial_cell_density;
end

hold on;
plot(obsFC(:),simFC_this_paper(:),'ko','MarkerSize',8,'MarkerEdgeColor','r');

%% Simulate pairwise coculture using parameters given by Mee et al. 2014
coefs_Mee2014 = readtable('data/parameters.xls','Sheet','coefficients of Mee 2014 model');
k_Mee2014       = 1e9; % carrying capactiy
beta_Mee2014    = 1;   % Michaelis constant
C12_Mee2014     = coefs_Mee2014{:,'C12'};
C21_Mee2014     = coefs_Mee2014{:,'C21'};

simFC_Mee2014 = zeros(size(obsFC));
option = odeset('RelTol',1e-6,'AbsTol',1e-6*ones(1,2),'NonNegative',ones(1,2));
for k=1:size(obsFC,1)
    [~,y] = ode15s(@two_member_gLV_model, [0,84], [1e7,1e7], option, [C12_Mee2014(k),C21_Mee2014(k)], beta_Mee2014, k_Mee2014);
    simFC_Mee2014(k,1) = y(end,1)/1e7;
    simFC_Mee2014(k,2) = y(end,2)/1e7;
end

plot(obsFC(:),simFC_Mee2014(:),'ko','MarkerSize',8,'MarkerEdgeColor','b');

%% Simulate pairwise coculture using parameters estimated from MCMC
load('output/Fig5B_pcorr.mat');
[~,max_pcorr_index] = max(pcorr_MCMC);
load('output/Fig5B_MCMC_run_100000_steps_smape.mat');
k_MCMC       = chain(max_pcorr_index,1); % carrying capactiy
beta_MCMC    = chain(max_pcorr_index,2); % Michaelis constant
C12_MCMC     = chain(max_pcorr_index,3:2:end);
C21_MCMC     = chain(max_pcorr_index,4:2:end);

simFC_MCMC = zeros(size(obsFC));
for k=1:size(obsFC,1)
    [~,y] = ode15s(@two_member_gLV_model, [0,84], [1e7,1e7], option, [C12_MCMC(k),C21_MCMC(k)], beta_MCMC, k_MCMC);
    simFC_MCMC(k,1) = y(end,1)/1e7;
    simFC_MCMC(k,2) = y(end,2)/1e7;
end

plot(obsFC(:),simFC_MCMC(:),'ko','MarkerSize',8,'MarkerEdgeColor','g');
plot([1e-2,1e2],[1e-2,1e2],'k--');
axis([1e-2,1e2,1e-2,1e2]);
xlabel('Observed fold change');
ylabel('Predicted fold change');
axis square;
box on;
set(gca,'XScale','log');
set(gca,'YScale','log');
legend('Our model','Mee et al. 2014','MCMC','Location','SouthWest');