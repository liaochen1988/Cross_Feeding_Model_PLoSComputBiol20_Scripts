% By Chen Liao, Memorial Sloan Kettering Cancer Center, May 16, 2020
%
% This file draws barplot to compare accuracy of data fitting between
% different models/parameters

%% Load observed data
tblPairwiseCoculture = readtable('data/experiment_data.xls','Sheet','pairwise coculture');
obsFC = tblPairwiseCoculture{:,[2,4]};

%% Simulate pairwise coculture using our manually curated parameters
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
auxotroph = {'C';'F';'G';'H';'I';'K';'L';'M';'P';'R';'S';'T';'W';'Y'};

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

pcorr_this_paper = corrcoef(obsFC(:), simFC_this_paper(:)); % Pearson correlation
pcorr_this_paper = pcorr_this_paper(1,2);

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

pcorr_Mee2014 = corrcoef(obsFC(:), simFC_Mee2014(:)); % Pearson correlation
pcorr_Mee2014 = pcorr_Mee2014(1,2);

%% Simulate pairwise coculture using parameters optimized by MCMC (MCMC sample that gives the highest PCC)
load('output/Fig5B_MCMC_run_100000_steps_smape.mat');
k_MCMC       = chain(:,1); % carrying capactiy
beta_MCMC    = chain(:,2); % Michaelis constant
C12_MCMC     = chain(:,3:2:end);
C21_MCMC     = chain(:,4:2:end);
assert(size(k_MCMC,1)==100000);

pcorr_MCMC = zeros(size(k_MCMC));
parfor i=1:size(k_MCMC,1)
    simFC_MCMC = zeros(size(obsFC));
    for k=1:size(obsFC,1)
        [~,y] = ode15s(@two_member_gLV_model, [0,84], [1e7,1e7], option, [C12_MCMC(i,k),C21_MCMC(i,k)], beta_MCMC(i), k_MCMC(i));
        simFC_MCMC(k,1) = y(end,1)/1e7;
        simFC_MCMC(k,2) = y(end,2)/1e7;
    end
    pcorr_tmp = corrcoef(obsFC(:), simFC_MCMC(:)); % Pearson correlation
    pcorr_MCMC(i) = pcorr_tmp(1,2);
end

save('output/pcorr.mat', 'pcorr_this_paper', 'pcorr_Mee2014', 'pcorr_MCMC');

%% Barplot
load('output/pcorr.mat');
addpath('../Violinplot-Matlab-master');

figure();

subplot(1,3,1);
box on;
bar(1, pcorr_this_paper);
ylim([0,1]);
set(gca,'YTick',[0:0.5:1]);
set(gca,'XTickLabel',{'Our model'});
xlabel('');
ylabel('Pearson correlation');

subplot(1,3,2);
box on;
bar(2, pcorr_Mee2014);
ylim([0,1]);
set(gca,'YTick',[0:0.5:1]);
set(gca,'XTickLabel',{'Mee et al. 2014'});
xlabel('');
ylabel('Pearson correlation');

subplot(1,3,3);
box on;
bar(3, max(pcorr_MCMC));
ylim([0,1]);
set(gca,'YTick',[0:0.5:1]);
set(gca,'XTickLabel',{'MCMC'});
xlabel('');
ylabel('Pearson correlation');

