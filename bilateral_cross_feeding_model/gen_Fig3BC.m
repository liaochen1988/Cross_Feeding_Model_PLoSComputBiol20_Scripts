% By Chen Liao, Memorial Sloan Kettering Cancer Center, May. 16, 2020

% This file implements a coarse-grained ecology model for bilateral
% cross-feeding between two engineered Escherichia coli auxotrophs, one
% unable to synthesize leucine and the other unable to synthesize lysine.

% Briefly, the two mutants cooperate by compensating for each other?s
% metabolic deficiency: the lysine auxotroph (K) secretes leucine that can
% be utilized by the leucine auxotroph (L), which in return facilitate
% growth of the lysine auxotroph by secreting lysine to the environment

% The model is simulated in both monoculture and coculture conditions and
% comapared to experimental data. Each simulation-experiment comparison is
% shown in a separate graph.

% Experiment/Figure 1: Growth of the lysine auxotroph in glucose batch monoculture
% Experiment/Figure 2: Growth of the leucine auxotroph in glucose batch monoculture
% Experiment/Figure 3: Growth of the lysine auxotroph and the leucine auxotroph in glucose batch coculture

% The original data in Figures 1-3 were publised in Zhang & Reed, 2014, PLoS ONE

addpath('../auxiliary_functions');

%% Constants and options
Density_per_OD600 = 5e8.*[1.6,1]; % linear coefficient between cell density (number of cells /mL) of [lysine,leucine] auxotroph and OD600
options_ode15s    = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(7,1),'NonNegative',[1:7]); % Matlab options
ms = 8; % marker size
lw = 1; % line width
mfc_red = [240,128,128]/255;  % marker face color for red
mfc_blue = [100,149,237]/255; % marker face color for blue
mfc_gray = [128,128,128]/255; % marker face color for gray
fs = 16; % font size

%% Model parameters 
parameters = [...
    3.61e-12; ... Maximum glucose uptake rate, mmol/h
    1.75;     ... Michaelis constant for glucose uptake, mM
    8.34e-14; ... Maximum lysine uptake rate by the lysine auxotroph, mmol/h
    1.22e-13; ... Maximum leucine uptake rate by the leucine axutroph, mmol/h
    5.00e-3;  ... Michaelis constant for lysine uptake, mM
    1.07e-3;  ... Michaelis constant for leucine uptake, mM
    1.00;     ... Number of lysine molecules produced per molecule of glucose is consumed
    1.00;     ... Number of leucine molecules produced per molecule of glucose is consumed
    3.20e-3;  ... Fraction of glucose leaked in form of leucine by the lysine auxotroph
    1.39e-2;  ... Fraction of glucose leaked in form of lysine by the leucine auxotroph
    3.00e11;  ... Biomass yield of E. coli grown on glucose, 1/mmol
    5.72e12;  ... Biomass yield of E. coli grown on lysine, 1/mmol
    2.53e12;  ... Biomass yield of E. coli grown on leucine, 1/mmol
    0.10;     ... Mortality rate constant of the lysine auxotroph, 1/h
    4e-4      ... Mortality rate constant of the leusine auxotroph, 1/h
    ];
         
%% Load in experimental data for comparison to simulations
Exp1Data_DeltaK     = readtable('data/Experiment1_LysAux_monoculture_growth.xlsx');
Exp2Data_DeltaL     = readtable('data/Experiment2_LeuAux_monoculture_growth.xlsx');
Exp3Data_DeltaKL    = readtable('data/Experiment3_LysLeuAux_coculture_growth.xlsx');

figure('Name','Fig3BC');

%% Experiment 1: Growth of the lysine auxotroph in glucose-minimum batch monoculture
D       = 0;    % Dilution rate, 1/h
G_fm    = 0;    % Feed medium glucose concentration, mM
K_fm    = 0;    % Feed medium lysine concentration, mM
L_fm    = 0;    % Feed medium leucine concentration, mM
G_init  = Exp1Data_DeltaK(1,4).Variables;                   % Initial glucose concentration, mM
OD_init = 0.01;                                             % Initial OD
DK_init = OD_init*Density_per_OD600(1)*1e3;                 % Initial cell density of lysine auxotroph, number of cells/L
K_mm    = 146.19;                                           % Molecular mass of lysine, g/mol
K_init  = (Exp1Data_DeltaK(1,6).Variables*1e-3)/K_mm*1e3;   % Initial lysine concentration, mM
[tK_s,yK_s] = ode15s(@lys_leu_bilateral_cf_simplified_model,           ... Model
    [0,15],                                ... Time interval of integration (unit: h)
    [G_init,DK_init,0,DK_init,0,K_init,0], ... Initial condition
    options_ode15s,                        ... Matlab options
    D,G_fm,K_fm,L_fm,parameters);

subplot(3,2,1);
hold on;
plot(Exp1Data_DeltaK(:,1).Variables,Exp1Data_DeltaK(:,2).Variables,'ko','MarkerFaceColor', mfc_red, 'MarkerSize', ms);
plot(tK_s,yK_s(:,2)/Density_per_OD600(1)/1e3,'r-','LineWidth',lw);
axis([0,15,0,0.5]);
set(gca,'XTick',[0:5:15]);
set(gca,'YTick',[0:0.1:0.5]);
set(gca,'FontSize',fs);
axis square;
box on;
xlabel('Time (h)');
ylabel('OD_{600}');

subplot(3,2,3);
hold on;
plot(Exp1Data_DeltaK(:,3).Variables,Exp1Data_DeltaK(:,4).Variables,'ko','MarkerFaceColor', mfc_red, 'MarkerSize', ms);
plot(tK_s,yK_s(:,1),'r-','LineWidth',lw);
axis([0,15,0,15]);
set(gca,'XTick',[0:5:15]);
set(gca,'YTick',[0:5:15]);
set(gca,'FontSize',fs);
axis square;
box on;
xlabel('Time (h)');
ylabel('Glucose (mM)');

subplot(3,2,5);
hold on;
plot(Exp1Data_DeltaK(:,5).Variables,Exp1Data_DeltaK(:,6).Variables,'ko','MarkerFaceColor', mfc_red, 'MarkerSize', ms);
plot(tK_s,yK_s(:,6)*K_mm,'r-','LineWidth',lw);
axis([0,15,0,10]);
set(gca,'XTick',[0:5:15]);
set(gca,'YTick',[0:2:10]);
set(gca,'FontSize',fs);
axis square;
box on;
xlabel('Time (h)');
ylabel('Lysine (mg/L)');
     
%% Experiment 2: Growth of the leucine auxotroph in glucose-minimum batch monoculture
D       = 0;    % Dilution rate, 1/h
G_fm    = 0;    % Feed medium glucose concentration, mM
K_fm    = 0;    % Feed medium lysine concentration, mM
L_fm    = 0;    % Feed medium leucine concentration, mM
G_init  = Exp2Data_DeltaL(1,4).Variables;                   % Initial glucose concentration, mM
OD_init = 0.01;                                             % Initial OD
DL_init = OD_init*Density_per_OD600(2)*1e3;                 % Initial cell density of leucine auxotroph, number of cells/L
L_mm    = 131.17;                                           % Molecular mass of leucine, g/mol
L_init  = (Exp2Data_DeltaL(1,6).Variables*1e-3)/L_mm*1e3;   % Initial leucine concentration, mM
[tL_s,yL_s] = ode15s(@lys_leu_bilateral_cf_simplified_model,           ... Model
    [0,15],                                ... Time interval of integration (unit: h)
    [G_init,0,DL_init,0,DL_init,0,L_init], ... Initial condition
    options_ode15s,                        ... Matlab options
    D,G_fm,K_fm,L_fm,parameters);

subplot(3,2,1);
hold on;
plot(Exp2Data_DeltaL(:,1).Variables,Exp2Data_DeltaL(:,2).Variables,'ko','MarkerFaceColor', mfc_blue, 'MarkerSize', ms);
plot(tL_s,yL_s(:,3)/Density_per_OD600(2)/1e3,'b-','LineWidth',lw);
%h=legend('Lysine auxotroph: exp','Lysine auxotroph: sim','Leucine auxotroph: exp','Leucine auxotroph: sim','Location','Southoutside');
%h.FontSize=12;
axis([0,15,0,0.5]);
set(gca,'XTick',[0:5:15]);
set(gca,'YTick',[0:0.1:0.5]);
set(gca,'FontSize',fs);
axis square;
box on;
xlabel('Time (h)');
ylabel('OD_{600}');

subplot(3,2,3);
hold on;
plot(Exp2Data_DeltaL(:,3).Variables,Exp2Data_DeltaL(:,4).Variables,'ko','MarkerFaceColor', mfc_blue, 'MarkerSize', ms);
plot(tL_s,yL_s(:,1),'b-','LineWidth',lw);
%h=legend('Lysine auxotroph: exp','Lysine auxotroph: sim','Leucine auxotroph: exp','Leucine auxotroph: sim','Location','Southoutside');
%h.FontSize=12;
axis([0,15,0,15]);
set(gca,'XTick',[0:5:15]);
set(gca,'YTick',[0:5:15]);
set(gca,'FontSize',fs);
axis square;
box on;
xlabel('Time (h)');
ylabel('Glucose (mM)');

subplot(3,2,5);
hold on;
plot(Exp2Data_DeltaL(:,5).Variables,Exp2Data_DeltaL(:,6).Variables,'ko','MarkerFaceColor', mfc_blue, 'MarkerSize', ms);
plot(tL_s,yL_s(:,7)*L_mm,'b-','LineWidth',lw);
h=legend('Lysine auxotroph: exp','Lysine auxotroph: sim','Leucine auxotroph: exp','Leucine auxotroph: sim','Location','bestoutside');
h.FontSize=12;
axis([0,15,0,15]);
set(gca,'XTick',[0:5:15]);
set(gca,'YTick',[0:5:15]);
set(gca,'FontSize',fs);
axis square;
box on;
xlabel('Time (h)');
ylabel('Leucine (mg/L)');

%% Experiment 3: Growth of the lysine and leucine auxotrophs in glucose-minimum batch coculture
G_init  = Exp3Data_DeltaKL(1,4).Variables;                   % Initial glucose concentration, mM
[tKL_s,yKL_s] = ode15s(@lys_leu_bilateral_cf_simplified_model,           ... Model
    [0,80],                                ... Time interval of integration (unit: h)
    [G_init,DK_init,DL_init,DK_init,DL_init,0,0], ... Initial condition
    options_ode15s,                        ... Matlab options
    D,G_fm,K_fm,L_fm,parameters);

subplot(3,2,2);
hold on;
plot(Exp3Data_DeltaKL(:,1).Variables,Exp3Data_DeltaKL(:,2).Variables,'ko','MarkerFaceColor', mfc_gray, 'MarkerSize', ms);
plot(tKL_s,yKL_s(:,2)/Density_per_OD600(1)/1e3+yKL_s(:,3)/Density_per_OD600(2)/1e3,'k-','LineWidth',lw);
%h=legend('Coculture: exp','Coculture: sim','Location','Southoutside');
%h.FontSize=12;
axis([0,80,0,0.5]);
set(gca,'XTick',[0:20:80]);
set(gca,'YTick',[0:0.1:0.5]);
set(gca,'FontSize',fs);
axis square;
box on;
xlabel('Time (h)');
ylabel('OD_{600}');

subplot(3,2,4);
hold on;
plot(Exp3Data_DeltaKL(:,3).Variables,Exp3Data_DeltaKL(:,4).Variables,'ko','MarkerFaceColor', mfc_gray, 'MarkerSize', ms);
plot(tKL_s,yKL_s(:,1),'k-','LineWidth',lw);
%h=legend('Coculture: exp','Coculture: sim','Location','Southoutside');
%h.FontSize=12;
axis([0,80,0,15]);
set(gca,'XTick',[0:20:80]);
set(gca,'YTick',[0:5:15]);
set(gca,'FontSize',fs);
axis square;
box on;
xlabel('Time (h)');
ylabel('Glucose (mM)');

subplot(3,2,6);
hold on;
plot(Exp3Data_DeltaKL(:,5).Variables,Exp3Data_DeltaKL(:,6).Variables,'ko','MarkerFaceColor', mfc_gray, 'MarkerSize', ms);
plot(tKL_s,yKL_s(:,4)./yKL_s(:,5),'k-','LineWidth',lw);
h=legend('Coculture: exp','Coculture: sim','Location','bestoutside');
h.FontSize=12;
axis([0,80,0,3]);
set(gca,'XTick',[0:20:80]);
set(gca,'YTick',[0:1:3]);
set(gca,'FontSize',fs);
axis square;
box on;
xlabel('Time (h)');
ylabel('\DeltaLys/\DeltaLeu');

y_obs = [y_obs;...
         remove_NaN(Exp3Data_DeltaKL(:,2).Variables);...
         remove_NaN(Exp3Data_DeltaKL(:,4).Variables);...
         remove_NaN(Exp3Data_DeltaKL(:,6).Variables)];
y_sim_simplified = [y_sim_simplified;...
         pchip(tKL_s, yKL_s(:,2)/Density_per_OD600(1)/1e3+yKL_s(:,3)/Density_per_OD600(2)/1e3, remove_NaN(Exp3Data_DeltaKL(:,1).Variables));...
         pchip(tKL_s, yKL_s(:,1), remove_NaN(Exp3Data_DeltaKL(:,3).Variables));...
         pchip(tKL_s, yKL_s(:,4)./yKL_s(:,5), remove_NaN(Exp3Data_DeltaKL(:,5).Variables))];
y_sim_full = [y_sim_full;...
         pchip(tKL_f, yKL_f(:,2)/Density_per_OD600(1)/1e3+yKL_f(:,3)/Density_per_OD600(2)/1e3, remove_NaN(Exp3Data_DeltaKL(:,1).Variables));...
         pchip(tKL_f, yKL_f(:,1), remove_NaN(Exp3Data_DeltaKL(:,3).Variables));...
         pchip(tKL_f, yKL_f(:,4)./yKL_f(:,5), remove_NaN(Exp3Data_DeltaKL(:,5).Variables))];
          