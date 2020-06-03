% By Chen Liao, Memorial Sloan Kettering Cancer Center, Oct. 14, 2019

% This file implements a coarse-grained ecology model for acetate-mediated
% cross-feeding between two Escherichia coli mutants, a glucose specialist
% (CV103) and an acetate specialist (CV101).

% Briefly, the CV103 mutant has faster glucose uptake rate than the CV101
% mutant. However, it cannot utilize acetate, while CV101 can grow on
% acetate and co-utilize both carbon sources. By secreting acetate, CV103
% creates an additional carbon-source niche for CV101 and the two species
% are thus involved in a one-way cross-feeding interaction.

% The model is simulated in three realistic experimental conditions and
% comapared to measured data. Each simulation-experiment comparison is
% shown in a separate graph.

% Experiment/Figure 1: Growth of CV101 and CV103 in 0.1% glucose batch monoculture
% Experiment/Figure 2: Growth of CV101 and CV103 in 0.0125% glucose batch monoculture at varied acetate concentration
% Experiment/Figure 3: Growth of CV101 and CV103 in 0.00625% glucose chemostat coculture w/ (1 mM) and w/o feeding acetate

% The original data in Figure 1 and 3 were publised in Rosenzweig et al., 1994, Genetics
% The original data in Figure 2 were published in Gudelj et al., 2016, PLoS Computational Biology

%% Constants and options
Density_per_OD420 = 1e9; % linear coefficient between cell density (cell number/mL) and OD420
options_ode15s    = odeset('RelTol', 1e-10, 'AbsTol', 1e-10 * ones(4,1), 'NonNegative', 1:4); % Matlab options
ms = 10; % marker size
lw = 2;  % line width
fs = 16; % font size

%% Model parameters
parameters = [5.52e-13, ... Maximum glucose uptake rate of CV101, mmol/h
              6.28e-13, ... Maximum glucose uptake rate of CV103, mmol/h
              0.01,     ... Michaelis constant for glucose uptake, mM
              26.30,    ... Half-inhibition constant for CV101 growth by acetate, mM
              112.70,   ... Half-inhibition constant for CV103 growth by acetate, mM
              0.87,     ... Half-inhibition constant for acetate uptake by glucose, mM
              4.26e-13, ... Maximum acetate uptake rate of CV101, mmol/h
              0.20,     ... Michaelis constant for acetate uptake by CV101, mM
              2.00e11,  ... Biomass yield of E. coli grown on acetate
              0.33,     ... Acetate leakage fraction (percentage of carbon loss)
              3.00      ... Number of acetate molecules produced per molecule of glucose is consumed
             ];

%% Load in experimental data for comparison to simulations
Exp1Data_CV101   = readtable('data/Experiment1_CV101_monoculture_growth.xlsx');
Exp1Data_CV103   = readtable('data/Experiment1_CV103_monoculture_growth.xlsx');
Exp2Data_CV101   = readtable('data/Experiment2_CV101_acetate_inhibition.xlsx');
Exp2Data_CV103   = readtable('data/Experiment2_CV103_acetate_inhibition.xlsx');
Exp3Data_All     = readtable('data/Experiment3_CV101CV103_coculture.xlsx');

%% Experiment 1: Growth of CV101 and CV103 in glucose-minimum batch monoculture
D               = 0;                    % Dilution rate, 1/h
G_fm            = 0;                    % Feed medium glucose concentration, mM
A_fm            = 0;                    % Feed medium acetate concentration, mM
G_mm            = 180.156;              % Molecular mass of glucose, g/mol
G_init          = (0.1*10)/G_mm*1e3;    % Initial glucose concentration (0.1%, 0.1 g in 100 ml water) in the culture vessel, mM
A_init          = 0;                    % Initial acetate concentration in the culture vessel, mM

% Growth of CV101 monoculture
init_CV101      = 0.024*Density_per_OD420*1e3; % initial population density of CV101, cell number/L
[t_exp1_CV101,y_exp1_CV101] = ode15s(@acetate_mediated_cf_simplified_model, ... Model
                                     [0,40],         ... Time interval of integration (unit: h)
                                     [G_init,init_CV101,0,A_init], ... Initial condition
                                     options_ode15s, ... Matlab options
                                     D,G_fm,A_fm,parameters ... Other parameters
                                    );
                          
% Growth of CV103 monoculture
init_CV103      = 0.012*Density_per_OD420*1e3; % initial population density of CV103, cell number/L
[t_exp1_CV103,y_exp1_CV103] = ode15s(@acetate_mediated_cf_simplified_model, ... Model
                                     [0,40],         ... Time interval of integration (unit: h)
                                     [G_init,0,init_CV103,A_init], ... Initial condition
                                     options_ode15s, ... Matlab options
                                     D,G_fm,A_fm,parameters ... Other parameters
                                    );
                                
figure('Name','Fig. 2B-E');
subplot(2,2,1);
hold on;
plot(Exp1Data_CV101(:,1).Variables, Exp1Data_CV101(:,2).Variables, ...
    'ko', 'MarkerFaceColor', 'r', 'MarkerSize', ms);
plot(t_exp1_CV101, y_exp1_CV101(:,2)/Density_per_OD420/1e3, 'r-', 'LineWidth', lw);
plot(Exp1Data_CV103(:,1).Variables, Exp1Data_CV103(:,2).Variables, ...
    'ko', 'MarkerFaceColor', 'b', 'MarkerSize', ms);
plot(t_exp1_CV103, y_exp1_CV103(:,3)/Density_per_OD420/1e3, 'b-', 'LineWidth', lw);
h=legend('CV101 monoculture: exp','CV101 monoculture: sim','CV103 monoculture: exp','CV103 monoculture: sim',...
         'Location','SouthEast');
h.FontSize = 12;
axis square;
box on;
set(gca,'XTick',[0:10:40]);
set(gca,'YTick',10.^[-2:1:1]);
set(gca,'FontSize',fs);
set(gca,'YScale','log');
xlabel('Time (h)');
ylabel('OD_{420}');
axis([0,40,0.01,10]);

subplot(2,2,2);
hold on;
plot(Exp1Data_CV101(:,3).Variables, Exp1Data_CV101(:,4).Variables,...
    'ko', 'MarkerFaceColor', 'r', 'MarkerSize', ms);
plot(t_exp1_CV101, y_exp1_CV101(:,4), 'r-', 'LineWidth', lw);
plot(Exp1Data_CV103(:,3).Variables, Exp1Data_CV103(:,4).Variables,...
    'ko', 'MarkerFaceColor', 'b', 'MarkerSize', ms);
plot(t_exp1_CV103, y_exp1_CV103(:,4), 'b-', 'LineWidth', lw);
h=legend('CV101 monoculture: exp','CV101 monoculture: sim','CV103 monoculture: exp','CV103 monoculture: sim',...
         'Location','NorthWest');
h.FontSize = 12;
axis square;
box on;
set(gca,'XTick',[0:10:40]);
set(gca,'YTick',[0:3:9]);
set(gca,'FontSize',fs);
xlabel('Time (h)');
ylabel('Acetate (mM)');
axis([0,40,0,9]);

%% Experiment 2: Growth inhibition of CV101 and CV103 by acetate in batch monoculture
G_init          = (0.0125*10)/G_mm*1e3;     % Initial glucose concentration (0.0125%, 0.0125 g in 100 ml water) in the culture vessel, mM
A_init          = [0:0.1:100];              % Initial acetate concentration in the culture vessel, mM

J_grow = zeros(2,length(A_init)); % growth rate at each acetate concentration
for i=1:length(A_init)
    [~,J_grow(:,i)] = acetate_mediated_cf_simplified_model(0, [G_init,0,0,A_init(i)], 0, 0, 0, parameters);
end

subplot(2,2,3);
hold on;
plot(Exp2Data_CV101(:,1).Variables / 60.052 * 1e3, ... convert g/L to mM
     Exp2Data_CV101(:,2).Variables/max(Exp2Data_CV101(:,2).Variables), ...
    'ko', 'MarkerFaceColor', 'r', 'MarkerSize', ms);
plot(A_init, J_grow(1,:)/max(J_grow(1,:)), 'r-', 'LineWidth', lw);
plot(Exp2Data_CV103(:,1).Variables / 60.052 * 1e3, ... convert g/L to mM
     Exp2Data_CV103(:,2).Variables/max(Exp2Data_CV103(:,2).Variables), ...
    'ko', 'MarkerFaceColor', 'b', 'MarkerSize', ms);
plot(A_init, J_grow(2,:)/max(J_grow(2,:)), 'b-', 'LineWidth', lw);
h=legend('CV101: exp','CV101: sim','CV103: exp','CV103: sim','Location','SouthEast');
h.FontSize = 12;
axis square;
box on;
set(gca,'XTick',[0:10:40]);
set(gca,'YTick',[0:0.2:1]);
set(gca,'FontSize',fs);
xlabel('Acetate (mM)');
ylabel('Relative growth rate');
axis([0,40,0,1]);

%% Experiment 3: Growth of CV101 and CV103 in chemostat coculture w/ and w/o acetate supplementation
D               = 0.2;                      % Dilution rate, 1/h
G_fm            = (0.00625*10)/G_mm*1e3;    % Feed medium glucose concentration (0.00625%, 0.00625 g in 100 ml water), mM
V_env           = 200*1e-3;                 % Culture medium volume (in the range between 145-210 mL), L
G_init          = 0;                        % Initial glucose concentration (0.00625%) in the culture vessel, mM
A_init          = 0;                        % Initial acetate concentration in the culture vessel, mM
init_total_cell_number = 2e8;               % Initial number of total cells

% Growth of coculture w/o acetate supplementation
A_fm            = 0;                               % Feed medium acetate concentration, mM
init_frac_CV101 = 0.5;                             % Initial fraction of CV101 cells
init_CV101      = init_total_cell_number*init_frac_CV101/V_env;         % initial population density of CV101, cell number/L
init_CV103      = init_total_cell_number*(1-init_frac_CV101)/V_env;     % initial population density of CV103, cell number/L
[t_exp3_wo_acetate,y_exp3_wo_acetate] = ode15s(...
    @acetate_mediated_cf_simplified_model,  ... Model
    [0,1500],         ... Time interval of integration (unit: h)
    [G_init,init_CV101,init_CV103,A_init], ... Initial condition
    options_ode15s, ... Matlab options
    D,G_fm,A_fm,parameters ... Other parameters
   );

% Growth of coculture w acetate supplementation
A_fm            = 1.0;                             % Feed medium acetate concentration, mM
init_frac_CV101 = 0.63;                            % Initial fraction of CV101 cells
init_CV101      = init_total_cell_number*init_frac_CV101/V_env;         % initial population density of CV101, cell number/L
init_CV103      = init_total_cell_number*(1-init_frac_CV101)/V_env;     % initial population density of CV103, cell number/L
[t_exp3_w_acetate,y_exp3_w_acetate] = ode15s(...
    @acetate_mediated_cf_simplified_model,  ... Model
    [0,1500],         ... Time interval of integration (unit: h)
    [G_init,init_CV101,init_CV103,A_init], ... Initial condition
    options_ode15s, ... Matlab options
    D,G_fm,A_fm,parameters ... Other parameters
   );
    
subplot(2,2,4);
hold on;
plot(Exp3Data_All(:,1).Variables, Exp3Data_All(:,2).Variables, ...
    'ko', 'MarkerFaceColor', 'b', 'MarkerSize', ms);
plot(t_exp3_wo_acetate/3.5, y_exp3_wo_acetate(:,3)./(y_exp3_wo_acetate(:,2)+y_exp3_wo_acetate(:,3)), 'b-', 'LineWidth', lw);
plot(Exp3Data_All(:,3).Variables, Exp3Data_All(:,4).Variables, ...
    'ko', 'MarkerFaceColor', 'r', 'MarkerSize', ms);
plot(t_exp3_w_acetate/3.5, y_exp3_w_acetate(:,3)./(y_exp3_w_acetate(:,2)+y_exp3_w_acetate(:,3)), 'r-', 'LineWidth', lw);
h=legend('Coculture w/o acetate: exp','Coculture w/o acetate: sim','Coculture w/ acetate: exp','Coculture w/ acetate: sim',...
         'Location','East');
h.FontSize = 12;
axis square;
box on;
set(gca,'XTick',[0:10:40]);
set(gca,'YTick',[0:0.2:1]);
set(gca,'FontSize',fs);
xlabel('Generation');
ylabel('Relative frequency of CV103');
axis([0,40,0,1]);
