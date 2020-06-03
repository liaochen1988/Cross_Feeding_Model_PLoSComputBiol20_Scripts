% By Chen Liao, Memorial Sloan Kettering Cancer Center, June 2, 2019
%
% The file simulates CV101-CV103 coculture up to 200 generations

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
Exp3Data_All     = readtable('data/Experiment3_CV101CV103_coculture.xlsx');

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
    
figure('Name','FigS2');
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
set(gca,'XTick',[0:50:200]);
set(gca,'YTick',[0:0.2:1]);
set(gca,'FontSize',fs);
xlabel('Generation');
ylabel('Relative frequency of CV103');
axis([0,200,0,1]);
