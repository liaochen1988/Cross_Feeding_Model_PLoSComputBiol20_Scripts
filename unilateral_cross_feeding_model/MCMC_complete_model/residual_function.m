function err = residual_function(params,data)

%% unpack params
parameters = [params(1), ... Maximum glucose uptake rate of CV101, mmol/h
              params(2), ... Maximum glucose uptake rate of CV103, mmol/h
              0.01,      ... Michaelis constant for glucose uptake, mM
              params(9), ... Half-inhibition constant for CV101 growth by acetate, mM
              params(10),... Half-inhibition constant for CV103 growth by acetate, mM
              params(4), ... Half-inhibition constant for acetate uptake by glucose, mM
              params(3), ... Maximum acetate uptake rate of CV101, mmol/h
              0.20,      ... Michaelis constant for acetate uptake by CV101, mM
              params(7), ... Biomass yield of E. coli grown on non-explicitly-modeled building blocks
              params(8), ... Biomass yield of E. coli grown on acetate
              params(5), ... Fraction of glucose influx allocated to biosynthesis of acetate
              params(6), ... Fraction of acetate leakage
              3.00       ... Number of acetate molecules produced per molecule of glucose is consumed
             ];

%% Get data
Exp1Data_CV101 = data.Exp1Data_CV101;
Exp1Data_CV103 = data.Exp1Data_CV103;
Exp2Data_CV101 = data.Exp2Data_CV101;
Exp2Data_CV103 = data.Exp2Data_CV103;
Exp3Data_All   = data.Exp3Data_All;

%% vectorized difference between simulated and experimental data 
y_obs = [];
y_sim = [];
options_ode15s    = odeset('RelTol', 1e-10, 'AbsTol', 1e-10 * ones(4,1), 'NonNegative', 1:4); % Matlab options

%% Experiment 1: Growth of CV101 and CV103 in glucose-minimum batch monoculture
D               = 0;                    % Dilution rate, 1/h
G_fm            = 0;                    % Feed medium glucose concentration, mM
A_fm            = 0;                    % Feed medium acetate concentration, mM
G_mm            = 180.156;              % Molecular mass of glucose, g/mol
G_init          = (0.1*10)/G_mm*1e3;    % Initial glucose concentration (0.1%, 0.1 g in 100 ml water) in the culture vessel, mM
A_init          = 0;                    % Initial acetate concentration in the culture vessel, mM
Density_per_OD420 = 1e9;                % linear coefficient between cell density (cell number/mL) and OD420

% Growth of CV101 monoculture
CV101_init      = 0.024*Density_per_OD420*1e3; % initial population density of CV101, cell number/L
[t_exp1_CV101,y_exp1_CV101] = ode15s(@acetate_mediated_cf_complete_model, ... Model
                                     [0,40],         ... Time interval of integration (unit: h)
                                     [G_init,CV101_init,0,A_init], ... Initial condition
                                     options_ode15s, ... Matlab options
                                     D,G_fm,A_fm,parameters ... Other parameters
                                    );
                          
% Growth of CV103 monoculture
CV103_init      = 0.012*Density_per_OD420*1e3; % initial population density of CV103, cell number/L
[t_exp1_CV103,y_exp1_CV103] = ode15s(@acetate_mediated_cf_complete_model, ... Model
                                     [0,40],         ... Time interval of integration (unit: h)
                                     [G_init,0,CV103_init,A_init], ... Initial condition
                                     options_ode15s, ... Matlab options
                                     D,G_fm,A_fm,parameters ... Other parameters
                                    );
 
y_obs = [y_obs;...
         Exp1Data_CV101(:,2).Variables;...
         Exp1Data_CV103(:,2).Variables;...
         Exp1Data_CV101(:,4).Variables;...
         Exp1Data_CV103(:,4).Variables];
y_sim = [y_sim;...
         pchip(t_exp1_CV101, y_exp1_CV101(:,2)/Density_per_OD420/1e3, Exp1Data_CV101(:,1).Variables);...
         pchip(t_exp1_CV103, y_exp1_CV103(:,3)/Density_per_OD420/1e3, Exp1Data_CV103(:,1).Variables);...
         pchip(t_exp1_CV101, y_exp1_CV101(:,4), Exp1Data_CV101(:,3).Variables);...
         pchip(t_exp1_CV103, y_exp1_CV103(:,4), Exp1Data_CV103(:,3).Variables)];

ndata_EXP1 = length(y_obs);

%% Experiment 2: Growth inhibition of CV101 and CV103 by acetate in batch monoculture
G_init          = (0.0125*10)/G_mm*1e3;     % Initial glucose concentration (0.0125%, 0.0125 g in 100 ml water) in the culture vessel, mM
A_init          = [0:0.1:100];              % Initial acetate concentration in the culture vessel, mM

J_grow = zeros(2,length(A_init)); % growth rate at each acetate concentration
for i=1:length(A_init)
    [~,J_grow(:,i)] = acetate_mediated_cf_complete_model(0, [G_init,0,0,A_init(i)], 0, 0, 0, parameters);
end

y_obs = [y_obs;...
         Exp2Data_CV101(:,2).Variables/max(Exp2Data_CV101(:,2).Variables);...
         Exp2Data_CV103(:,2).Variables/max(Exp2Data_CV103(:,2).Variables)];
y_sim = [y_sim;...
         pchip(A_init, J_grow(1,:)/max(J_grow(1,:)), Exp2Data_CV101(:,1).Variables / 60.052 * 1e3);...
         pchip(A_init, J_grow(2,:)/max(J_grow(2,:)), Exp2Data_CV103(:,1).Variables / 60.052 * 1e3)];

%% Experiment 3: Growth of CV101 and CV103 in chemostat coculture w/ and w/o acetate supplementation
D               = 0.2;                      % Dilution rate, 1/h
G_fm            = (0.00625*10)/G_mm*1e3;    % Feed medium glucose concentration (0.00625%, 0.00625 g in 100 ml water), mM
V_env           = 200*1e-3;                 % Culture medium volume (in the range between 145-210 mL), L
G_init          = 0;                        % Initial glucose concentration (0.00625%) in the culture vessel, mM
A_init          = 0;                        % Initial acetate concentration in the culture vessel, mM
TCN_init        = 2e8;                      % Initial number of total cells

% Growth of coculture w/o acetate supplementation
A_fm            = 0;                               % Feed medium acetate concentration, mM
CV101_if        = 0.5;                             % Initial fraction of CV101 cells
CV101_init      = TCN_init*CV101_if/V_env;         % initial population density of CV101, cell number/L
CV103_init      = TCN_init*(1-CV101_if)/V_env;     % initial population density of CV103, cell number/L
[t_exp3_wo_acetate,y_exp3_wo_acetate] = ode15s(...
    @acetate_mediated_cf_complete_model,  ... Model
    [0,150],         ... Time interval of integration (unit: h)
    [G_init,CV101_init,CV103_init,A_init], ... Initial condition
    options_ode15s, ... Matlab options
    D,G_fm,A_fm,parameters ... Other parameters
   );

% Growth of coculture w acetate supplementation
A_fm            = 1.0;                             % Feed medium acetate concentration, mM
CV101_if        = 0.63;                            % Initial fraction of CV101 cells
CV101_init      = TCN_init*CV101_if/V_env;         % initial population density of CV101, cell number/L
CV103_init      = TCN_init*(1-CV101_if)/V_env;     % initial population density of CV103, cell number/L
[t_exp3_w_acetate,y_exp3_w_acetate] = ode15s(...
    @acetate_mediated_cf_complete_model,  ... Model
    [0,150],         ... Time interval of integration (unit: h)
    [G_init,CV101_init,CV103_init,A_init], ... Initial condition
    options_ode15s, ... Matlab options
    D,G_fm,A_fm,parameters ... Other parameters
   );
    
y_obs = [y_obs;...
         Exp3Data_All(:,2).Variables;...
         Exp3Data_All(:,4).Variables];
y_sim = [y_sim;...
         pchip(t_exp3_wo_acetate/3.5, y_exp3_wo_acetate(:,3)./(y_exp3_wo_acetate(:,2)+y_exp3_wo_acetate(:,3)), Exp3Data_All(:,1).Variables);...
         pchip(t_exp3_w_acetate/3.5, y_exp3_w_acetate(:,3)./(y_exp3_w_acetate(:,2)+y_exp3_w_acetate(:,3)), Exp3Data_All(:,3).Variables)];
     
%% compute difference between observed and simulated data

% SMAPE: symmetric mean absolute percentage error
weights = [0.1*ones(ndata_EXP1,1);ones(length(y_obs)-ndata_EXP1,1)];
err = sum(weights.*(abs(y_obs-y_sim)./(abs(y_obs)+abs(y_sim))))/length(y_obs);

% MSE: mean squared error
% err = sum(weights.*((y_obs-y_sim).^2))/length(y_obs);

end

