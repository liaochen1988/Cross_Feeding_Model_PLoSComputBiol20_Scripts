function err = residual_function(params,data)

%% unpack parameters
parameters = [...
    3.61e-12; ... Maximum glucose uptake rate, mmol/h
    params(1);... Michaelis constant for glucose uptake, mM
    params(2);... Maximum lysine uptake rate by the lysine auxotroph, mmol/h
    params(3);... Maximum leucine uptake rate by the leucine axutroph, mmol/h
    5.00e-3;  ... Michaelis constant for lysine uptake, mM
    1.07e-3;  ... Michaelis constant for leucine uptake, mM
    1.00;     ... Number of lysine molecules produced per molecule of glucose is consumed
    1.00;     ... Number of leucine molecules produced per molecule of glucose is consumed
    params(4);... Fraction of glucose leaked in form of leucine by the lysine auxotroph
    params(5);... Fraction of glucose leaked in form of lysine by the leucine auxotroph
    3.00e11;  ... Biomass yield of E. coli grown on glucose, 1/mmol
    params(6);... Biomass yield of E. coli grown on lysine, 1/mmol
    params(7);... Biomass yield of E. coli grown on leucine, 1/mmol
    params(8);... Mortality rate constant of the lysine auxotroph, 1/h
    params(9) ... Mortality rate constant of the leusine auxotroph, 1/h
    ];

%% Get data
Exp1Data_DeltaK = data.Exp1Data_DeltaK;
Exp2Data_DeltaL = data.Exp2Data_DeltaL;
Exp3Data_DeltaKL = data.Exp3Data_DeltaKL;

%% vectorized difference between simulated and experimental data 
y_obs = [];
y_sim = [];
Density_per_OD600 = 5e8.*[1.6,1]; % linear coefficient between cell density (number of cells /mL) of [lysine,leucine] auxotroph and OD600
options_ode15s    = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(7,1),'NonNegative',[1:7]); % Matlab options
         
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

[tK,yK] = ode15s(@lys_leu_bilateral_cf_simplified_model,... Model
                 [0,15],                                ... Time interval of integration (unit: h)
                 [G_init,DK_init,0,DK_init,0,K_init,0], ... Initial condition
                 options_ode15s,                        ... Matlab options
                 D,G_fm,K_fm,L_fm,parameters);

y_obs = [y_obs;...
         remove_NaN(Exp1Data_DeltaK(:,2).Variables);...
         remove_NaN(Exp1Data_DeltaK(:,4).Variables);...
         remove_NaN(Exp1Data_DeltaK(:,6).Variables)];
y_sim = [y_sim;...
         pchip(tK, yK(:,2)/Density_per_OD600(1)/1e3, remove_NaN(Exp1Data_DeltaK(:,1).Variables));...
         pchip(tK, yK(:,1), remove_NaN(Exp1Data_DeltaK(:,3).Variables));...
         pchip(tK, yK(:,6)*K_mm, remove_NaN(Exp1Data_DeltaK(:,5).Variables))];

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

[tL,yL] = ode15s(@lys_leu_bilateral_cf_simplified_model,... Model
                 [0,15],                                ... Time interval of integration (unit: h)
                 [G_init,0,DL_init,0,DL_init,0,L_init], ... Initial condition
                 options_ode15s,                        ... Matlab options
                 D,G_fm,K_fm,L_fm,parameters);

y_obs = [y_obs;...
         remove_NaN(Exp2Data_DeltaL(:,2).Variables);...
         remove_NaN(Exp2Data_DeltaL(:,4).Variables);...
         remove_NaN(Exp2Data_DeltaL(:,6).Variables)];
y_sim = [y_sim;...
         pchip(tL, yL(:,3)/Density_per_OD600(2)/1e3, remove_NaN(Exp2Data_DeltaL(:,1).Variables));...
         pchip(tL, yL(:,1), remove_NaN(Exp2Data_DeltaL(:,3).Variables));...
         pchip(tL, yL(:,7)*L_mm, remove_NaN(Exp2Data_DeltaL(:,5).Variables))];

%% Experiment 3: Growth of the lysine and leucine auxotrophs in glucose-minimum batch coculture
G_init  = Exp3Data_DeltaKL(1,4).Variables;                   % Initial glucose concentration, mM
[tKL,yKL] = ode15s(@lys_leu_bilateral_cf_simplified_model,           ... Model
                   [0,80],                                ... Time interval of integration (unit: h)
                   [G_init,DK_init,DL_init,DK_init,DL_init,0,0], ... Initial condition
                   options_ode15s,                        ... Matlab options
                   D,G_fm,K_fm,L_fm,parameters);
             
y_obs = [y_obs;...
         remove_NaN(Exp3Data_DeltaKL(:,2).Variables);...
         remove_NaN(Exp3Data_DeltaKL(:,4).Variables);...
         remove_NaN(Exp3Data_DeltaKL(:,6).Variables)];
y_sim = [y_sim;...
         pchip(tKL, yKL(:,2)/Density_per_OD600(1)/1e3+yKL(:,3)/Density_per_OD600(2)/1e3, remove_NaN(Exp3Data_DeltaKL(:,1).Variables));...
         pchip(tKL, yKL(:,1), remove_NaN(Exp3Data_DeltaKL(:,3).Variables));...
         pchip(tKL, yKL(:,4)./yKL(:,5), remove_NaN(Exp3Data_DeltaKL(:,5).Variables))];
                  
%% compute difference between observed and simulated data

% SMAPE: symmetric mean absolute percentage error
err = sum(abs(y_obs-y_sim)./(abs(y_obs)+abs(y_sim)))/length(y_obs);

% MSE: mean squared error
% err = sum((y_obs-y_sim).^2)/length(y_obs);

end

