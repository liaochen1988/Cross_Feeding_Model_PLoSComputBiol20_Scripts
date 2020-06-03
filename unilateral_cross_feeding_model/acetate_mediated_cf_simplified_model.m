function [dxdt,J_grow] = acetate_mediated_cf_simplified_model(t, x, D, G_fm, A_fm, params)

% Input variables:
%       t:          time
%       x:          state variables
%       D:          dilution rate
%       G_fm:       feed medium glucose concentration
%       A_fm:       feed medium acetate concentration
%       params: 	cell parameters for CV101 and CV103

% Output variables:
%       dxdt:       derivatives
%       J_grow:     specific growth rates of CV101 and CV103

% State variables
Glu       = x(1); % glucose concentration in the culture vessel
CV101_a   = x(2); % active population density of CV101
CV103_a   = x(3); % active population density of CV103
Ace       = x(4); % acetate concentration in the culture vessel

% Unpack cell parameter values
% Suffix '1' and '3' represent CV101- and CV103-specific parameters respectively
% V_1_g:     maximum glucose uptake rate of CV101, mmol/h
% V_3_g:     maximum glucose uptake rate of CV103, mmol/h
% K_g:       Michaelis constant for glucose uptake, mM
% I_1_a:     half-inhibition constant for CV101 growth by acetate, mM
% I_3_a:     half-inhibition constant for CV103 growth by acetate, mM
% C_1_g:     half-inhibition constant for acetate uptake by glucose, mM
% V_1_a:     maximum acetate uptake rate of CV101, mmol/h
% K_1_a:     Michaelis constant for acetate uptake by CV101, mM
% gamma_a:   biomass yield of E. coli grown on acetate
% varphi_a:  acetate leakage fraction (percentage of carbon loss)
% delta_a:   number of acetate molecules produced per molecule of glucose is consumed

params_cell = num2cell(params);
[V_1_g, V_3_g, K_g, I_1_a, I_3_a, C_1_g, V_1_a, K_1_a, gamma_a, varphi_a, delta_a] = params_cell{:};

% Glucose uptake rate
J_upt_1_g = V_1_g*Glu/(K_g+Glu);
J_upt_3_g = V_3_g*Glu/(K_g+Glu);

% Acetate uptake rate is inhibited in presence of glucose
J_upt_1_a = V_1_a*Ace/(K_1_a+Ace)* C_1_g/(C_1_g+Glu);

% Specific growth rate is inhibited by acetate
J_grow_1 = gamma_a*(delta_a*(1-varphi_a)*J_upt_1_g+J_upt_1_a)*I_1_a/(I_1_a+Ace);
J_grow_3 = gamma_a*(delta_a*(1-varphi_a)*J_upt_3_g)*I_3_a/(I_3_a+Ace);
J_grow = [J_grow_1;J_grow_3];

% Differential equations
dxdt = zeros(4,1);
dxdt(1) = D*(G_fm-Glu)-J_upt_1_g*CV101_a-J_upt_3_g*CV103_a;  % glucose
dxdt(2) = (J_grow_1-D)*CV101_a;       % population density of CV101
dxdt(3) = (J_grow_3-D)*CV103_a;       % population density CV103
dxdt(4) = D*(A_fm-Ace)+...
          (delta_a*varphi_a*J_upt_1_g-J_upt_1_a)*CV101_a+...
          delta_a*varphi_a*J_upt_3_g*CV103_a;  % acetate
end

