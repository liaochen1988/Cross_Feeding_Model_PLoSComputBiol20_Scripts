function [dxdt,J_grow_dk_arr,J_grow_dl_arr] = lys_leu_bilateral_cf_simplified_model(t, x, D, G_fm, K_fm, L_fm, params)

% Input variables:
%       t:          time
%       x:          state variables
%       D:          dilution rate
%       G_fm:       feed medium glucose concentration
%       K_fm:       feed medium lysine concentration
%       L_fm:       feed medium leucine concentration
%       Cell_p: 	cell parameters for CV101 and CV103

% Output variables:
%       dxdt:       derivatives
%       J_grow:     specific growth rates of lysine and leucine auxotroph

% State variables
Glu = x(1);       % glucose concentration in the culture vessel
DeltaK_a = x(4);  % active populatiton density of lysine auxotroph
DeltaL_a = x(5);  % active populatiton density of leucine auxotroph
Lys = x(6);       % lysine conconcentration in the culture vessel
Leu = x(7);       % leucine concentration in the culture vessel

% Unpack cell parameter values
% Suffix 'dk' and 'dl' represent lysine auxotroph- and leucine auxotroph-specific parameters respectively
% Suffix 'k' and 'l' represent lysine- and leucine-specific parameters respectively
% Vmax_g:                       Maximum glucose uptake rate, mmol/h
% Km_g:                         Michaelis constant for glucose uptake, mM
% Vmax_dk_k:                    Maximum lysine uptake rate by the lysine auxotroph, mmol/h
% Vmax_dl_l:                    Maximum leucine uptake rate by the leucine axutroph, mmol/h
% Km_k:                         Michaelis constant for lysine uptake, mM
% Km_l:                         Michaelis constant for leucine uptake, mM
% delta_k:                      Number of lysine molecules produced per molecule of glucose is consumed
% delta_l:                      Number of leucine molecules produced per molecule of glucose is consumed
% phi_dl_k:                     Fraction of glucose uptake leaked in form of lysine production
% phi_dk_l:                     Fraction of glucose uptake leaked in form of leucine production
% gamma_g:                      Biomass yield of E. coli grown on glucose
% gamma_k:                      Biomass yield of E. coli grown on lysine
% gamma_l:                      Biomass yield of E. coli grown on leucine
% eta_dk:                       Mortality rate constant of the lysine auxotroph, 1/h
% eta_dl:                       Mortality rate constant of the leusine auxotroph, 1/h
 
params_cell = num2cell(params(:,1));
[Vmax_g,Km_g,Vmax_dk_k,Vmax_dl_l,Km_k,Km_l,delta_k,delta_l,...
 phi_dk_l,phi_dl_k,gamma_g,gamma_k,gamma_l,eta_dk,eta_dl] = params_cell{:};
         
% Glucose uptake rate
J_upt_dk_g = Vmax_g*Glu/(Km_g+Glu);
J_upt_dl_g = J_upt_dk_g;

% Leakage flux
J_leak_dk_l = J_upt_dk_g * delta_l * phi_dk_l;
J_leak_dl_k = J_upt_dl_g * delta_k * phi_dl_k;

% Amino acid uptake rate
J_upt_dk_k = Vmax_dk_k*Lys/(Km_k+Lys);
J_upt_dl_l = Vmax_dl_l*Leu/(Km_l+Leu);

% Specific growth rate
J_grow_dk_arr = [gamma_k*J_upt_dk_k, gamma_g*(1-phi_dk_l)*J_upt_dk_g];
J_grow_dk = min(J_grow_dk_arr);
J_grow_dl_arr = [gamma_l*J_upt_dl_l, gamma_g*(1-phi_dl_k)*J_upt_dl_g];
J_grow_dl = min(J_grow_dl_arr);

% Differential equations
dxdt = zeros(7,1);
dxdt(1) = D*(G_fm-Glu)-J_upt_dk_g*DeltaK_a-J_upt_dl_g*DeltaL_a; % glucose
dxdt(2) = DeltaK_a*(J_grow_dk-D); % total population density of lysine auxotroph
dxdt(3) = DeltaL_a*(J_grow_dl-D); % total population density of leucine auxotroph
dxdt(4) = DeltaK_a*(J_grow_dk-eta_dk-D); % active population density of lysine auxotroph
dxdt(5) = DeltaL_a*(J_grow_dl-eta_dl-D); % active population density of leucine auxotroph
dxdt(6) = D*(K_fm-Lys)+(J_leak_dl_k*DeltaL_a-J_upt_dk_k*DeltaK_a); % lysine
dxdt(7) = D*(L_fm-Leu)+(J_leak_dk_l*DeltaK_a-J_upt_dl_l*DeltaL_a); % leucine

end