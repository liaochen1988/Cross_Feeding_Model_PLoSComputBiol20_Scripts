function [dxdt, J_grow] = multi_aa_auxotroph_cf_model_w_cheater(t, x, Vmax_g, Km_g, Vmax_aa, Km_aa, Yield_g, Yield_aa, Rcarbon, Byp_frac, Cdr)

% Input variables:
%       t:          time
%       x:          state variables
%       Vmax_g:     maximum glucose uptake rates, umol/h
%       Km_g:       michaelis constants for glucose uptake, uM
%       Vmax_aa:    maximum amino acid uptake rates, umol/h
%       Km_aa:      michaelis constants for amino acid uptake, uM
%       Yield_g:    biomass yields of E. coli auxotrophs growing on glucose, 1/umol
%       Yield_aa:   biomass yields of E. coli auxotrophs growing on auxotrophic amino acids, 1/umol
%       Rcarbon:    number of amino acid molecules produced per molecule of glucose is consumed
%       Byp_frac:   amino acid leakge fraction (percentage of carbon loss)
%       Cdr:        cell death rate, 1/h

% Output variables:
%       dxdt:       derivatives
%       J_grow:     specific growth rates of all auxotrophs

Glucose             = x(1);         % glucose concentration in the culture vessel
Aux_active          = x(2:15);      % population density of active cells of amino acid auxotrophs
Aux_cheater_active  = x(16:29);     % population density of active cheaters
AA                  = x(30:43);     % amino acid concentration in the culture vessel

% Glucose uptake rate, umol/h
Jg_upt  = Vmax_g*Glucose./(Km_g+Glucose);

% Amino acid uptake rate, umol/h
% We assume that each cell only uptakes the amino acid it cannot produce
Ja_upt  = Vmax_aa.*AA./(Km_aa+AA);

% Specific growth rate, 1/h
J_grow  = zeros(14,1);
for i=1:14
    J_grow(i) = min(Yield_aa(i) * Ja_upt(i), Yield_g(i) * Jg_upt(i) * (1-sum(Byp_frac(:,i))));
end

J_grow_cheater  = zeros(14,1);
for i=1:14
    J_grow_cheater(i) = min(Yield_aa(i) * Ja_upt(i), Yield_g(i) * Jg_upt(i));
end

% Differential equations
dxdt            = zeros(1+5*14,1);
dxdt(1)         = -sum(Jg_upt.*(Aux_active+Aux_cheater_active));  % glucose
dxdt(2:15,1)    = Aux_active.*(J_grow-Cdr);                       % active population of amino acid auxotrophs
dxdt(16:29,1)   = Aux_cheater_active.*(J_grow_cheater-Cdr);       % active population of amino acid auxotroph cheaters
for k=1:14
    dxdt(k+29)  = Rcarbon(k)*sum(Jg_upt.*Byp_frac(k,:)'.*Aux_active)-(Aux_active(k)+Aux_cheater_active(k))*Ja_upt(k); % amino acids
end
dxdt(44:57,1)   = Aux_active.*J_grow;                   % total population of amino acid auxotrophs
dxdt(58:71,1)   = Aux_cheater_active.*J_grow_cheater;   % total population of amino acid auxotroph cheaters

end

