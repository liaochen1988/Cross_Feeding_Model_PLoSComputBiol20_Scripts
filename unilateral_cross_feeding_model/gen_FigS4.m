% By Chen Liao, Memorial Sloan Kettering Cancer Center, May 17, 2020
%
% This file plots phase diagram of Fig.2G by varing C_1_g and I_3_a
 
addpath('../auxiliary_functions');

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

tol  = 1e-6;
options_ode15s  = odeset('RelTol', tol, 'AbsTol', tol * ones(4,1), 'NonNegative', 1:4);
tspan           = [0,1e6];

D               = 0.2;              % h, Dilution rate
G_fm            = 10.^[-2:0.05:2];  % mM, Feed medium glucose concentration
A_fm            = 0;                % mM, Feed medium acetate concentration
G_init          = 0;                % mM, Initial glucose concentration
A_init          = 0;                % mM, Initial acetate concentration
TD_init         = 2e8;              % Initial total cell density
DensityPerA420  = 1e9;              % Density per OD
varphi_a        = 0:0.01:0.5;       % Acetate leakage fraction (percentage of carbon loss)
C_1_g           = 10.^[-3:0.5:1];   % Half-inhibition constant for acetate uptake by glucose, mM
I_3_a           = 10.^[-1:0.5:3];    % Half-inhibition constant for CV103 growth by acetate, mM

coexistence_region = zeros(length(C_1_g), length(I_3_a), length(G_fm), length(varphi_a));
for m=1:length(C_1_g)
    for n=1:length(I_3_a)
        [m,n]
        curr_coexistence_region  = zeros(length(G_fm),length(varphi_a));
        parfor i=1:length(G_fm)
            CV101_freq_tmp = zeros(1,length(varphi_a));
            is_coexistence  = zeros(1,length(varphi_a));
            for j=1:length(varphi_a)
                curr_parameters = parameters;
                curr_parameters(5) = I_3_a(n);
                curr_parameters(6) = C_1_g(m);
                curr_parameters(10) = varphi_a(j);
                
                % CV101 : CV103 = 1:1
                CV101_frac_init     = 0.50;
                init_cond   = [G_init, TD_init*CV101_frac_init, TD_init*(1-CV101_frac_init), A_init];
                [~,y_coculture] = ode15s(@acetate_mediated_cf_simplified_model, tspan, init_cond, options_ode15s, D, G_fm(i), A_fm, curr_parameters);
                
                % calculate relative fraction of CV101
                CV101_coculture_ss  = y_coculture(end,2);
                CV103_coculture_ss  = y_coculture(end,3);
                if (CV101_coculture_ss < 1)
                    CV101_coculture_ss = 0;
                end
                if (CV103_coculture_ss < 1)
                    CV103_coculture_ss = 0;
                end
                
                if (CV101_coculture_ss>0 && CV103_coculture_ss>0)
                    is_coexistence(j) = 1;
                else
                    is_coexistence(j) = NaN;
                end
            end
            curr_coexistence_region(i,:)  = is_coexistence;
        end
        
        coexistence_region(m,n,:,:) = curr_coexistence_region;
    end
end

save('output/figS4_phase_diagram_sensitivity.mat');

%% plot phase diagram
load('output/figS4_phase_diagram_sensitivity.mat');
figure();
for m=1:length(C_1_g)
    for n=1:length(I_3_a)
        subplot(length(C_1_g),length(I_3_a),(m-1)*length(I_3_a)+n);
        curr_coexistence_region = squeeze(coexistence_region(m,n,:,:));
        C = [[curr_coexistence_region(:,1:end-1) zeros(size(curr_coexistence_region(:,1:end-1),1),1)] ; zeros(1,size(curr_coexistence_region(:,1:end-1),2)+1)];
        h=pcolor(C);
        colormap('hsv');
        set(h, 'EdgeColor', 'none');
        axis square;
        box on;
        set(gca,'XTicklabel','');
        set(gca,'YTicklabel','');
    end
end