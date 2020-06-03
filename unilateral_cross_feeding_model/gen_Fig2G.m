% By Chen Liao, Memorial Sloan Kettering Cancer Center, June 02, 2020
%
% This file plots phase diagram that spans on feed medium glucose
% concentration and acetate leakage fraction

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

%% Vary glucose supply level and acetate leakage fraction
D               = 0.2;              % h, Dilution rate
G_fm            = 10.^[-2:0.02:2];  % mM, Feed medium glucose concentration
A_fm            = 0;                % mM, Feed medium acetate concentration
G_init          = 0;                % mM, Initial glucose concentration
A_init          = 0;                % mM, Initial acetate concentration
TD_init         = 2e8;              % Initial total cell density
DensityPerA420  = 1e9;              % Cell density per OD420
varphi_a        = 0:0.005:1.0;                              % Acetate leakage fraction
CV101_frac      = zeros(length(G_fm),length(varphi_a));     % Fraction of CV101
eco_relat        = zeros(length(G_fm),length(varphi_a));     % Ecological relationship

parfor i=1:length(G_fm)
    CV101_freq_tmp = zeros(1,length(varphi_a));
    eco_relat_tmp  = zeros(1,length(varphi_a));
    for j=1:length(varphi_a)
        [i,j]
        curr_parameters = parameters;
        curr_parameters(10) = varphi_a(j);
        
        % CV101 alone
        CV101_frac_init = 1.00;
        init_cond   = [G_init, TD_init*CV101_frac_init, TD_init*(1-CV101_frac_init), A_init];
        [~,y_CV101_alone] = ode15s(@acetate_mediated_cf_simplified_model, tspan, init_cond, options_ode15s, D, G_fm(i), A_fm, curr_parameters);
        
        % CV103 alone
        CV101_frac_init = 0.00;
        init_cond   = [G_init, TD_init*CV101_frac_init, TD_init*(1-CV101_frac_init), A_init];
        [~,y_CV103_alone] = ode15s(@acetate_mediated_cf_simplified_model, tspan, init_cond, options_ode15s, D, G_fm(i), A_fm, curr_parameters);
        
        % CV101 : CV103 = 1:1
        CV101_frac_init     = 0.50;
        init_cond   = [G_init, TD_init*CV101_frac_init, TD_init*(1-CV101_frac_init), A_init];
        [~,y_coculture] = ode15s(@acetate_mediated_cf_simplified_model, tspan, init_cond, options_ode15s, D, G_fm(i), A_fm, curr_parameters);
        
        % calculate relative fraction of CV101
        CV101_alone_ss = y_CV101_alone(end,2);
        CV103_alone_ss = y_CV103_alone(end,3);
        CV101_coculture_ss  = y_coculture(end,2);
        CV103_coculture_ss  = y_coculture(end,3);
        if (CV101_alone_ss <1)
            CV101_alone_ss = 0;
        end
        if (CV103_alone_ss < 1)
            CV103_alone_ss = 0;
        end
        if (CV101_coculture_ss < 1)
            CV101_coculture_ss = 0;
        end
        if (CV103_coculture_ss < 1)
            CV103_coculture_ss = 0;
        end
        if (CV101_coculture_ss==0 && CV103_coculture_ss==0)
            CV101_freq_tmp(j) = NaN;
        else
            CV101_freq_tmp(j) = CV101_coculture_ss/(CV101_coculture_ss+CV103_coculture_ss);
        end
        
        % infer ecological relationship
        % 0: collapse
        % 1: competitative exclusion
        % 2: mutualism
        % 3: competition
        % 4; commensalism
        % 5: amensalism
        % 6: parasitism
        % 7: no effect
        if (CV101_coculture_ss==0 && CV103_coculture_ss==0)
            eco_relat_tmp(j) = 0;
        else
            if (CV101_coculture_ss ==0 || CV103_coculture_ss==0)
                eco_relat_tmp(j) = 1;
            else    
                sign_change_CV101 = sign(CV101_coculture_ss-CV101_alone_ss);
                sign_change_CV103 = sign(CV103_coculture_ss-CV103_alone_ss);
                if (sign_change_CV101==1 && sign_change_CV103==1)
                    eco_relat_tmp(j)=2;
                end
                if (sign_change_CV101==-1 && sign_change_CV103==-1)
                    eco_relat_tmp(j)=3;
                end
                if ((sign_change_CV101==1 && sign_change_CV103==0)||(sign_change_CV101==0 && sign_change_CV103==1))
                    eco_relat_tmp(j)=4;
                end
                if ((sign_change_CV101==-1 && sign_change_CV103==0)||(sign_change_CV101==0 && sign_change_CV103==-1))
                    eco_relat_tmp(j)=5;
                end
                
                if ((sign_change_CV101==1 && sign_change_CV103==-1)||(sign_change_CV101==-1 && sign_change_CV103==1))
                    eco_relat_tmp(j)=6;
                end
                
                if (sign_change_CV101==0 && sign_change_CV103==0)
                    eco_relat_tmp(j)=7;
                end
            end
        end
    end
    
    CV101_frac(i,:) = CV101_freq_tmp;
    eco_relat(i,:)  = eco_relat_tmp;
end

save('output/Fig2G_phase_diagram.mat','parameters','D','G_fm','A_fm','G_init','A_init','DensityPerA420','varphi_a','CV101_frac','eco_relat')

%% plot relative fraction of CV101 and ecological relationship
load(sprintf('output/fig2G_phase_diagram.mat'));
figure('Name','Fig.2G');

ax1=subplot(1,2,1);
hold on;
C = [[CV101_frac(:,1:end-1) zeros(size(CV101_frac(:,1:end-1),1),1)] ; zeros(1,size(CV101_frac(:,1:end-1),2)+1)];
h=pcolor(C);
set(h, 'EdgeColor', 'none');
axis square;
box on;
colormap(ax1,redblue(50));
c = colorbar;
caxis([0,1]);
c.Location = 'southoutside';
axis([1,length(varphi_a)*0.5,1,length(G_fm)]);
set(gca,'XTick',1:(length(varphi_a)*0.5-1)/5:length(varphi_a)*0.5);
set(gca,'YTick',1:(length(G_fm)-1)/4:length(G_fm));
set(gca,'Xticklabel',{'0','0.1','0.2','0.3','0.4','0.5'});
set(gca,'Yticklabel',{'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}'});
xlabel('Fraction of acetate leakage');
ylabel('Feed medium glucose (mM)');

% theoretical model
% figure();
% hold on;
% varphi_a_finer = [0:0.0001:1];
% params_cell = num2cell(parameters);
% [V_1_g, V_3_g, K_g, I_1_a, I_3_a, C_1_g, V_1_a, K_1_a, gamma_a, ~, delta_a] = params_cell{:};
% DV_g = (V_3_g-V_1_g)/V_3_g;
% lower_boundary = D*(K_g./V_3_g./gamma_a./delta_a./(1-varphi_a_finer)+K_1_a*DV_g./(V_1_a*gamma_a*delta_a*varphi_a_finer));
% upper_boundary = D*(K_g./V_3_g./gamma_a./delta_a./(1-varphi_a_finer)+K_1_a*DV_g*(1-DV_g)./(V_1_a*gamma_a*delta_a*(varphi_a_finer-DV_g)));
% plot(varphi_a_finer,lower_boundary,'k-');
% plot(varphi_a_finer,upper_boundary,'k-');
% axis([0,0.5,1e-2,1e2]);
% set(gca,'XTick',[0:0.1:0.5]);
% set(gca,'YTick',10.^[-2:1:2]);
% set(gca,'YScale','log');
% xlabel('Fraction of acetate leakage');
% ylabel('Feed medium glucose (mM)');
% axis square;
% box on;
% 

ax2=subplot(1,2,2);
C = [[eco_relat(:,1:end-1) zeros(size(eco_relat(:,1:end-1),1),1)] ; zeros(1,size(eco_relat(:,1:end-1),2)+1)];
h=pcolor(C);
set(h, 'EdgeColor', 'none');
axis square;
box on;
colormap(ax2,distinguishable_colors(8));
c = colorbar;
caxis([-0.5,7.5]);
c.Ticks = linspace(0, 7, 8); 
c.Location = 'southoutside';
axis([1,length(varphi_a)*0.5,1,length(G_fm)]);
set(gca,'XTick',1:(length(varphi_a)*0.5-1)/5:length(varphi_a)*0.5);
set(gca,'YTick',1:(length(G_fm)-1)/4:length(G_fm));
set(gca,'Xticklabel',{'0','0.1','0.2','0.3','0.4','0.5'});
set(gca,'Yticklabel',{'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}'});
xlabel('Fraction of acetate leakage');
ylabel('Feed medium glucose (mM)');
