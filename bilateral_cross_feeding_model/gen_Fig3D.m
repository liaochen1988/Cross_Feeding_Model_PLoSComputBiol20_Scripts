% By Chen Liao, Memorial Sloan Kettering Cancer Center, May. 14, 2019

% This file plots phase diagram of obligate cross-feeding between two 
% amino acid auxotrophs by scanning medium glucose concentration and amino
% acid leakage fraction.

addpath('../auxiliary_functions');

%% Default Parameters
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

tol  = 1e-6;
options_ode15s  = odeset('RelTol', tol, 'AbsTol', tol * ones(7,1), 'NonNegative', [1:7]);
tspan           = [0,1e6];

%% scan parameters and simulate steady state
D               = 0.1;            % h, Dilution rate, 1/h
G_fm            = 10;    % mM, Feed medium glucose concentration
K_fm            = 0;    % mM, Feed medium lysine concentration
L_fm            = 0;    % mM, Feed medium leucine concentration
G_init          = 0;    % mM, Initial glucose concentration
K_init          = 0;    % mM, Initial lysine concentration
L_init          = 0;    % mM, Initial leucine concentration
TD_init      	= 1e9;  % Initial total cell density

phi_dk_l = 0:0.005:1;
phi_dl_k = 0:0.005:1;

DK_frac   = zeros(length(phi_dk_l),length(phi_dl_k));  % fraction of lysine auxotroph
eco_relat = zeros(length(phi_dk_l),length(phi_dl_k));  % ecological relationship

parfor i=1:length(phi_dk_l)
    DK_freq_tmp   = zeros(1,length(phi_dl_k));
    eco_relat_tmp = zeros(1,length(phi_dl_k));  
    for j=1:length(phi_dl_k)
        [i,j]
        
        curr_parameters = parameters;
        curr_parameters(9) = phi_dk_l(i);
        curr_parameters(10) = phi_dl_k(j);

        % only lysine auxotroph
        DK_freq_init  = 1.00;
        init_cond     = [G_init, TD_init*DK_freq_init, TD_init*(1-DK_freq_init),TD_init*DK_freq_init, TD_init*(1-DK_freq_init), K_init, L_init];
        [~,y_DK_alone] = ode15s(@lys_leu_bilateral_cf_simplified_modelv2, tspan, init_cond, options_ode15s, D, G_fm, K_fm, L_fm, curr_parameters);
        
        % only leucine auxotroph
        DK_freq_init     = 0.00;
        init_cond     = [G_init, TD_init*DK_freq_init, TD_init*(1-DK_freq_init),TD_init*DK_freq_init, TD_init*(1-DK_freq_init), K_init, L_init];
        [~,y_DL_alone] = ode15s(@lys_leu_bilateral_cf_simplified_modelv2, tspan, init_cond, options_ode15s, D, G_fm, K_fm, L_fm, curr_parameters);
        
        % lysine auxotroph : leucine auxotroph = 1:1
        DK_freq_init  = 0.50;
        init_cond     = [G_init, TD_init*DK_freq_init, TD_init*(1-DK_freq_init),TD_init*DK_freq_init, TD_init*(1-DK_freq_init), K_init, L_init];
        [~,y_coculture] = ode15s(@lys_leu_bilateral_cf_simplified_modelv2, tspan, init_cond, options_ode15s, D, G_fm, K_fm, L_fm, curr_parameters);
        
        % calculate frequency of lysine auxotroph
        DK_alone_ss = y_DK_alone(end,4);
        DL_alone_ss = y_DL_alone(end,5);
        DK_coculture_ss  = y_coculture(end,4);
        DL_coculture_ss  = y_coculture(end,5);
        if (DK_alone_ss <1)
            DK_alone_ss = 0;
        end
        if (DL_alone_ss < 1)
            DL_alone_ss = 0;
        end
        if (DK_coculture_ss < 1)
            DK_coculture_ss = 0;
        end
        if (DL_coculture_ss < 1)
            DL_coculture_ss = 0;
        end
        if (DK_coculture_ss==0 && DL_coculture_ss==0)
            DK_freq_tmp(j) = NaN;
        else
            DK_freq_tmp(j) = DK_coculture_ss/(DK_coculture_ss+DL_coculture_ss);
        end
        
        % get ecological relationship
        % 0: collapse
        % 1: competitative exclusion
        % 2: mutualism
        % 3: competition
        % 4; commensalism
        % 5: amensalism
        % 6: parasitism
        % 7: no effect
        if (DK_coculture_ss==0 && DL_coculture_ss==0)         
            eco_relat_tmp(j) = 0;
        else
            if (DK_coculture_ss ==0 || DL_coculture_ss==0)
                eco_relat_tmp(j) = 1;
            else
                sign_change_DK = sign(DK_coculture_ss-DK_alone_ss);
                sign_change_DL = sign(DL_coculture_ss-DL_alone_ss);
                if (sign_change_DK==1 && sign_change_DL==1)
                    eco_relat_tmp(j)=2;
                end
                if (sign_change_DK==-1 && sign_change_DL==-1)
                    eco_relat_tmp(j)=3;
                end
                if ((sign_change_DK==1 && sign_change_DL==0)||(sign_change_DK==0 && sign_change_DL==1))
                    eco_relat_tmp(j)=4;
                end
                if ((sign_change_DK==-1 && sign_change_DL==0)||(sign_change_DK==0 && sign_change_DL==-1))
                    eco_relat_tmp(j)=5;
                end
                if ((sign_change_DK==1 && sign_change_DL==-1)||(sign_change_DK==-1 && sign_change_DL==1))
                    eco_relat_tmp(j)=6;
                end
                if (sign_change_DK==0 && sign_change_DL==0)
                    eco_relat_tmp(j)=7;
                end
            end
        end
    end
    
    DK_frac(i,:)    = DK_freq_tmp;
    eco_relat(i,:)  = eco_relat_tmp;
end

save('output/Fig3D_phase_diagram.mat','parameters','D','G_fm','K_fm','L_fm','G_init','K_init','L_init','TD_init','phi_dk_l','phi_dl_k','DK_frac','eco_relat')

%% plot relative fraction of lysine auxotroph and ecological relationship
load('output/Fig3D_phase_diagram.mat');
figure('Name','Fig.3D');

% fraction of DK
ax1=subplot(1,2,1);
C = [[DK_frac(:,1:end-1) zeros(size(DK_frac(:,1:end-1),1),1)] ; zeros(1,size(DK_frac(:,1:end-1),2)+1)];
h=pcolor(C);
set(h, 'EdgeColor', 'none');
axis square;
box on;
hold on;
colormap(redblue(50));
c = colorbar;
caxis([0,1]);
axis([1,length(phi_dk_l),1,length(phi_dl_k)]);
c.Location = 'southoutside';
set(gca,'XTick',1:(length(phi_dk_l)-1)/5:length(phi_dk_l));
set(gca,'YTick',1:(length(phi_dl_k)-1)/5:length(phi_dl_k));
set(gca,'Xticklabel',{'0','0.2','0.4','0.6','0.8','1.0'});
set(gca,'Yticklabel',{'0','0.2','0.4','0.6','0.8','1.0'});
xlabel('Fraction of glucose leaked in form of leucine \phi_{\Deltak,l}');
ylabel('Fraction of glucose leaked in form of lysine \phi_{\Deltal,k}');
axis square;
box on;

% ecological relationship
ax2 = subplot(1,2,2);
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
axis([1,length(phi_dk_l),1,length(phi_dl_k)]);
set(gca,'XTick',1:(length(phi_dk_l)-1)/5:length(phi_dk_l));
set(gca,'XTick',1:(length(phi_dl_k)-1)/5:length(phi_dl_k));
set(gca,'Xticklabel',{'0','0.2','0.4','0.6','0.8','1.0'});
set(gca,'Yticklabel',{'0','0.2','0.4','0.6','0.8','1.0'});
xlabel('Fraction of glucose leaked in form of leucine \phi_{\Deltak,l}');
ylabel('Fraction of glucose leaked in form of lysine \phi_{\Deltal,k}');
axis square;
box on;
