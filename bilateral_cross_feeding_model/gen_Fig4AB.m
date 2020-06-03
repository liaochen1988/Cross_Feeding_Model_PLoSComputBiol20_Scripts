% By Chen Liao, Memorial Sloan Kettering Cancer Center, May. 14, 2019

% This file plots phase diagram of obligate cross-feeding between two 
% amino acid auxotrophs by scanning medium lysine and leucine concentration

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
 
%% vary medium concentration of lysine and leucine
D             = 0.1;    % h, Dilution rate, 1/h
G_fm          = 10;     % mM, Feed medium glucose concentration
G_init        = 0;      % mM, Initial glucose concentration
K_init        = 0;      % mM, Initial lysine concentration
L_init        = 0;      % mM, Initial leucine concentration
TD_init       = 1e9;    % Initial total cell density

K_fm          = 10.^[-4:0.02:1];   % mM, Feed medium lysine concentration
L_fm          = 10.^[-4:0.02:1];   % mM, Feed medium leucine concentration
DK_frac       = zeros(length(K_fm),length(L_fm)); % fraction of lysine auxotroph
eco_relat     = zeros(length(K_fm),length(L_fm)); % ecological relationship

parfor i=1:length(K_fm)
    DK_freq_tmp = zeros(1,length(L_fm));
    eco_relat_tmp = zeros(1,length(L_fm));
    for j=1:length(L_fm)
        [i,j]
        
        % only lysine auxotroph
        DK_freq_init  = 1.00;
        init_cond     = [G_init, TD_init*DK_freq_init, TD_init*(1-DK_freq_init),TD_init*DK_freq_init, TD_init*(1-DK_freq_init), K_init, L_init];
        [~,y_DK_alone] = ode15s(@lys_leu_bilateral_cf_model, tspan, init_cond, options_ode15s, D, G_fm, K_fm(i), L_fm(j), parameters);
        
        % only leucine auxotroph
        DK_freq_init     = 0.00;
        init_cond     = [G_init, TD_init*DK_freq_init, TD_init*(1-DK_freq_init),TD_init*DK_freq_init, TD_init*(1-DK_freq_init), K_init, L_init];
        [~,y_DL_alone] = ode15s(@lys_leu_bilateral_cf_model, tspan, init_cond, options_ode15s, D, G_fm, K_fm(i), L_fm(j), parameters);
        
        % lysine auxotroph : leucine auxotroph = 1:1
        DK_freq_init     = 0.50;
        init_cond     = [G_init, TD_init*DK_freq_init, TD_init*(1-DK_freq_init),TD_init*DK_freq_init, TD_init*(1-DK_freq_init), K_init, L_init];
        [~,y_coculture] = ode15s(@lys_leu_bilateral_cf_model, tspan, init_cond, options_ode15s, D, G_fm, K_fm(i), L_fm(j), parameters);
        
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
        
        % infer ecological relationship
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

save('output/Fig4A_phase_diagram.mat','parameters','D','G_fm','K_fm','L_fm','G_init','K_init','L_init','TD_init','DK_frac','eco_relat')

%% plot relative fraction of lysine auxotroph and ecological relationship
load('output/fig4A_phase_diagram.mat');
figure('Name','Fig.4A');

ax1=subplot(1,2,1);
C = [[DK_frac(:,1:end-1) zeros(size(DK_frac(:,1:end-1),1),1)] ; zeros(1,size(DK_frac(:,1:end-1),2)+1)];
h=pcolor(C);
set(h, 'EdgeColor', 'none'); 
colormap(redblue(50));
c = colorbar;
caxis([0,1]);
c.Location = 'southoutside';
axis([1,length(L_fm),1,length(K_fm)]);
set(gca,'XTick',1:(length(L_fm)-1)/2:length(L_fm));
set(gca,'YTick',1:(length(K_fm)-1)/2:length(K_fm));
set(gca,'Xticklabel',{'10^{-4}','10^{-1.5}','10^{1}'});
set(gca,'Yticklabel',{'10^{-4}','10^{-1.5}','10^{1}'});
xlabel('Feed medium leusine (mM)');
ylabel('Feed medium lysine (mM)');
axis square;
box on;

ax2=subplot(1,2,2);
C = [[eco_relat(:,1:end-1) zeros(size(eco_relat(:,1:end-1),1),1)] ; zeros(1,size(eco_relat(:,1:end-1),2)+1)];
h=pcolor(C);
set(h, 'EdgeColor', 'none');
colormap(ax2,distinguishable_colors(8));
c = colorbar;
caxis([-0.5,7.5]);
c.Ticks = linspace(0, 7, 8); 
c.Location = 'southoutside';
axis([1,length(L_fm),1,length(K_fm)]);
set(gca,'XTick',1:(length(L_fm)-1)/2:length(L_fm));
set(gca,'YTick',1:(length(K_fm)-1)/2:length(K_fm));
set(gca,'Xticklabel',{'10^{-4}','10^{-1.5}','10^{1}'});
set(gca,'Yticklabel',{'10^{-4}','10^{-1.5}','10^{1}'});
xlabel('Feed medium leusine (mM)');
ylabel('Feed medium lysine (mM)');
axis square;
box on;

figure('Name','Fig.4B');
K_fm_rep  = [0.01, 0.01, 10];   % mM, Feed medium lysine concentration
L_fm_rep  = [0.01, 10, 0.01];   % mM, Feed medium leucine concentration
for i=1:length(K_fm_rep)
    DK_frac_init  = 0.50;
    init_cond     = [G_init, TD_init*DK_frac_init, TD_init*(1-DK_frac_init),TD_init*DK_frac_init, TD_init*(1-DK_frac_init), K_init, L_init];
    [t,y_coculture] = ode15s(@lys_leu_bilateral_cf_simplified_model, tspan, init_cond, options_ode15s, D, G_fm, K_fm_rep(i), L_fm_rep(i), parameters);
                   
    % get growth rate determined by different factors
    J_grow_dk_arr = zeros(length(t),2);
    J_grow_dl_arr = zeros(length(t),2);
    for k=1:length(t)
        [~,J_grow_dk_arr(k,:),J_grow_dl_arr(k,:)] = lys_leu_bilateral_cf_simplified_model(t(k), y_coculture(k,:), D, G_fm, K_fm_rep(i), L_fm_rep(i), parameters);
    end
          
    % Nutrient
    subplot(length(K_fm_rep),3,i);
    hold on;
    plot(t/(log(2)/D), y_coculture(:,1), 'k-'); % glucose
    plot(t/(log(2)/D), y_coculture(:,6), 'r-'); % lysine
    plot(t/(log(2)/D), y_coculture(:,7), 'b-'); % leucine
    legend('glucose','lysine','leucine');
    xlim([0,15]);
    set(gca,'XTick',[0,5,10,15]);
    set(gca,'XTickLabels',{});
    set(gca,'XminorTick','off')
    ylabel('Nutrient (mM)');
    box on;
    set(gca,'ticklength',3*get(gca,'ticklength'));
    ylim([0.0005,20]);
    set(gca,'YScale','log');
    set(gca,'YTick',[1e-3,1e-1,1e1]);
    set(gca,'YTickLabels',{'10^{-3}';'10^{-1}';'10^{1}'});
    set(gca,'YminorTick','off')
    
    % log10(Cell density)
    subplot(length(K_fm_rep),3,3+i);
    hold on;
    plot(t/(log(2)/D), log10(y_coculture(:,4)/1e3), 'r-'); % lysine auxotroph
    plot(t/(log(2)/D), log10(y_coculture(:,5)/1e3), 'b-'); % leucine auxotroph
    legend('\DeltaK','\DeltaL');
    xlim([0,15]);
    set(gca,'XTick',[0,5,10,15]);
    set(gca,'XTickLabels',{});
    ylabel('log_{10} (Cell density) (cells/mL)');
    box on;
    set(gca,'ticklength',3*get(gca,'ticklength'));
    ylim([5,10]);
    set(gca,'YTick',[5,7.5,10]);
    set(gca,'YTickLabels',{'5';'7.5';'10'});
    
    % growth rate difference when limited by glucose and limited by 
    % auxotrophic amino acids
    subplot(length(K_fm_rep),3,6+i);
    hold on;
    plot(t/(log(2)/D), J_grow_dk_arr(:,1)-J_grow_dk_arr(:,2),'r-');
    plot(t/(log(2)/D), J_grow_dl_arr(:,1)-J_grow_dl_arr(:,2),'b-');
    legend('\DeltaK','\DeltaL');
    xlim([0,15]);
    set(gca,'XTick',[0,5,10,15]);
    set(gca,'XTickLabels',{'0';'5';'10';'15'});
    xlabel('Generation');
    ylabel('\DeltaGR (h^{-1})');
    box on;
    set(gca,'ticklength',3*get(gca,'ticklength'));
    ylim([-1,1]);
    plot([0,15],[0,0],'k--');
 end