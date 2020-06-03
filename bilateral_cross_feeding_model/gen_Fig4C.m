% By Chen Liao, Memorial Sloan Kettering Cancer Center, June 2, 2020
%
% This file plots isosurface of lysine and leucine auxotroph in the
% nutritional space

addpath('../auxiliary_functions');

%% Model parameters
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

%% vary feed medidum nutrient concentration
Glu_conc = 10.^[-1:0.01:2]; % mM
Lys_conc = 10.^[-3:0.01:0];  % External lysine concentration, mM
Leu_conc = 10.^[-4:0.01:0];  % External leucine concentration, mM

growth_rate_DeltaK = zeros(length(Lys_conc),length(Leu_conc),length(Glu_conc));
growth_rate_DeltaL = zeros(length(Lys_conc),length(Leu_conc),length(Glu_conc));
min_index_DeltaK = zeros(length(Lys_conc),length(Leu_conc),length(Glu_conc));
min_index_DeltaL = zeros(length(Lys_conc),length(Leu_conc),length(Glu_conc));

parfor k=1:length(Glu_conc)
    growth_rate_DeltaK_tmp = zeros(length(Lys_conc),length(Leu_conc));
    growth_rate_DeltaL_tmp = zeros(length(Lys_conc),length(Leu_conc));
    min_index_DeltaK_tmp = zeros(length(Lys_conc),length(Leu_conc));
    min_index_DeltaL_tmp = zeros(length(Lys_conc),length(Leu_conc));
    for j=1:length(Lys_conc)
        for i=1:length(Leu_conc)
            [~,tmp_1,tmp_2] = lys_leu_bilateral_cf_simplified_model(0, [Glu_conc(k),0,0,0,0,Lys_conc(j),Leu_conc(i)], 0, 0, 0, 0, parameters);
            growth_rate_DeltaK_tmp(j,i) = min(tmp_1)-parameters(end-1);
            growth_rate_DeltaL_tmp(j,i) = min(tmp_2)-parameters(end);
            if (growth_rate_DeltaK_tmp(j,i)<0)
                growth_rate_DeltaK_tmp(j,i) = NaN;
            end
            if (growth_rate_DeltaL_tmp(j,i)<0)
                growth_rate_DeltaL_tmp(j,i) = NaN;
            end 
            
            [~,min_index_DeltaK_tmp(j,i)] = min(tmp_1);
            [~,min_index_DeltaL_tmp(j,i)] = min(tmp_2);
            
        end
    end
    growth_rate_DeltaK(:,:,k) = growth_rate_DeltaK_tmp;
    growth_rate_DeltaL(:,:,k) = growth_rate_DeltaL_tmp;
    min_index_DeltaK(:,:,k) = min_index_DeltaK_tmp;
    min_index_DeltaL(:,:,k) = min_index_DeltaL_tmp;
end

%% get equal growht rate = 0.1
D_coor_1 = [];
D_coor_color_1 = [];

for k=1:length(Glu_conc)
    for j=1:length(Lys_conc)
        for i=1:length(Leu_conc)
            gr_dk = growth_rate_DeltaK(j,i,k);
            gr_dl = growth_rate_DeltaL(j,i,k);
            if (abs(gr_dk-0.1)<1e-3 && abs(gr_dl-0.1)<1e-3)
                D_coor_1 = [D_coor_1;j,i,k];
                min_index_dk = min_index_DeltaK(j,i,k);
                min_index_dl = min_index_DeltaL(j,i,k);
                if (min_index_dk==1 && min_index_dl==1)
                    D_coor_color_1 = [D_coor_color_1;1];
                end
                if (min_index_dk==1 && min_index_dl==2)
                    D_coor_color_1 = [D_coor_color_1;2];
                end
                if (min_index_dk==2 && min_index_dl==1)
                    D_coor_color_1 = [D_coor_color_1;3];
                end
                if (min_index_dk==2 && min_index_dl==2)
                    D_coor_color_1 = [D_coor_color_1;4];
                end
            end
        end
    end
end
            
%% plot isosurface
figure('Name','Fig4C');
hold on;
axis square;
box on;

growth_rate_ratio = log10(growth_rate_DeltaK./growth_rate_DeltaL);
isosurface(growth_rate_ratio,0);
alpha(0.5);
plot3(D_coor_1(:,2), D_coor_1(:,1), D_coor_1(:,3), 'k');
scatter3(D_coor_1(:,2), D_coor_1(:,1), D_coor_1(:,3),5,D_coor_color_1);

set(gca,'YTick',1:100:length(Lys_conc));
set(gca,'XTick',1:100:length(Leu_conc));
set(gca,'ZTick',1:100:length(Glu_conc));
set(gca,'Yticklabel',{'10^{-3}','10^{-2}','10^{-1}','10^{0}'});
set(gca,'Xticklabel',{'10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'});
set(gca,'Zticklabel',{'10^{-1}','10^{0}','10^{1}','10^{2}'});
xlabel('Leucine (mM)');
ylabel('Lysine (mM)');
zlabel('Glucose (mM)');
axis([1,length(Leu_conc),1,length(Lys_conc),1,length(Glu_conc)]);

view([-43,30]);