% By Chen Liao, Memorial Sloan Kettering Cancer Center, June 02, 2020
%
% This file plots the growth rate ratio of CV101 to CV103 in nutrient space
% (spanned by glucose and acetate concentration)

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

%% Vary glucose and acetate
G_conc = 10.^[-3:0.01:3];  % External glucose concentration, mM
A_conc = 10.^[-3:0.01:3];  % Externaal acetate concentration, mM
growth_rate_CV101 = zeros(length(G_conc),length(A_conc));
growth_rate_CV103 = zeros(length(G_conc),length(A_conc));

% CV101 and CV103 have equal growth rate for certain glucose and acetate concentrations
isoline_acetate_index = zeros(length(A_conc),1);
isoline_growth_rate = zeros(length(A_conc),1);
 
for i=1:length(A_conc)
    for j=1:length(G_conc)
        [~,tmp] = acetate_mediated_cf_simplified_model(0, [G_conc(j),0,0,A_conc(i)], 0, 0, 0, parameters);
        growth_rate_CV101(j,i) = tmp(1);
        growth_rate_CV103(j,i) = tmp(2);
    end
    [~,isoline_acetate_index(i)] = min(abs(log10(growth_rate_CV101(:,i)./growth_rate_CV103(:,i))));
    if (isoline_acetate_index(i)==1)
        isoline_acetate_index(i) = nan;
        isoline_growth_rate(i) = nan;
    else
        isoline_growth_rate(i) = growth_rate_CV101(isoline_acetate_index(i),i);
    end
end

figure('Name','Fig. 2F');
hold on;
axis square;
box on;

% plot growth rate ratio
growth_rate_ratio = log10(growth_rate_CV101./growth_rate_CV103);
C = [[growth_rate_ratio(:,1:end-1) zeros(size(growth_rate_ratio(:,1:end-1),1),1)] ; zeros(1,size(growth_rate_ratio(:,1:end-1),2)+1)];
h=pcolor(C);
set(h, 'EdgeColor', 'none');

% plot isoline of equal growth rate
plot(1:length(A_conc),isoline_acetate_index,'k-','LineWidth',2);

% plot contour of growth rate = 0.2 for both CV101 and CV103
contour(growth_rate_CV101, [0,0.2], 'r--', 'ShowText', 'on', 'labelspacing', 700);
contour(growth_rate_CV103, [0,0.2], 'b--', 'ShowText', 'on', 'labelspacing', 700);

% draw intersections
[~,index] = sort(abs(isoline_growth_rate-0.2));
plot(index(1),isoline_acetate_index(index(1)),'ko','MarkerSize',12,'MarkerFaceColor',[0.5,0.5,0.5]);
plot(index(2),isoline_acetate_index(index(2)),'ko','MarkerSize',12,'MarkerFaceColor',[0.5,0.5,0.5]);

colormap(redblue(50));
c = colorbar;
caxis([-0.5,0.5]);
c.Location = 'southoutside';
set(gca,'XTick',[1:200:length(A_conc)]);
set(gca,'YTick',[1:200:length(G_conc)]);
set(gca,'Xticklabel',{'10^{-3}','10^{-1}','10^{1}','10^{3}'});
set(gca,'Yticklabel',{'10^{-3}','10^{-1}','10^{1}','10^{3}'});
xlabel('Acetate (mM)');
ylabel('Glucose (mM)');
axis([1,length(A_conc),1,length(G_conc)]);
