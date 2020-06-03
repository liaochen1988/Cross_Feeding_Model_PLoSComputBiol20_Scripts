% This script tests how variation of external conditions
% (1) substrate concentraiton
% (2) antibiotic concentration
% affect the ratio of metabolite overflow flux to the total metabolite influx
%
% Last updated by Chen Liao, June 02, 2020

%% model parameters
alpha = 3e6;        % uM, total amino acid concentration
k_e1 = 1.2e3;       % 1/h, maximum metabolite biosynthesis rate
m_e1 = 325;         % number of amino acids in a metabolic protein
m_e2 = 7336*1.6;    % number of amino acids in a ribosomal protein
k_e2 = 7.56e4;      % 1/h, maximum rate of peptide elongation
Km_s = 100;         % uM, Michaelis constant
Km_m = 20;          % uM, Michaelis constant
Ki_r = 60;          % uM, half-inhibition constant
k_r = 3.6e3;        % 1/h, biosynthesis rate of transcriptional regulator
d_r = 1.26e2;       % 1/h, degradation rate constant of transcriptional regulator
Ym = 1/alpha;       % 1/uM, yield of amino acid
K_ia = 1;           % uM, half-inhibition constant of ribosome synthesis by antibiotic

k_m = [0,1e4,5e4,1e5,5e5,1e6,5e6,1e7]; % 1/h, diffusion rate constant
my_colors = jet(length(k_m)); % line color for each k_m value
tspan = [0,1e6]; % h, time span
options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6 1e-6],'NonNegative',[1,2,3]);

figure('Name','Fig. S11');

%% vary external substrate concentration
substrate = 10.^[-2:0.2:3]; % uM, concentration of external substrate
antibiotic = 0; % uM, concentration of external antibiotic
varphi = zeros(length(k_e1), length(k_m)); % fraction of metabolite influx that leak into environment
Mhat = zeros(length(k_e1), length(k_m)); % metabolite concentration
J_grow = zeros(length(k_e1), length(k_m)); % specific growth rate
J_in = zeros(length(k_e1), length(k_m)); % influx of metabolite
for i=1:length(substrate)
    y0 = [0,0.5,0];
    for j=1:length(k_m)
        [t,y] = ode15s(@ecoli_growth_model,tspan,y0,options,substrate(i),k_m(j),m_e1,m_e2,k_e1,k_e2,Km_s,Km_m,Ki_r,k_r,d_r,alpha,Ym,K_ia,antibiotic);
        J_in(i,j) = k_e1*(alpha-m_e2*y(end,2))/m_e1*substrate(i)/(substrate(i)+Km_s);
        varphi(i,j) = (k_m(j)*y(end,1))/J_in(i,j); 
        Mhat(i,j) = y(end,1);
        J_grow(i,j) = k_e2*y(end,2)*y(end,1)/(Km_m+y(end,1))*1/(1+antibiotic)*Ym;
        y0 = y(end,:);
    end
end

subplot(4,2,1);
hold on;
for i=1:length(k_m)
    plot(substrate/1e3, J_grow(:,i), 'Color', my_colors(i,:), 'LineWidth', 1);
end

% Experimental data from "Bacterial growth laws reflect the evolutionary
% importance of energy efficiency" Arijit Maitra and Ken A. Dill 2015
Glucose_concentration = [0.008,0.011,0.015,0.038,0.054,0.057,0.070,0.049,...
                         0.085,0.123,0.168,0.154,0.174,0.209,0.285,0.301,...
                         0.369,0.388,0.412,0.501,0.479,0.543,0.618,0.635,...
                         0.644,0.836,0.737,0.337];
Specific_growth_rate =  [0.030,0.084,0.139,0.287,0.376,0.448,0.400,0.536,...
                         0.783,0.849,0.860,0.789,0.811,0.848,0.816,0.818,...
                         0.815,0.815,0.814,0.818,0.913,0.898,0.834,0.842,...
                         0.855,0.871,0.921,0.957];
plot(Glucose_concentration, Specific_growth_rate, 'ko', 'MarkerSize', 8);

axis square;
box on;
xlabel('Substrate (S) concentration (mM)');
ylabel('Specific growth rate (1/h)');
set(gca,'XScale','log');
xlim([1e-3,1e0]);
set(gca,'XTick',[1e-3,1e-2,1e-1,1e0]);
ylim([0,1.0]);
set(gca,'YTick',[0:0.2:1]);
legend(strcat('k_m=',string(num2cell(k_m))), 'Location', 'bestoutside');
title('Vary Substrate (S)');

subplot(4,2,3);
hold on;
for i=1:length(k_m)
    plot(J_in(:,i)/1e3, Mhat(:,i), 'Color', my_colors(i,:), 'LineWidth', 1);
end
axis square;
box on;
xlabel('Metabolite influx (mM/h)');
ylabel('Metabolite (M) concentration (\muM)');
set(gca,'XScale','log');
xlim([1e1,1e3]);
set(gca,'XTick',[1e1,1e2,1e3]);
ylim([0,2.0]);
set(gca,'YTick',[0:1:2]);

subplot(4,2,5);
hold on;
for i=1:length(k_m)
    plot(substrate/1e3, Mhat(:,i), 'Color', my_colors(i,:), 'LineWidth', 1);
end
axis square;
box on;
xlabel('Substrate (S) concentration (mM)');
ylabel('Metabolite (M) concentration (\muM)');
set(gca,'XScale','log');
xlim([1e-3,1e0]);
set(gca,'XTick',[1e-3,1e-2,1e-1,1e0]);
ylim([0,5.0]);
set(gca,'YTick',[0:1:5]);

subplot(4,2,7);
hold on;
for i=1:length(k_m)
    plot(substrate/1e3, varphi(:,i), 'Color', my_colors(i,:), 'LineWidth', 1);
end
axis square;
box on;
set(gca,'XScale','log');
xlabel('Substrate (S) concentration (mM)');
ylabel('Metabolite leakage fraction');
xlim([1e-3,1e0]);
set(gca,'XTick',[1e-3,1e-2,1e-1,1e0]);
ylim([0,1]);
set(gca,'YTick',[0:0.2:1]);

%% vary external antibiotic concentration
substrate = 100; % uM, external substrate concentration
antibiotic = 0:1:100; % uM, external antibiotic concentration
varphi = zeros(length(antibiotic), length(k_m)); % fraction of metabolite influx that leak into environment
Mhat = zeros(length(antibiotic), length(k_m)); % metabolite concentration
J_grow = zeros(length(antibiotic), length(k_m)); % specific growth rate
J_in = zeros(length(antibiotic), length(k_m)); % influx of metabolite
J_overflow = zeros(length(antibiotic), length(k_m)); % influx of metabolite
for i=1:length(antibiotic)
    y0 = [0,0.5,0];
    for j=1:length(k_m)
        [t,y] = ode15s(@ecoli_growth_model,tspan,y0,options,substrate,k_m(j),m_e1,m_e2,k_e1,k_e2,Km_s,Km_m,Ki_r,k_r,d_r,alpha,Ym,K_ia,antibiotic(i));
        J_in(i,j) = k_e1*(alpha-m_e2*y(end,2))/m_e1*substrate/(substrate+Km_s);
        J_overflow(i,j) = k_m(j)*y(end,1);
        varphi(i,j) = (k_m(j)*y(end,1))/J_in(i,j); 
        Mhat(i,j) = y(end,1);
        J_grow(i,j) = k_e2*y(end,2)*y(end,1)/(Km_m+y(end,1))*Ym*1/(1+antibiotic(i));
        y0 = y(end,:);
    end
end

subplot(4,2,2);
hold on;
for i=1:length(k_m)
    plot(antibiotic, J_grow(:,i), 'Color', my_colors(i,:), 'LineWidth', 1);
end
axis square;
box on;
xlabel('Antibiotic (A) concentration (\muM)');
ylabel('Specific growth rate (1/h)');
xlim([0,100]);
set(gca,'XTick',[0:25:100]);
ylim([0,0.6]);
set(gca,'YTick',[0:0.2:0.6]);
title('Vary Antibiotic (A)');

subplot(4,2,4);
hold on;
for i=1:length(k_m)
    plot(J_in(:,i)/1e3, Mhat(:,i), 'Color', my_colors(i,:), 'LineWidth', 1);
end
axis square;
box on;
xlabel('Metabolite influx (mM/h)');
ylabel('Metabolite (M) concentration (\muM)');
xlim([0,2000]);
set(gca,'XTick',[0:500:2000]);
ylim([0,300]);
set(gca,'YTick',[0:100:300]);

subplot(4,2,6);
hold on;
for i=1:length(k_m)
    plot(antibiotic, Mhat(:,i), 'Color', my_colors(i,:), 'LineWidth', 1);
end
axis square;
box on;
xlabel('Antibiotic (A) concentration (\muM)');
ylabel('Metabolite (M) concentration (\muM)');
xlim([0,100]);
set(gca,'XTick',[0:25:100]);
ylim([0,100]);
set(gca,'YTick',[0:20:100]);

subplot(4,2,8);
hold on;
for i=1:length(k_m)
    plot(antibiotic, varphi(:,i), 'Color', my_colors(i,:), 'LineWidth', 1);
end
axis square;
box on;
xlabel('Antibiotic (A) concentration (\muM)');
ylabel('Metabolite leakage fraction');
xlim([0,100]);
set(gca,'XTick',[0:25:100]);
ylim([0,1]);
set(gca,'YTick',[0:0.2:1]);