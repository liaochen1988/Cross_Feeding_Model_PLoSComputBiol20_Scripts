% By Chen Liao, Memorial Sloan Kettering Cancer Center, June 2, 2019

% This file plots distribution of MCMC estimated parameter values using (1)
% Violin plot and (2) Pairwise scatter plot

addpath('../auxiliary_functions/Violinplot-Matlab-master/');
load('output/Example2_simplified_MCMC_100000_smape');

figure('Name','FigS8A');
hold on;
violinplot(chain,[],'ShowData',false);
set(gca,'YScale','log');
ylim([1e-15,1e15]);
ylabel('Values');
set(gca,'YTick',10.^[-10,0,10]);
set(gca,'XTickLabel',{'K_g','V_{\Deltak, k}','V_{\Deltal,l}','\phi_{\Deltak,l}',...
                      '\phi_{\Deltal,k}','\gamma_k','\gamma_l','\eta_{\Deltak}','\eta_{\Deltal}'})           
manual_curated_values = [1.75, 8.35e-14, 1.22e-13, 3.2e-3, 1.39e-2, 5.72e12, 2.53e12, 0.102, 4.00e-4];
plot(1:length(manual_curated_values),manual_curated_values,'rx');
xtickangle(90);

figure('Name','FigS8B');
[S,AX,BigAx,H,HAx] = plotmatrix(chain);