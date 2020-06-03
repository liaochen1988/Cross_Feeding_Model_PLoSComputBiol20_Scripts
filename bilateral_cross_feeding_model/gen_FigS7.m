% By Chen Liao, Memorial Sloan Kettering Cancer Center, June, 2, 2019

% This file plots distribution of MCMC estimated parameter values of the
% complete model

addpath('../auxiliary_functions/Violinplot-Matlab-master/');

load('output/Example2_complete_MCMC_100000_smape');

figure();
hold on;
violinplot(chain,[],'ShowData',false);
set(gca,'YScale','log');
ylim([1e-25,1e15]);
ylabel('Values');
set(gca,'YTick',10.^[-20,-10,0,10]);
set(gca,'XTickLabel',{'K_g','V_{\Deltak, k}','V_{\Deltal,k}','V_{\Deltak,l}','V_{\Deltal,l}','\phi_{\Deltal,k}',...
                      '\phi_{\Deltak,l}','\varphi_{\Deltal,k}','\varphi_{\Deltak,l}','\gamma_k','\gamma_l','\eta_{\Deltak}','\eta_{\Deltal}'})           
xtickangle(90);

