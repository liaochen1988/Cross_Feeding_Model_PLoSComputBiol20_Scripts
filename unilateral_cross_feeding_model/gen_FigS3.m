% By Chen Liao, Memorial Sloan Kettering Cancer Center, June, 2, 2019

% This file plots distribution of MCMC estimated parameter values using (1)
% Violin plot and (2) Pairwise scatter plot

addpath('../auxiliary_functions/Violinplot-Matlab-master/');

load('output/Example1_simplified_MCMC_20000_smape_part1');
Markov_chain_part1 = chain;
load('output/Example1_simplified_MCMC_20000_smape_part2');
Markov_chain_part2 = chain;
load('output/Example1_simplified_MCMC_20000_smape_part3');
Markov_chain_part3 = chain;
load('output/Example1_simplified_MCMC_20000_smape_part4');
Markov_chain_part4 = chain;
load('output/Example1_simplified_MCMC_20000_smape_part5');
Markov_chain_part5 = chain;
Markov_chain = [Markov_chain_part1;...
                Markov_chain_part2;...
                Markov_chain_part3;...
                Markov_chain_part4;...
                Markov_chain_part5];

figure('Name','FigS3A');
hold on;
violinplot(Markov_chain,[],'ShowData',false);
set(gca,'YScale','log');
ylim([1e-15,1e15]);
ylabel('Values');
set(gca,'YTick',10.^[-10,0,10]);
set(gca,'XTickLabel',{'V_{1,g}','V_{3,g}','V_{1,aa}','C_{1,g}',...
                      '\varphi_{a}','\gamma_{a}','I_{1,a}','I_{3,a}'})           
manual_curated_values = [5.52e-13, 6.28e-13, 4.26e-13, 0.87, 0.33, 2e11, 26.3, 112.7];
plot(1:length(manual_curated_values),manual_curated_values,'rx');

figure('Name','FigS3B');
[S,AX,BigAx,H,HAx] = plotmatrix(Markov_chain);