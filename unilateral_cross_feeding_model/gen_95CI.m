% By Chen Liao, Memorial Sloan Kettering Cancer Center, June 2, 2020
%
% This file calculates 95% confidence interval of free parameters

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

N = size(Markov_chain,1);                   % Number of ?Experiments? In Data Set
yMean = mean(Markov_chain);           % Mean Of All Experiments At Each Value Of ?x?
ySEM = std(Markov_chain)/sqrt(N);     % Compute ?Standard Error Of The Mean? Of All Experiments At Each Value Of ?x?
CI95 = tinv([0.025 0.975], N-1);            % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, ySEM, CI95(:));      % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ?x?

parameters = {'V_{1,g}','V_{3,g}','V_{1,aa}','C_{1,g}','\varphi_{a}','\gamma_{a}','I_{1,a}','I_{3,a}'};
for i=1:size(yCI95,2)
    fprintf('%s: [%2.4e, %2.4e]\n',parameters{i}, yCI95(1,i)+yMean(i), yCI95(2,i)+yMean(i));
end
