% By Chen Liao, Memorial Sloan Kettering Cancer Center, June 2, 2020
%
% This file calculates 95% confidence interval of free parameters

load('output/Example2_simplified_MCMC_100000_smape.mat')

N = size(chain,1);                   % Number of ?Experiments? In Data Set
yMean = mean(chain);                 % Mean Of All Experiments At Each Value Of ?x?
ySEM = std(chain)/sqrt(N);           % Compute ?Standard Error Of The Mean? Of All Experiments At Each Value Of ?x?
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ?x?

parameters = {'K_g','V_{\Deltak, k}','V_{\Deltal,l}','\phi_{\Deltak,l}',...
              '\phi_{\Deltal,k}','\gamma_k','\gamma_l','\eta_{\Deltak}','\eta_{\Deltal}'};
for i=1:size(yCI95,2)
    fprintf('%s: [%2.4e, %2.4e]\n',parameters{i}, yCI95(1,i)+yMean(i), yCI95(2,i)+yMean(i));
end

