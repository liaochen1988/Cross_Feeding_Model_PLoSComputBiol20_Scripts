% By Chen Liao, Memorial Sloan Kettering Cancer Center, June 3, 2019
%
% This file fits amino acid leakage parameters from pairwise coculture data

Byp_frac = zeros(14,14);
initial_glucose = 0.002/180.156*1e6/1e-3;   % uM, glucose concentration 0.2% wt/vol (0.2 g in 100g water)
initial_cell_density = [1e7*1e3,1e7*1e3];   % 1e7 per mL (convert mL to L)

tblPairwiseCoculture = readtable('data/experiment_data.xls','Sheet','pairwise coculture');
obsFC = tblPairwiseCoculture{:,[2,4]};
simFC = zeros(size(obsFC));
for i=1:height(tblPairwiseCoculture)
    
    % use different initial conditions
    num_repeats=5;
    best_sol = [0,0];
    lowest_resnorm = 1e10;
    for j=1:num_repeats
        [sol,resnorm,~,exitflag] = lsqnonlin(@mse_2mer_ss,... % Model
            [rand,rand],... % Initial condition
            [-Inf,-Inf],... % Lower bounds
            [Inf,Inf],... % Upper bounds
            optimoptions('lsqnonlin', 'Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-6, 'MaxFunc', Inf, 'MaxIter', Inf), ... % Matlab optioins
            initial_cell_density, initial_glucose, [tblPairwiseCoculture.FC1(i),tblPairwiseCoculture.FC2(i)]);
        assert (exitflag>0);
        
        if (resnorm < lowest_resnorm)
            best_sol = sol;
            lowest_resnorm = resnorm;
        end
    end
    
    % 10^sol is multiplication of Yield_aa, Rcarbon and Byp_frac
    first_strain_index = find(strcmp(auxotroph, tblPairwiseCoculture.Strain1(i)));
    second_strain_index = find(strcmp(auxotroph, tblPairwiseCoculture.Strain2(i)));
    Byp_frac(first_strain_index,second_strain_index) = 10.^best_sol(1)/Yield_aa(first_strain_index)/Rcarbon(first_strain_index);
    Byp_frac(second_strain_index,first_strain_index) = 10.^best_sol(2)/Yield_aa(second_strain_index)/Rcarbon(second_strain_index);
    
    % get simulated values
    [~,N1ss,N2ss] = mse_2mer_ss(best_sol, initial_cell_density, initial_glucose, [tblPairwiseCoculture.FC1(i),tblPairwiseCoculture.FC2(i)]);
    simFC(i,:) = [N1ss/initial_cell_density(1), N2ss/initial_cell_density(2)];
end

% plot prediction vs data
figure('Name','Fit byproduct fraction using simplified model');
hold on;
plot(obsFC(:,1),simFC(:,1),'k.','MarkerSize',16);
plot(obsFC(:,2),simFC(:,2),'k.','MarkerSize',16);
plot([1e-2,1e2],[1e-2,1e2],'k--');
axis([1e-2,1e2,1e-2,1e2]);
xlabel('Observed fold change');
ylabel('Predicted fold change');
axis square;
box on;
set(gca,'XScale','log');
set(gca,'YScale','log');