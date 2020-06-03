% By Chen Liao, Memorial Sloan Kettering Cancer Center, May. 16, 2020
%
% This file simulates effects of network perturbations (downshift of
% nutrient, add antibiotics, add cheaters) on the steady state of a subset
% of the 14 amino acid auxotroph community

%% Read parameters
tblPara = readtable('data/parameters.xls','Sheet','kinetic rates and yields');
Vmax_g   = tblPara.Vmax_g;     % Maximum glucose uptake rates, umol/h
Km_g     = tblPara.Km_g;       % Michaelis constants for glucose uptake, uM
Vmax_aa  = tblPara.Vmax_aa;    % Maximum amino acid uptake rates, umol/h
Km_aa    = tblPara.Km_aa;      % Michaelis constants for amino acid uptake, uM
Yield_g  = tblPara.Yield_g;    % Biomass yields of E. coli auxotrophs growing on glucose, 1/umol
Yield_aa = tblPara.Yield_aa;   % Biomass yields of E. coli auxotrophs growing on auxotrophic amino acids, 1/umol
Rcarbon  = tblPara.Rcarbon;    % Number of amino acid molecules produced per molecule of glucose is consumed
Cdr      = tblPara.Cdr;       % Cell death rate, 1/h

tblPara = readtable('data/parameters.xls','Sheet','byproduct fraction');
Byp_frac = tblPara{:, 2:end};

figure();

%% We assume perturbation is applied at steady state of unperturbed community
G_init  = 0.002/180.156*1e6/1e-3; % uM, glucose concentration 0.2% wt/vol (0.2 g in 100g water)
tlen = 60; % 60 day is sufficient long for the community to reach steady state

N0 = ones(1,14)*1e7*1e3; % 1e7 per mL (convert mL to L)
y0 = [G_init, N0, N0, zeros(1,14)];
y_ss = zeros(length(y0), tlen+1);
y_ss(:,1) = y0';
for day=1:tlen
    [~,y] = ode15s(@multi_aa_auxotroph_cf_model,...
                   [0,24],...
                   y0,...
                   odeset('RelTol',1e-6,'AbsTol',1e-6*ones(1,1+3*14),'NonNegative',ones(1,1+3*14)),...
                   Vmax_g, Km_g, Vmax_aa, Km_aa, Yield_g, Yield_aa, Rcarbon, Byp_frac, Cdr);
    y_ss(:,day+1) = y(end,:);
    y0 = y(end,:)/100;
    y0(1) = G_init;
end
y_before_perturbation = y_ss(:,end);

%% Nutrient downshift
G_fm = 10.^[2:0.02:5]; % uM, Feed medium glucose concentration
total_density_fold_change = zeros(length(G_fm),1);
parfor j=1:length(G_fm)
    y0 = y_before_perturbation;
    y0(1) = G_fm(j);
    y0(1+[1,2,3,4,7,9,10,11,13,14]) = 0; % strains that do not belong to the stable subset are set to 0
    initial_total_cell_density = sum(y0(2:15));
    
    traj = zeros(tlen+1, 14);
    traj(1,:) = y0(2:14+1);
    for day=1:tlen
        [~,y] = ode15s(@multi_aa_auxotroph_cf_model,...
                       [0,24],...
                       y0,...
                       odeset('RelTol',1e-6,'AbsTol',1e-6*ones(1,1+3*14),'NonNegative',ones(1,1+3*14)),...
                       Vmax_g, Km_g, Vmax_aa, Km_aa, Yield_g, Yield_aa, Rcarbon, Byp_frac, Cdr);
        traj(day+1,:) = y(end,2:14+1);
        y0 = y(end,:)/100;
        y0(1) = G_fm(j);
    end
    total_density_fold_change(j) = sum(traj(end,:))/initial_total_cell_density;
end

subplot(1,3,1);
hold on;
plot(G_fm, total_density_fold_change, 'k-');
axis square;
box on;
axis([1e2,1e5,0,6]);
set(gca,'XScale','log');
set(gca,'YTick',[0,2,4,6]);
set(gca,'XTick',[1e2,1e3,1e4,1e5]);
xlabel('Glucose (uM)');
ylabel({'Fold change in population';'density of active cells'});
set (gca,'xdir','reverse');

%% Add antibiotics
index_of_antibiotic_target = [5,6,8,12]; % auxotrophs targeted by antibiotics
abx = [0:0.01:2]; % antibiotic concentration
total_density_fold_change = zeros(length(index_of_antibiotic_target),length(abx));
for i=1:length(index_of_antibiotic_target)
    parfor j=1:length(abx)
        y0 = y_before_perturbation;
        y0(1) = G_init;
        y0(1+[1,2,3,4,7,9,10,11,13,14]) = 0; % strains that are not belong to the stable subset are set to 0
        initial_total_cell_density = sum(y0(2:15));
        
        traj = zeros(tlen+1, 14);
        traj(1,:) = y0(2:14+1);
        for day=1:tlen
            [~,y] = ode15s(@multi_aa_auxotroph_cf_model_abx,...
                           [0,24],...
                           y0,...
                           odeset('RelTol',1e-6,'AbsTol',1e-6*ones(1,1+3*14),'NonNegative',ones(1,1+3*14)),...
                           Vmax_g, Km_g, Vmax_aa, Km_aa, Yield_g, Yield_aa, Rcarbon, Byp_frac, Cdr, abx(j), 1, 1, index_of_antibiotic_target(i));
            traj(day+1,:) = y(end,2:14+1);
            y0 = y(end,:)/100;
            y0(1) = G_init;
        end
        total_density_fold_change(i,j) = sum(traj(end,:))/initial_total_cell_density;
    end
end

subplot(1,3,2);
hold on;
cc = jet(length(index_of_antibiotic_target));
for i=1:length(index_of_antibiotic_target)
    plot(abx, total_density_fold_change(i,:), 'k-', 'Color', cc(i,:), 'LineWidth', 1);
end
legend({'\DeltaI';'\DeltaK';'\DeltaM';'\DeltaT'});
axis square;
box on;
axis([0,1,0,1.5]);
set(gca,'XTick',[0:0.2:1]);
set(gca,'YTick',[0:0.5:1.5]);
xlabel('Antibiotics (\muM)');
ylabel({'Fold change in population';'density of active cells'});

%% add non-producing cheater (added with amount equal to a certain fraction of total cells)
index_of_cheater = [5,6,8,12]; % index of auxotrophs whose non-producing cheater are introduced
cheater_fraction  = [0,10.^[-2:0.02:0]]; % cheater fraction
total_density_fold_change = zeros(length(index_of_cheater),length(cheater_fraction));

for i=1:length(index_of_cheater)
    parfor j=1:length(cheater_fraction)
        y0 = zeros(1,71);
        y0(1) = G_init; % initial glucose
        y0(2:15) = y_before_perturbation(2:15,end); % active cells of amino acid auxotrophs
        y0(1+[1,2,3,4,7,9,10,11,13,14]) = 0; % strains that do not belong to the stable subset are set to 0
        y0(14+1+index_of_cheater(i)) = sum(y_before_perturbation(2:15,end))*cheater_fraction(j); % active cells of amino acid auxotroph cheaters
        y0(30:43) = y_before_perturbation(30:43,end);   %  amino acids
        y0(44:57) = y_before_perturbation(16:29,end);   %  total cells of amino acid auxotrophs
        y0(56+1+index_of_cheater(i)) = sum(y_before_perturbation(2:15,end))*cheater_fraction(j); % total cells of amino acid auxotroph cheaters
        initial_total_cell_density = sum(y0(2:29));
                
        traj = zeros(tlen+1, 28);
        traj(1,:) = y0(2:28+1);
        for day=1:tlen
            [~,y] = ode15s(@multi_aa_auxotroph_cf_model_w_cheater,...
                           [0,24],...
                           y0,...
                           odeset('RelTol',1e-6,'AbsTol',1e-6*ones(1,1+5*14),'NonNegative',ones(1,1+5*14)),...
                           Vmax_g, Km_g, Vmax_aa, Km_aa, Yield_g, Yield_aa, Rcarbon, Byp_frac, Cdr);
            traj(day+1,:) = y(end,2:28+1);
            y0 = y(end,:)/100;
            y0(1) = G_init;
        end
        total_density_fold_change(i,j) = sum(traj(end,:))/initial_total_cell_density;
    end
end

subplot(1,3,3);
hold on;
cc = jet(length(index_of_cheater));
for i=1:length(index_of_cheater)
    plot(cheater_fraction*100, total_density_fold_change(i,:), 'k-', 'Color', cc(i,:), 'LineWidth', 1);
end
legend({'\DeltaI';'\DeltaK';'\DeltaM';'\DeltaT'});
axis square;
box on;
axis([0,50,0,1.5]);
set(gca,'XTick',[0:10:50]);
set(gca,'YTick',[0:0.5:1.5]);
xlabel('Cheater fraction (%)');
ylabel({'Fold change in population';'density of active cells'});
