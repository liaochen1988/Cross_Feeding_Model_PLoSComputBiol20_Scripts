% By Chen Liao, Memorial Sloan Kettering Cancer Center, Oct. 15, 2019
%
% This file implements a coarse-grained ecology model for multilateral
% cross-feeding between 14 engineered Escherichia coli auxotrophs, one
% unable to synthesize a particular amino acid.

% Briefly, these 14 mutants cooperate by compensating for each other's
% metabolic deficiency: each of the 14 auxotrophs can secrete all 14 amino
% acids to environment except for the one that it is auxotrophic for, and
% can also uptake the amino acid that it cannot produce from the environment.

addpath('../auxiliary_functions');

% All experimental data used in the codes were publised in Mee et al., 2014, PNAS
auxotroph = {'C';'F';'G';'H';'I';'K';'L';'M';'P';'R';'S';'T';'W';'Y'};

%% Load parameters

% kinetic rates and yields parameters
% these parameters were obtained either from literature or manual curation
tblPara = readtable('data/parameters.xls','Sheet','kinetic rates and yields');
Vmax_g   = tblPara.Vmax_g;     % Maximum glucose uptake rates, umol/h
Km_g     = tblPara.Km_g;       % Michaelis constants for glucose uptake, uM
Vmax_aa  = tblPara.Vmax_aa;    % Maximum amino acid uptake rates, umol/h
Km_aa    = tblPara.Km_aa;      % Michaelis constants for amino acid uptake, uM
Yield_g  = tblPara.Yield_g;    % Biomass yields of E. coli auxotrophs growing on glucose, 1/umol
Yield_aa = tblPara.Yield_aa;   % Biomass yields of E. coli auxotrophs growing on auxotrophic amino acids, 1/umol
Rcarbon  = tblPara.Rcarbon;    % Number of amino acid molecules produced per molecule of glucose is consumed
Cdr      = tblPara.Cdr;       % Cell death rate, 1/h

% Matrix of amino acid (byproduct) leakge fraction (percentage of carbon loss)
% Dimension: 14 secreted amino acids (rows) by 14 amino acid auxotrophs (columns)
% Each element (i,j) of the matrix is the leakage fraction of amino acid i by the auxotroph j
% these parameters were programmably fit from pairwise coculture data
tblPara = readtable('data/parameters.xls','Sheet','byproduct fraction');
Byp_frac = tblPara{:, 2:end};

%% Model training: serial dilution of all 14 auxotrophs and 4 selected 13-auxotroph communities
% Growth of 13- and 14-member cocultures was done in 3-mL cultures in a 30 °C 
% rotating drum and passaged without washing by 100-fold dilution every 24 h as 
% the cultures reach saturation. The growth and fold-growth metrics mentioned 
% throughout the text refer to the yield of the community calculated by final 
% cell density/initial density.

conditions = {'14-strain community (7 days)',...
              '13-strain community (no K)',...
              '13-strain community (no R)',...
              '13-strain community (no T)',...
              '13-strain community (no M)'};

duration        = 7;                        % 7-day experiment
initial_glucose = (0.2*10)/180.156*1e6;     % initial glucose concentration (0.2%, 0.2 g in 100 ml water), uM

figure('Name','Fig6A');

ms = 8; % marker size
lw = 2; % line width
fs = 12; % font size

titles = {'All 14';'No K';'No R';'No T';'No M'};
colororder = varycolor(14);
[ha, pos] = tight_subplot(14,5,[.02 .01],[.05 .05],[.25 .25]);
for k=1:length(conditions)
    
    % load data
    tbl_relative_abundance_rep1 = readtable('data/experiment_data.xls','Sheet',conditions{k});
    obs_cell_density = tbl_relative_abundance_rep1{:,2:end};
    sim_cell_density = zeros(size(obs_cell_density));
    
    initial_cell_density = obs_cell_density(:,1)/1e-3; % convert from per mL to per L
    sim_cell_density(:,1) = obs_cell_density(:,1);
    y0 = [initial_glucose;initial_cell_density;initial_cell_density;zeros(14,1)]; % initial condition  
    for day=1:duration
        [~,y] = ode15s(@multi_aa_auxotroph_cf_model, ... Model
                       [0,24],                       ... Time interval of integration (unit: h)          
                       y0,                           ... Initial condition
                       odeset('RelTol',1e-6,'AbsTol',1e-6*ones(1,1+3*14),'NonNegative',ones(1,1+3*14)), ... Matlab options
                       Vmax_g, Km_g, Vmax_aa, Km_aa, Yield_g, Yield_aa, Rcarbon, Byp_frac, Cdr ... Model parameters
                      );
        sim_cell_density(:,day+1) = (y(end,2:14+1)/1e3)'; % convert from per L to per mL
        y0 = y(end,:)/100; % 100 fold dilution
        y0(1) = initial_glucose; % replenish glucose to the initial level
    end
     
    for i=1:14
        axes(ha((i-1)*5+k));
        hold on;
        
        plot([0:duration],obs_cell_density(i,:),'ko','MarkerFaceColor',colororder(i,:),'MarkerSize',ms);
        plot([0:duration],abs(sim_cell_density(i,:)),'k-','Color',colororder(i,:),'LineWidth',lw);
        axis([0,duration+1,1e3,1e9]);
        box on;
        set(gca,'YScale','log');
        set(gca,'FontSize',fs);
        
        % set yaxis labels
        if(k==1)
            ylabel(auxotroph{i});
            set(gca,'Ytick',[1e3,1e6,1e9]);
            set(gca,'Yticklabel',{'10^3','10^6','10^9'});
        else
            set(gca,'yticklabel',{[]});
        end

        % set xaxis labels
        if(i==14)
            xlabel('Days');
            set(gca,'XTick',[0:8]);
            set(gca,'Xticklabel',{'0','1','2','3','4','5','6','7','8'});
        else
            set(gca,'xticklabel',{[]});
        end

        % set title of each subplot column
        if(i==1)
            title(titles{k});
        end
    end
end

%% Model testing: serial dilution of all 14 auxotrophs over 50 days

% load data
tbl_relative_abundance_rep1 = readtable('data/experiment_data.xls','Sheet','14-strain community 1 (50 days)');
tbl_relative_abundance_rep2 = readtable('data/experiment_data.xls','Sheet','14-strain community 2 (50 days)');

% simulation
days_measured = tbl_relative_abundance_rep1.Properties.VariableNames(2:end);
days_measured = str2double(cellfun(@(S) S(4:end), days_measured, 'UniformOutput', 0));
duration = max(days_measured); % days

initial_cell_density = ones(14,1)*1e7/1e-3; % convert from per mL to per L
initial_glucose = (0.2*10)/180.156*1e6;     % initial glucose concentration (0.2%, 0.2 g in 100 ml water), uM
y0 = [initial_glucose;initial_cell_density;initial_cell_density;zeros(14,1)];

sim_relative_abundance = zeros(14,length(days_measured));
sim_relative_abundance(:,1) = ones(14,1)/14;
for day=1:duration
    [~,y] = ode15s(@multi_aa_auxotroph_cf_model, ... Model
                   [0,24],                       ... Time interval of integration (unit: h)
                   y0,                           ... Initial condition
                   odeset('RelTol',1e-6,'AbsTol',1e-6*ones(1,1+3*14),'NonNegative',ones(1,1+3*14)), ... Matlab options
                   Vmax_g, Km_g, Vmax_aa, Km_aa, Yield_g, Yield_aa, Rcarbon, Byp_frac, Cdr ... Model parameters
                  );
    day_index = find(days_measured==day);
    if (~isempty(day_index))
        sim_relative_abundance(:,day_index) = (y(end,2:14+1)/1e3)';
        sim_relative_abundance(:,day_index) = sim_relative_abundance(:,day_index)/sum(sim_relative_abundance(:,day_index));
    end
    
    y0 = y(end,:)/100;
    y0(1) = initial_glucose;
end
    
% plot
figure('Name','Fig6B');
hold on;

colororder = varycolor(14);
dx = 12;
bw = 0.25; % barwitdh
h1=bar(1:dx:dx*length(days_measured)-1, sim_relative_abundance', 'stacked', 'BarWidth', bw);
h2=bar(4:dx:dx*length(days_measured)+2, tbl_relative_abundance_rep1{:,2:end}', 'stacked', 'BarWidth', bw);
h3=bar(7:dx:dx*length(days_measured)+5, tbl_relative_abundance_rep2{:,2:end}', 'stacked', 'BarWidth', bw);

for i=1:length(h1)
   h1(i).FaceColor=colororder(i,:);
   h2(i).FaceColor=colororder(i,:);
   h3(i).FaceColor=colororder(i,:);
end
set(gca,'XTick',4:dx:dx*length(days_measured)+2);
set(gca,'XTicklabel',days_measured);
box on;
xlabel('Days');
ylim([0,1]);
ylabel('Relative abundance');