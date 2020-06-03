% By Chen Liao, Memorial Sloan Kettering Cancer Center, June 3, 2020
%
% This file simulates the effect of removing cross feeding links on community integrity

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

%% Remove cross feeding links
auxotroph = {'I';'K';'M';'T'};
index_of_targets = [5,6,8,12];

tlen = 30;

figure('Name','FigS10');
for i=1:length(index_of_targets)
    for j=1:length(index_of_targets)
        [i,j]
        % turn off secretion of amino acid j by species i
        
        G_init  = 0.002/180.156*1e6/1e-3; % uM, glucose concentration 0.2% wt/vol (0.2 g in 100g water)
        N0 = ones(1,14)*1e7*1e3; % 1e7 per mL (convert mL to L)
        y0 = [G_init,N0,N0,zeros(1,14)];
        ys = zeros(14,tlen+1);
        ys(:,1) = N0';
        
        Byp_frac_copy = Byp_frac;
        Byp_frac_copy(index_of_targets(j),index_of_targets(i)) = 0;
        for day=1:tlen
            [~,y] = ode15s(@multi_aa_auxotroph_cf_model,...
                           [0,24],...
                           y0,...
                           odeset('RelTol',1e-6,'AbsTol',1e-6*ones(1,1+3*14),'NonNegative',ones(1,1+3*14)),...
                           Vmax_g, Km_g, Vmax_aa, Km_aa, Yield_g, Yield_aa, Rcarbon, Byp_frac_copy, Cdr);
            ys(:,day+1) = y(end,2:14+1)/1e3;
            y0 = y(end,:)/100;
            y0(1) = G_init;
        end
        
        subplot(4,4,(i-1)*4+j);
        title(sprintf('Delta%s->%s',auxotroph{i},auxotroph{j}));
        hold on;        
        colororder = jet(4);
        for k=1:length(index_of_targets)
            plot(0:tlen,log10(ys(index_of_targets(k),:)),'k-','Color',colororder(k,:),'LineWidth',1.5);
        end
        axis square;
        box on;
        ylim([0,10]);
        xlim([0,30]);
        set(gca,'XTick',[0,10,20,30]);
        set(gca,'YTick',[0,5,10]);
        set(gca,'ticklength',6*get(gca,'ticklength'));
        xlabel('Days');
        ylabel('log_{10} (Cell density) (cells/mL)');
    end
end