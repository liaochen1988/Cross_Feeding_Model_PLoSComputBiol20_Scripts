% By Chen Liao, Memorial Sloan Kettering Cancer Center, June 3, 2020
%
% This file plots glucose and amino acids kinetics for the coculture of
% 14-amino-acid-auxotroph community

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

%% Simulate all 14 auxotroph model
tlen = 7;
dT=0.01;
nT = 24/dT;
time = 0:dT/24:7; % day
aa_conc = zeros(7*nT+1, 14);
glu_conc = zeros(7*nT+1, 1);

G_init  = 0.002/180.156*1e6/1e-3; % uM, glucose concentration 0.2% wt/vol (0.2 g in 100g water)
N0 = ones(1,14)*1e7*1e3; % 1e7 per mL (convert mL to L)
y0 = [G_init,N0,N0,zeros(1,14)];
for day=1:tlen
    [t,y] = ode15s(@multi_aa_auxotroph_cf_model,...
                   0:dT:24,...
                   y0,...
                   odeset('RelTol',1e-6,'AbsTol',1e-6*ones(1,1+3*14),'NonNegative',ones(1,1+3*14)),...
                   Vmax_g, Km_g, Vmax_aa, Km_aa, Yield_g, Yield_aa, Rcarbon, Byp_frac, Cdr);
    glu_conc(1+(day-1)*nT:1+day*nT,:) = y(:,1);
    aa_conc(1+(day-1)*nT:1+day*nT,:) = y(:,30:43);
    y0 = y(end,:)/100;
    y0(1) = G_init;
end

%% plot comparison between our model and Mee 2014 model
figure('Name','FigS9');
colororder = varycolor(14);
hold on;
plot(time,glu_conc,'k--');
for i=1:14
    plot(time,aa_conc(:,i),'k-','Color',colororder(i,:));
end
box on;
xlabel('Day');
ylabel('Glucose/Amino acid (\muM)');
set(gca,'YScale','log');
set(gca,'YTick',[1e-5,1e0,1e5]);
axis([0,7,1e-5,1e5]);
legend({'Glu';'C';'F';'G';'H';'I';'K';'L';'M';'P';'R';'S';'T';'W';'Y'},'Location','bestoutside');