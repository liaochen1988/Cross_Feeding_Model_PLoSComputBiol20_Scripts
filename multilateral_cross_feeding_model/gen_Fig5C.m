clear all;
clc;

auxotroph = {'C';'F';'G';'H';'I';'K';'L';'M';'P';'R';'S';'T';'W';'Y'};
%cost = [24.7, 52.0, 11.7, 38.3, 32.3, 30.3, 27.3, 34.3, 20.3, 27.3, 11.7, 18.7, 74.3, 50.0];
tblPara = readtable('data/parameters.xls','Sheet','byproduct fraction');
Byp_frac = tblPara{:, 2:end};

figure();
axis square;
h=heatmap(strcat(char(hex2dec('0394')),auxotroph),auxotroph,Byp_frac*100);
h.Title = 'Fraction of amino acid (carbon) leakage';
h.XLabel = 'Auxotroph';
h.YLabel = 'Amino acid';
caxis([0,20]);

