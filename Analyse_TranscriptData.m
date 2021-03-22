% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %                    Analyse_TranscriptData.m                              % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% This script:
% -- runs the function CollectDGEData_fcn and saves the structure
%    collectExps

% Plots: 
% -- Sebastians plot :)

clear all; close all

%% load the data and set variables
inPath = "/home/isabel/Dokumente/P5_ExpEvol_Ngo/RNA/2021Feb_Run2/2b_Results_DGE/";

sampleName = ["NgoG4run2_Azi15min_vs_DMSO15min","NgoG4run2_Azi60min_vs_DMSO60min",...
    "NgoG4run2_Cef15min_vs_DMSO15min","NgoG4run2_Cef60min_vs_DMSO60min"];

conditions = ["Azi15min","Azi60min","Cef15min","Cef60min"];

% annotation data from blast pipeline
%     Ann_Ngo2Sibs = '/home/isabel/Dokumente/P5_ExpEvol_Ngo/dictionaries/MS11_ncbi_IR/Ann_Ngo2Sibs/Ms11_Annfrom_FA1090_Nmen.mat';
%     ann = load(Ann_Ngo2Sibs);
%     ann = ann.outList;

    
% the GO annotation, mainly to test for how many of the genes we have this
% information
NmenGO = load('/home/isabel/Dokumente/P5_ExpEvol_Ngo/dictionaries/Nmeningitidis_MC58/NmenGO_sorted.mat');
NmenGO = NmenGO.uniqList;
FA1090GO = load('/home/isabel/Dokumente/P5_ExpEvol_Ngo/dictionaries/FA1090/FA1090GO_sorted.mat');
FA1090GO = FA1090GO.uniqList;
GONames = [NmenGO(:).geneName FA1090GO(:).geneName];
GONames = GONames(GONames~="");

% where to save to??
outPath = '/home/isabel/Dokumente/P5_ExpEvol_Ngo/RNA/2021Feb_Run2/2b_Results_DGE/plots/';

% do you want to save figures to outPath?
plotFigs = "ON";
    
%% Run the function
[collectExps] = CollectDGEData_fcn(sampleName,inPath,conditions,GONames) 


%% collect the data of the genes that you are interested in

close all

% which genes are you interested in??
%genes_oi       = ["NGFG_01150","NGFG_01501"]';
genes_oi       = [collectExps(1).expData(1:10).locustag_in_Ngo]';

% which colors should the bars have- fc and which transparency fa?
fc = [0 100/255 0; 0 100/255 0; 0 0 139/255; 0 0 139/255];
fa = [0.3 0.7 0.3 0.7];

% which colors should the stars have?
fc_sig = [fc(1,:);0 0 0;fc(3,:);0 0 0];


genes_oi_where = [];
clear genes_oi_l2fC genes_oi_where genes_oi_pvalue sig
for i=1:numel(sampleName)
[~,idx] = ismember(genes_oi,[collectExps(i).expData.locustag_in_Ngo]);
sig = strings(1,numel(idx));
genes_oi_where(:,i) = idx';
genes_oi_l2fC(:,i) = [collectExps(1).expData(idx).log2FoldChange]';
genes_oi_pvalue(:,i) = [collectExps(1).expData(idx).pvalue]';

sig(genes_oi_pvalue(:,i) <= 0.1) = "*";
sig(genes_oi_pvalue(:,i) <= 0.05) = "**";

genes_oi_signic(:,i) = sig';
end

genes_NM_oi = [collectExps(numel(sampleName)).expData(idx).gene];
genes_pro_oi = string([collectExps(numel(sampleName)).expData(idx).product]);

set(groot,'defaultAxesTickLabelInterpreter','None');  

figure(1)
set(gcf,'Renderer', 'painters', 'Position', [10 10 1200 700])
title("The most interesting genes")
GOI = bar(genes_oi_l2fC);
xticklabels(genes_NM_oi' + " " + genes_oi)
set(gca, 'XTickLabelRotation', 15)
ylim([-1 1.8])

for i=1:numel(conditions)
GOI(i).FaceColor = fc(i,:); GOI(i).FaceAlpha = fa(i);

text(GOI(i).XEndPoints-0.05,repmat(1.5,1,10),genes_oi_signic(:,i),'FontWeight','bold','Color',fc_sig(i,:));

end

legend(conditions,'NumColumns',4)

%% Do you want to plot??

if plotFigs == "ON"
   print(figure(1),'-painters','-dpng',outPath + "10MostSigGenes" + "Azi15min" +".png") 
    
end
