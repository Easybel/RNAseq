% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %                    Analyse_TranscriptData.m                              % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% This script:
% -- runs the function CollectDGEData_fcn and saves the structure
%    collectExps

% Plots:
% -- Sebastian's plot :)

clear all; close all
addpath '/home/isabel/Documents/Doktorarbeit_Mai2022/SCRIPTS/RNA'
addpath '/home/isabel/Documents/Doktorarbeit_Mai2022/SCRIPTS/Plots_kleinesPaper'

%% load the data and set variables
%basePath = "/mnt/AA684269684233FB/sciebo/P5_onSciebo/";
basePath = "/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/";
inPath = basePath + "RNA/2022March_Run4/4_pilE_Expression/mapping/minO20/";

expName = "2022-11-03NgoG4Run4_pilEexpr_minO20_"+ ["G4pilEKO" "G4pilENG17" "G4pilENG24" "G4pilENG32"] + "_vs_G4wt";
conditions = ["pilEKO" "NG17" "NG24" "NG32"];%"Cef15min", ,"Azi15min", "DMSO60min"];
reference = "DG4 biofilm";

% we need the list of multimapper genes
MMGenes = basePath + "/dictionaries/1A_MS11_ManualLists/MULTIMAPPER_Ms11toMs11/MultiMapper_MS11_2_MS11.mat";
% don't add any additional genes!
AddMMGenes = [];

% get phages!
phagePath = "/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/dictionaries/1A_MS11_ManualLists/MS11_PhageList.mat";
MS11BedPath = "/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/dictionaries/1A_MS11_ManualLists/MS11_onlyCDS.bed.mat"; 
% % stuff that is needed for GO enrichtment
% the GO annotation, mainly to test for how many of the genes we have this
% information
NmenGO = load(basePath + "/dictionaries/4_Nmeningitidis_MC58/NmenGO_sorted.mat");
FA1090GO = load(basePath +"/dictionaries/3_FA1090/FA1090GO_sorted.mat");


excMM           = "ON";
saveCollectExps = "OFF";  % do you want to save collectExps to outPath?
saveFigs        = "OFF";  % do you want to save figures to outPlots?
useGO           = "OFF"; % Go terms are added but ignored

% % predefined colors
% colors
fc = [179 225 172; 0 82 33; 254 179 67; 250 69 10; 150 150 150]/255; % set
fa = [0.7 0.7 0.7 0.7 0.7];
gray = fc(5,:);
red  = [165 0 38]/255;
blue = [0 0 204]/255;
green = [27 120 55]/255;
colors = [255 255 255; 221 170 51; 0 68 136; 187 85 102; 0 0 0]/255;

% % set the cut-offs
upL2FC  = 0.5;
dwL2FC = -0.5;
signifi = 0.05;


%% handle the data

NmenGO = NmenGO.uniqList;
FA1090GO = FA1090GO.uniqList;
GONames = [NmenGO(:).geneName FA1090GO(:).geneName];
GONames = GONames(GONames~="");

if excMM == "ON"
    MMGenes = load(MMGenes);
    MMGenes = MMGenes.blastMatch_mm;
    List_MMGenes = [MMGenes.locustag_in_Query];
    List_MMGenes = [List_MMGenes AddMMGenes];
end

% load the phages
phageList = load(phagePath);
phageList = phageList.phageList;

% find the genes that belong to the phages!
MS11Bed = load(MS11BedPath);
MS11Bed = MS11Bed.uniqList_CDS;
MS11Start = [MS11Bed.start];
MS11Ende  = [MS11Bed.ende];

for i=1:numel(phageList)
    mask_inPhage = phageList(i).startMS11 <= MS11Start & phageList(i).endeMS11 >= MS11Ende;
    phageList(i).Genes = MS11Bed(mask_inPhage);
clear mask_inPhage
end


%% Run the function
[collectExps] = T_CollectDGEData_fcn(expName,inPath,conditions,GONames);

%% exclude Genes
%%% -- exclude the multimapper genes that were identified before, because
%%% the multimapper reads are also excluded :)
if excMM == "ON"
    collectExps_woMM     = rmfield(collectExps,{'expData'});
    collectExps_MMGenes  = rmfield(collectExps,{'expData'});

    for i=1:numel(collectExps)

        mask_mmGenes  =  contains([collectExps(i).expData(:).locustag_in_Ngo],List_MMGenes);

        collectExps_woMM(i).expData    = collectExps(i).expData(~mask_mmGenes);
        collectExps_MMGenes(i).expData = collectExps(i).expData(mask_mmGenes);

    end

    % it is less confusing if you work in the following with the list wo MM
    % Genes

    collectExps_all = collectExps;
    collectExps     = collectExps_woMM;
end


%%  collect the expression of the phages!
collectPhageExp = [];
for i=1:numel(collectExps_woMM)
    collectPhageExp(i).condition = collectExps_woMM(i).condition;
    collectPhageExp(i).phage     = phageList;
    for j=1:numel(collectPhageExp(i).phage)

        for g = 1:numel([collectPhageExp(i).phage(j).Genes])
            collectPhageExp(i).phage(j).lt(g) = [collectPhageExp(i).phage(j).Genes(g).locus_tag];
            idx_lt = find([collectExps_woMM(i).expData.locustag_in_Ngo] == collectPhageExp(i).phage(j).lt(g)) ;
            if isempty(idx_lt)
                collectPhageExp(i).phage(j).log2FC(g)   = nan;
                collectPhageExp(i).phage(j).padj(g)     = nan;
            else
                collectPhageExp(i).phage(j).log2FC(g)   = collectExps_woMM(i).expData(idx_lt).log2FoldChange;
                collectPhageExp(i).phage(j).padj(g)     = collectExps_woMM(i).expData(idx_lt).padj;
            end
        end
        mask_woMM = ~isnan(collectPhageExp(i).phage(j).log2FC);
        collectPhageExp(i).phage(j).ltwoMM   = collectPhageExp(i).phage(j).lt(mask_woMM);
        collectPhageExp(i).phage(j).log2FCwoMM   = collectPhageExp(i).phage(j).log2FC(mask_woMM);
        collectPhageExp(i).phage(j).padjwoMM     = collectPhageExp(i).phage(j).padj(mask_woMM);
    end
end

%% prepare the data for a box plot!

for i=1:numel([collectPhageExp.condition])
    collectPhageExp(i).BoxPlot.dataLog  = [collectPhageExp(i).phage.log2FCwoMM];
    collectPhageExp(i).BoxPlot.ID       = [];
    collectPhageExp(i).BoxPlot.IDrandom = [];
    for j=1:numel(collectPhageExp(i).phage)
        collectPhageExp(i).BoxPlot.ID       = [collectPhageExp(i).BoxPlot.ID repmat(j,1,numel(collectPhageExp(i).phage(j).ltwoMM))];
        collectPhageExp(i).BoxPlot.IDrandom = [collectPhageExp(i).BoxPlot.IDrandom ];
    end
    collectPhageExp(i).BoxPlot.padj  = [collectPhageExp(i).phage.padjwoMM];
end


for i=1:numel([collectPhageExp.condition])

figure(i); hold on 
pMask = collectPhageExp(i).BoxPlot.padj <= 0.05;

plot([collectPhageExp(i).BoxPlot.ID(pMask)-0.2] + rand(1,sum(pMask))/3,collectPhageExp(i).BoxPlot.dataLog(pMask),'Marker','.',...
    'MarkerEdgeColor',red,'MarkerFaceColor',red,'LineStyle','none','MarkerSize',8)
plot([collectPhageExp(i).BoxPlot.ID(~pMask)-0.2] + rand(1,sum(~pMask))/3,collectPhageExp(i).BoxPlot.dataLog(~pMask),'Marker','.',...
    'MarkerEdgeColor',gray,'MarkerFaceColor',gray,'LineStyle','none','MarkerSize',8)
boxplot(collectPhageExp(i).BoxPlot.dataLog,collectPhageExp(i).BoxPlot.ID)
ax = gca;
ax.XTickLabel = ["NgoPhi 1" "NgoPhi 2" "NgoPhi 3" "NgoPhi 4" "NgoPhi 5"];
    
title(collectPhageExp(i).condition)


end