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
set(groot, 'defaultAxesTickLabelInterpreter','none'); set(groot, 'defaultLegendInterpreter','none');
%% load the data and set variables
basePath = "/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/";

collectedData = basePath + "RNA/2022March_Run4/1b_Analysis_naiv_woMS11_DpilT/3_Results_Transk/matlabOutputs/";
dataName      = "20221104_collectExps_wCat";

conditions = ["pilEKO" "NG17" "NG24" "NG32"];%"Cef15min", ,"Azi15min", "DMSO60min"];
reference = "DG4 biofilm";

% % set the cut-offs
upL2FC  = 0.5;
dwL2FC = -0.5;
signifi = 0.05;

% which title to plot
catName = "Genetic Information Processing";

 % which category are you interested in?
catOI = "Genetic Information Processing"; %"Phage Associated";%
% "Genetic Information Processing"; 
% "Cellular Processes/Organismal Systems/Human Disease"; %
% "Energy Metabolism"; %
% "Carbohydrate Metabolism" ;
% or give a list of genes that should be plotted!!

% category that is not KEGG but rather "selfmade", with FA1090 IDs
excelIn  = "OFF";
if excelIn == "ON"
    IntGenesType = "pilGenes";% "Glycolysis genes";% "Oxidative phosphorylation"; ;% "Oxidative phosphorylation";
    % % which single genes do you want to plot as bars?
    excelSheet = basePath +  "dictionaries/1A_MS11_ManualLists/pilGenes.xlsx";% "dictionaries/3_FA1090/ "dictionaries/3_FA1090/KEGG_Categories.xlsx";% 
    sheet = "pilGenes"; %"C-Glycolysis"; %% %"pilGenes"; %"C-Oxidative_phosphorylation";% "C-Oxidative_phosphorylation";%"pilGenes";
    x1Range = 'A1:ZZ165';
    [num,txt,raw] = xlsread(excelSheet,sheet,x1Range);
    clear path; clear name; clear sheet; clear x1Range;
    listOI_ltFA = string(txt(:,1))';

end
% or give a list
listOI_ltNgo = [];%listOI_ltFA;%[];%"NGFG_" + ["00826" "01150" "00584" "01491" "00464"]; 
%listOI_ltFA = [];
%excelIn  = "OFF";


%% load data

collectExps_wCat = load(collectedData + dataName);
collectExps_wCat = collectExps_wCat.collectExps_wCat;


%% cluster the data

% automatically capture interesting genes
% what criteria??

% get the gene names from one of the conditions, egal which one ... 
allLT_Ngo     = [collectExps_wCat(1).expData.locustag_in_Ngo];
allLT_FA      = string([collectExps_wCat(1).expData.locustag_in_FA1090]);
allName       = [collectExps_wCat(1).expData.gene];
allCategories = [collectExps_wCat(1).expData.KEGGCat];

collect4Heat = struct('Samples',[],'Log2Fold',[],'locustag',[]);
collect4Heat.locustag = cellstr(allLT_Ngo)';
collect4Heat.gene = cellstr(allName)';

for i=1:numel(conditions)
    collect4Heat.Samples{i,1}    = char(conditions(i));

    [~,idx_Cl] = ismember(allLT_Ngo,[collectExps_wCat(i).expData.locustag_in_Ngo]);
    collect4Heat.Log2Fold(:,i) = single([collectExps_wCat(i).expData(idx_Cl).log2FoldChange])';

    check = all(allLT_Ngo == [collectExps_wCat(i).expData(idx_Cl).locustag_in_Ngo]);
    if check ~=1
        error('something went wrong!!')
    end

end

% if we are only interested in 1 category, then we do this here:

if ~isempty(catOI)
    % create a mask for the category of interest
    mask_catOI = allCategories == catOI;

    collect4Heat_OI = collect4Heat;
    collect4Heat_OI.locustag  = collect4Heat_OI.locustag(mask_catOI,:);
    collect4Heat_OI.Log2Fold  = collect4Heat_OI.Log2Fold(mask_catOI,:);
    collect4Heat_OI.gene      = collect4Heat_OI.gene(mask_catOI,:);


end

if ~isempty(listOI_ltNgo)
    % create a mask for the category of interest
    mask_listOI = (sum(allLT_Ngo' == listOI_ltNgo,2)>0);

    collect4Heat_OI = collect4Heat;
    collect4Heat_OI.locustag  = collect4Heat_OI.locustag(mask_listOI,:);
    collect4Heat_OI.Log2Fold  = collect4Heat_OI.Log2Fold(mask_listOI,:);
    collect4Heat_OI.gene      = collect4Heat_OI.gene(mask_listOI,:);

end

if excelIn == "ON"
    % create a mask for the category of interest
    mask_listOI = (sum(allLT_FA' == listOI_ltFA,2)>0);

    collect4Heat_OI = collect4Heat;
    collect4Heat_OI.locustag  = collect4Heat_OI.locustag(mask_listOI,:);
    collect4Heat_OI.Log2Fold  = collect4Heat_OI.Log2Fold(mask_listOI,:);
    collect4Heat_OI.gene      = collect4Heat_OI.gene(mask_listOI,:);

end
%% plot
% plot all the data
% cg_all = clustergram(collect4Heat.Log2Fold, 'RowLabels', collect4Heat.locustag,...
%     'ColumnLabels', collect4Heat.Samples,...
%     'RowPdist', 'correlation',...
%     'ColumnPdist', 'correlation');%,...
%    % 'ImputeFun', @knnimpute);
% 
% cg_all.Colormap = redbluecmap;

% plot only one category
cg_catOI= clustergram(collect4Heat_OI.Log2Fold, 'RowLabels', [[string(collect4Heat_OI.gene)] + " " + [string(collect4Heat_OI.locustag)]],...
    'ColumnLabels', collect4Heat_OI.Samples,...
    'RowPdist', 'correlation',...
    'ColumnPdist', 'correlation');%,...
   % 'ImputeFun', @knnimpute);


cg_catOI.Colormap = parula(30);
%title = addTitle(cg_catOI,catOI);
title = addTitle(cg_catOI,catName);


print(figure(1),'-dpng','-painters',"/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/RNA/2022March_Run4/1b_Analysis_naiv_woMS11_DpilT/3_Results_Transk/plots/ClusterHeatmap.png")


k=4
sum([collectExps_wCat(k).expData.log2FoldChange] <= -0.5 & [collectExps_wCat(k).expData.padj] <= 0.05 &...
    [collectExps_wCat(k).expData.KEGGCat]=="Genetic Information Processing" &...
    contains([collectExps_wCat(k).expData.product],"ribosomal"))

sum([collectExps_wCat(k).expData.log2FoldChange] <= -0.5 & [collectExps_wCat(k).expData.padj] <= 0.05 &...
    [collectExps_wCat(k).expData.KEGGCat]=="Genetic Information Processing")