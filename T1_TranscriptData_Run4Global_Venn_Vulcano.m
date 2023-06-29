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
set(groot,'defaultAxesFontname','CMU Serif')
set(0,'DefaultTextFontname', 'CMU Serif')
set(groot,'defaultAxesFontSize',11)
%% load the data and set variables
%basePath = "/mnt/AA684269684233FB/sciebo/P5_onSciebo/";
basePath = "/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/";
inPath = basePath + "RNA/2022March_Run4/1b_Analysis_naiv_woMS11_DpilT/2b_Results_DGE/";

expName = "2022-11-03NgoG4Run4_"+ ["G4pilEKO" "G4pilENG17" "G4pilENG24" "G4pilENG32"] + "_vs_G4wt";
conditions = ["pilEKO" "NG17" "NG24" "NG32"];
reference = "DG4 biofilm";

% we need the list of multimapper genes
MMGenes = basePath + "/dictionaries/1A_MS11_ManualLists/MULTIMAPPER_Ms11toMs11/MultiMapper_MS11_2_MS11.mat";
% add additional mm genes - here, pilE and pilS
AddMMGenes = ["NGFG_01821" "NGFG_02431" "NGFG_02484" "NGFG_02482" "NGFG_02405" "NGFG_02481" "NGFG_02253" "NGFG_00014" "NGFG_02485"...
    "NGFG_02487" "NGFG_01819" "NGFG_01818"];

% % stuff that is needed for GO enrichtment
% the GO annotation, mainly to test for how many of the genes we have this
% information
NmenGO = load(basePath + "/dictionaries/4_Nmeningitidis_MC58/NmenGO_sorted.mat");
FA1090GO = load(basePath +"/dictionaries/3_FA1090/FA1090GO_sorted.mat");


excMM           = "ON";
saveCollectExps = "OFF";  % do you want to save collectExps to outPath?
saveFigs        = "OFF";  % do you want to save figures to outPlots?
useGO           = "OFF"; % Go terms are added but ignored

plotVulcano     = "ON";

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

%% are the GOHit genes unique? -- otherwise there is a problem with TopGO
if useGO == "ON"
    collectExps_wGO = rmfield(collectExps,'expData');
    UniquenessTest = zeros(4,1);

    for i=1:numel(collectExps)

        mask_GOHit = [collectExps(i).expData(:).GOHit] > 0;
        collectExps_wGO(i).expData = collectExps(i).expData(mask_GOHit);
        collectExps_wGO(i).expData = rmfield(collectExps_wGO(i).expData,{'padj','locustag_in_Nmen','locustag_in_FA1090'});
        % now still we have to test, if the GOHit info2GO names are unique!
        UniquenessTest(i) = numel(unique([collectExps_wGO(i).expData(:).info2GO])) == numel([collectExps_wGO(i).expData(:).info2GO]);
    end
end

%% now if you want save the data separately

if saveCollectExps == "ON"

    for i=1:numel(collectExps)
        T = [];
        T = struct2table(collectExps(i).expData);

        writetable(T,[outPath + datestr(now,'yyyymmdd') + "_DGEGenes_woMM_" +  reference + "_vs_" + string(conditions(i))  + ".csv"],'Delimiter','\t')
    end
end



%% collect data for the Venn diagram

for i=1:numel(conditions)

    Venndata(i).condition = conditions(i);

    if collectExps(i).condition == conditions(i)
        Venndata(i).mask_up = [collectExps(i).expData.log2FoldChange] >= upL2FC & [collectExps(i).expData.padj] <= signifi;
        Venndata(i).mask_dw = [collectExps(i).expData.log2FoldChange] <= dwL2FC & [collectExps(i).expData.padj] <= signifi;

        Venndata(i).ltag_up = [collectExps(i).expData(Venndata(i).mask_up).locustag_in_Ngo];
        Venndata(i).ltag_dw = [collectExps(i).expData(Venndata(i).mask_dw).locustag_in_Ngo];

        Venndata(i).upFrac  = numel(Venndata(i).ltag_up)/numel(Venndata(i).mask_dw);
        Venndata(i).dwFrac  = numel(Venndata(i).ltag_dw)/numel(Venndata(i).mask_dw);
    else
        error('something is wrong ...')
    end

end

% % vulcano plot as a test

if plotVulcano == "ON"

    for i=1:numel(expName)
        figure(100+i)
        %subplot(2,2,i)
        set(gcf,'Renderer', 'painters', 'Position', [10 10 800 600])
        set(gca,'YScale','lin');
        hold on
        title("Condition:" + " " + conditions(i) + " vs " + reference)
        % problem, if genes have pvalue 0, the log10(0) = Inf, then we have
        % a problem ...

        pV_log10    = -log10([collectExps(i).expData(:).padj]);
        pV_log10MAX = max(pV_log10(~isinf(pV_log10)));
        YLim    = pV_log10MAX *1.1;
        pV_log10(find(isinf(pV_log10))) = YLim;

        s(i) = scatter([collectExps(i).expData(:).log2FoldChange],pV_log10,4,'Marker','o',...
            'MarkerFaceColor',fc(5,:),'MarkerEdgeColor',fc(5,:),'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0,'DisplayName','all genes');

        %red for up-regulated
        u(i) = scatter([collectExps(i).expData(Venndata(i).mask_up).log2FoldChange],pV_log10(Venndata(i).mask_up),4,'Marker','o',...
            'MarkerFaceColor',red,'MarkerEdgeColor','none','DisplayName','up');
        % blue for down regulated
        d(i) = scatter([collectExps(i).expData(Venndata(i).mask_dw).log2FoldChange],pV_log10(Venndata(i).mask_dw),4,'Marker','o',...
            'MarkerFaceColor',green,'MarkerEdgeColor','none','DisplayName','down');
        %also plot the genes of interest as defined in the paragraph abov

        plot([dwL2FC dwL2FC],[0 YLim],'k--','HandleVisibility','off')
        plot([upL2FC upL2FC],[0 YLim],'k--','HandleVisibility','off')
        plot([-5 5],-log10([signifi signifi]),'k--','HandleVisibility','off')

        upFrac(i)   = round(100 * Venndata(i).upFrac,1);
        downFrac(i) = round(100 * Venndata(i).dwFrac,1);

        % find symmetric xlim
        padj_NAN_mask{i} = ~isnan(s(i).YData);
        XLim = max(abs(min(s(i).XData(padj_NAN_mask{i}))),abs(max(s(i).XData(padj_NAN_mask{i}))));
        text((XLim + 0.5)/2,YLim*0.9 ,string(upFrac(i)) + " " + "%",'Color',red)
        text((-XLim - 0.5)/2 - 0.4,YLim*0.9 ,string(downFrac(i)) + " " + "%",'Color',blue)

        %legend('Location','north','NumColumns',1)

        ylim([0 YLim])
        xlim([-XLim XLim])
        %xlabel('log_2 Fold Change')
        %ylabel('-log_{10} adjusted p-value')

        if saveFigs == "ON"
            print(figure(10+i),'-painters','-dpng',fileName_Vul(i))
            % saveas(figure(10+i),fileName_Vul(i),'.emf') -> only works on windows
        end
    end
end

%% make the Venn diagram

% first prepare the data
compareVenn =  ["pilEKO" "NG17" "NG32"];

for i=1:numel(compareVenn)

    getVennUp(i).refCond = compareVenn(i);
    getVennDw(i).refCond = compareVenn(i);
    mask_ref =  [Venndata.condition] == getVennUp(i).refCond;

    getVennUp(i).totalUp = numel([Venndata(mask_ref).ltag_up]);
    getVennDw(i).totalDw = numel([Venndata(mask_ref).ltag_dw]);


    for j=1:numel(compareVenn)
        mask_compareto = [Venndata.condition] == compareVenn(j);
        
        getVennUp(i).(["sharedw"] + compareVenn(j)) = sum(contains(Venndata(mask_ref).ltag_up,Venndata(mask_compareto).ltag_up));
        getVennDw(i).(["sharedw"] + compareVenn(j)) = sum(contains(Venndata(mask_ref).ltag_dw,Venndata(mask_compareto).ltag_dw));

        % get the all three Info
        maskUp{j}  = contains(Venndata(mask_ref).ltag_up,Venndata(mask_compareto).ltag_up);
        maskDw{j}= contains(Venndata(mask_ref).ltag_dw,Venndata(mask_compareto).ltag_dw);
        
        clear mask_compareto
    end

    allThreeUp(i) = sum((maskUp{1} & maskUp{2}) & maskUp{3});
    allThreeDw(i) = sum((maskDw{1} & maskDw{2}) & maskDw{3});
    clear mask_ref maskUp maskDown
end


figure(1)
venn([getVennUp.totalUp],...
    [[getVennUp(1).("sharedw" + compareVenn(2))] [getVennUp(1).("sharedw" + compareVenn(3))] [getVennUp(2).("sharedw" + compareVenn(3))],...
    allThreeUp(1)],'FaceColor',{colors(1,:),colors(2,:),colors(4,:)},'EdgeColor','black','FaceAlpha',{0.7,0.7,0.7});
legend(compareVenn)
title('upregulated genes compared to DG4 biofilm ...')


figure(2)
venn([getVennDw.totalDw],...
    [[getVennDw(1).("sharedw" + compareVenn(2))] [getVennDw(1).("sharedw" + compareVenn(3))] [getVennDw(2).("sharedw" + compareVenn(3))],...
    allThreeDw(1)],'FaceColor',{colors(1,:),colors(2,:),colors(4,:)},'EdgeColor','black','FaceAlpha',{0.7,0.7,0.7});
legend(compareVenn)
title('downregulated genes compared to DG4 biofilm ...')

%% save the data for the other script
% outPut = "/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/RNA/2022March_Run4/1b_Analysis_naiv_woMS11_DpilT/3_Results_Transk/matlabOutputs/";
% save(outPut + datestr(now,'yyyymmdd') + "_" + "collectExps.mat",'collectExps')
% save(outPut + datestr(now,'yyyymmdd') + "_" + "Venndata.mat" ,'Venndata')

%% save the plots in the right size!
% x_width=7 ;y_width=6;
% savePath = "/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/RNA/2022March_Run4/1b_Analysis_naiv_woMS11_DpilT/3_Results_Transk/plots/";
% for i=1:4
%     figure(100 + i)
%     set(gcf,'PaperUnits','centimeters')
%     set(gcf, 'PaperPosition', [0 0 x_width y_width]);
%     print(figure(100+i),'-dsvg','-painters',savePath + "20221104_Volcano_" + conditions(i) + "_vs_DG4Biofilm.svg")
% end

%% get pilE expression

pilEexp = collectExps_all;
for i=1:4
pilEexp(i).expData = pilEexp(i).expData([pilEexp(i).expData.locustag_in_Ngo] == "NGFG_01821");
end
