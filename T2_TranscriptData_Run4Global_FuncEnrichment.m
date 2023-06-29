% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %                    T_TranscriptData_Run4_FuncEnrichment.m           % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% This script:
% -- takes the output from T_TranscriptData_Run4_Venn_Vulcano.m: Venndata &
% collectExp and does the functional enrichment


clear all; close all

%% load the data and set variables
%basePath = "/mnt/AA684269684233FB/sciebo/P5_onSciebo/";
basePath = "/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/RNA/2022March_Run4/1b_Analysis_naiv_woMS11_DpilT/";

conditions = ["pilEKO" "NG17" "NG24" "NG32"];%"Cef15min", ,"Azi15min", "DMSO60min"];
reference = "DG4 biofilm";

collectPath = basePath + "3_Results_Transk/matlabOutputs/20221104_collectExps.mat";
VennPath = basePath + "3_Results_Transk/matlabOutputs/20221104_Venndata.mat";

% where we get the categories from:
GeneCatIn   = "/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/dictionaries/GeneCategories/msystems.00729-19-st003Genco2020_2MATLAB.xlsx";
T_GeneCat = readtable(GeneCatIn);

% % predefined colors
% colors
fc = [179 225 172; 0 82 33; 254 179 67; 250 69 10; 150 150 150]/255; % set
fa = [0.7 0.7 0.7 0.7 0.7];
gray = fc(5,:);
red  = [139 0 0]/255;
blue = [0 0 204]/255;
colors = [255 255 255; 221 170 51; 0 68 136; 187 85 102; 0 0 0]/255;

%% load data

collectExps = load(collectPath);
collectExps = collectExps.collectExps;

Venndata = load(VennPath);
Venndata = Venndata.Venndata;

%% connect the collectExp List with Categories and evaluate
collectExps_wCat = collectExps;

% are the conditions the same and in the right order?
if any(conditions ~= [collectExps.condition])
    error("conditions dont match!!! check this")
end

for x=1:numel(conditions)                                      % go over all experiments
    for g=1:numel([collectExps(x).expData.log2FoldChange])  % and over all entries
        FA1090LT = []; splitUp = []; FA1090LT_edit = [];
        
        FA1090LT      = string(collectExps(x).expData(g).locustag_in_FA1090{1});
        
        % can be that there is no ortholog
        if FA1090LT == ""
            collectExps_wCat(x).expData(g).KEGGCat = "no_FA1090_Ortholog";
            
        else
            splitUp = split(FA1090LT,'_');
            FA1090LT_edit = splitUp(1) + splitUp(2);
            LT_idx = find(strcmp(FA1090LT_edit,T_GeneCat.LocusTag)); % search the gene in the table with GeneCat
            
            if isempty(LT_idx)
                collectExps_wCat(x).expData(g).KEGGCat = "no_GeneCat_found";
            else
                collectExps_wCat(x).expData(g).KEGGCat = string(T_GeneCat.KEGGCategory(LT_idx));
            end
        end
        
    end
end

%% genes for which to do enrichment

% for example the intersection of up regulated genes in the 3 planktonic conditions 
sample_cond = ["pilEKO"];% "pilEKO" "NG17" "NG32"];
regulation = ["dw" "up" "both"]; % there is dw, up and both is ["dw" "up"]

Venndata_All = Venndata; % overwrite the actual variable but save it here
Venndata = Venndata(find(contains([Venndata.condition],sample_cond)));

% find the list of genes that intersect here!

% when intersect between 3 conditions

geneList.sample = sample_cond;
for i=1:2
    geneList.(regulation(i) + "_lt_Ngo") = Venndata.("ltag_" + regulation(i));
end
geneList.(regulation(3) + "_lt_Ngo") = [geneList.(regulation(1) + "_lt_Ngo") geneList.(regulation(2) + "_lt_Ngo")];

% we have the lt in Ngo, what is the index in FA1090, for which we have the
% categories??
for i=1:numel(regulation)
geneList_mask        = contains([collectExps_wCat(1).expData.locustag_in_Ngo],geneList.(regulation(i) + "_lt_Ngo")');
geneList.(regulation(i) + "_lt_FA1090") = collectExps_wCat(1).expData(geneList_mask);
clear geneList_mask
end

%% collect the category information

% rRNA is not included, because no genes are detected for this category
% new cat:
% % -- no_FA1090_Ortholog
% % -- no_GeneCat_found


% write information in the table
T_CollectCat                  = table(unique([collectExps_wCat(1).expData.KEGGCat]'),'VariableNames',{'KEGGCat'});
T_CollectCat.Genes_in_Cat     = sum([collectExps_wCat(1).expData.KEGGCat]==T_CollectCat.KEGGCat,2);
for i=1:numel(regulation)
    GeneIntersect = [geneList.(regulation(i) + "_lt_FA1090")];
    T_CollectCat.(regulation(i) + "_in_Cat") = sum(T_CollectCat.KEGGCat==[GeneIntersect.KEGGCat],2);
end
clear GeneIntersect

%% make pie charts

% pie([T_CollectCat.Genes_in_Cat])
% lgd = legend([T_CollectCat.KEGGCat],'Interpreter','none');


%% fisher's exact test

% we already implement Bonferroni:
Bonf_factor = 19; %numel(T_CollectGeneCat_Azi60.Genes_in_Cat);

%fixed values:
allGenes = numel([collectExps_wCat(1).expData]);

for r=1:numel(regulation)
    allDiffGenes    = sum(T_CollectCat.(regulation(r) + "_in_Cat"));
    allNotDiffGenes = allGenes - allDiffGenes;

    % for down regulated genes
    for i=1:numel([T_CollectCat.Genes_in_Cat])
        a = T_CollectCat.(regulation(r) + "_in_Cat")(i);
        b = T_CollectCat.Genes_in_Cat(i) - a;
        c = allDiffGenes - a;
        d = allNotDiffGenes - b;

        if a+b ~= T_CollectCat.Genes_in_Cat(i)
            error("your contingency table is not contingent!!!!")
        end

        tt = table([a;c],[b;d],'VariableNames',{'diff','not diff'},'RowNames',{'in Cat','not in Cat'});
        [h,p_both] = fishertest(tt,'Tail','both');
        T_CollectCat.((regulation(r)) + "_padj_both")(i) = p_both * Bonf_factor;
        [h,p_left] = fishertest(tt,'Tail','left');
        T_CollectCat.((regulation(r)) + "_padj_left")(i) = p_left * Bonf_factor;
        [h,p_right] = fishertest(tt,'Tail','right');
        T_CollectCat.((regulation(r)) + "_padj_right")(i) = p_right * Bonf_factor;
        clear a b c d tt h
    end
    clear allDiffGenes allNotDiffGenes
end


% give an output what is significant!!
if ~isempty(find(T_CollectCat.dw_padj_left < 0.05))
disp("In " + T_CollectCat.KEGGCat(find(T_CollectCat.dw_padj_left < 0.05)) + ", sig. little genes are downregulated!")
end
if ~isempty(find(T_CollectCat.dw_padj_right < 0.05))
disp("In " + T_CollectCat.KEGGCat(find(T_CollectCat.dw_padj_right < 0.05)) + ", sig. many genes are downregulated!")
end
if ~isempty(find(T_CollectCat.up_padj_left < 0.05))
disp("In " + T_CollectCat.KEGGCat(find(T_CollectCat.up_padj_left < 0.05)) + ", sig. little genes are upregulated!")
end
if ~isempty(find(T_CollectCat.up_padj_right < 0.05))
disp("In " + T_CollectCat.KEGGCat(find(T_CollectCat.up_padj_right < 0.05)) + ", sig. many genes are upregulated!")
end
if ~isempty(find(T_CollectCat.both_padj_left < 0.05))
disp("In " + T_CollectCat.KEGGCat(find(T_CollectCat.both_padj_left < 0.05)) + ", sig. little genes are both down- and upregulated!")
end
if ~isempty(find(T_CollectCat.both_padj_right < 0.05))
disp("In " + T_CollectCat.KEGGCat(find(T_CollectCat.both_padj_right < 0.05)) + ", sig. many genesare both down- and upregulated!")
end