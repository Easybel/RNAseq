% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %                   GOList_to_formR.m                              % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% 

%Q59623	sodC	GO:0005507	Superoxide dismutase [Cu-Zn]	NCBITaxon:122586
%Q59623	sodC	GO:0004784	Superoxide dismutase [Cu-Zn]	NCBITaxon:122586
%Q59623	sodC	GO:0042597	Superoxide dismutase [Cu-Zn]	NCBITaxon:122586
%Q59623	sodC	GO:0019430	Superoxide dismutase [Cu-Zn]	NCBITaxon:122586
%P0C277	folD	GO:0035999	Bifunctional protein FolD	NCBITaxon:122586
%P0C277	folD	GO:0005829	Bifunctional protein FolD	NCBITaxon:122586
%P0C277	folD	GO:0004488	Bifunctional protein FolD	NCBITaxon:122586

clear all; close all


GOList{1} = '/home/isabel/Dokumente/P5_ExpEvol_Ngo/dictionaries/Nmeningitidis_MC58/Nmen_GO_editted.txt';
species{1} = 'Nmen';


GOList{2} = '/home/isabel/Dokumente/P5_ExpEvol_Ngo/dictionaries/FA1090/FA1090_GO_editted.txt';
species{2} = 'FA1090';

outPath = '/home/isabel/Dokumente/P5_ExpEvol_Ngo/dictionaries/';

prompt = 'Do you want to get the GO list for sepcies 1, 2 or [1 2]?'
variant = input(prompt);

if numel(variant)==1
    fid = fopen(GOList{variant}); go_collect = textscan(fid,'%s %s %s %s %s','delimiter','\t');
    fclose(fid);
elseif numel(variant)>1
    fid = fopen(GOList{variant(1)}); go_collect = textscan(fid,'%s %s %s %s %s','delimiter','\t');
    fclose(fid);
    fid = fopen(GOList{variant(2)});
    go = textscan(fid,'%s %s %s %s %s','delimiter','\t');
    for i=1:5
        go_collect{i} = vertcat(go_collect{i},go{i});
    end
    fclose(fid);
end

%%
% go through all lines of go and crate a uniq list in the form 
% geneName \t GOID1,GOID2,GOID3...

uniqList = struct('geneName',{},'GO',[],'product',[]);
searchFields = fieldnames(uniqList); 

m=0; 
for i=1:numel(go_collect{1})
    % the first condition always gives the first blast hit
    if i==1 || ~any(strcmp(go_collect{2}{i},{go_collect{2}{1:i-1}}))
        m = m+1; k = 1;
        uniqList(m).geneName = string(go_collect{2}{i});
        uniqList(m).GO(k) = string(go_collect{3}{i});
        uniqList(m).product(k) = string(go_collect{4}{i});
    elseif any(strcmp(go_collect{2}{i},{go_collect{2}{1:i-1}}))
        idx = min( find (cellfun(@(x) strcmp(go_collect{2}{i},x), {uniqList(:).geneName})));
        uniqList(idx).GO(end+1) = string(go_collect{3}{i});
        uniqList(idx).product(end+1) = string(go_collect{4}{i});
    end
    
end

% etwas nachjustieren: 
% ---- die GO terms sollen sich nicht wiederholen
% ---- die GO terms sollen mit Kommas separiert sein
fileID = fopen([outPath species{1} '_' species{2} '_' 'GO_sorted.txt'],'w');
for i=1:size(uniqList,2)
    l = {uniqList(i).GO};
    uni = unique(l{1});
    uniqList(i).GO = uni;
    
    % save the data
    formatSpec = ['%s\t' repmat('%s,',1,numel(uni)-1) '%s\n'];
    fprintf(fileID,formatSpec,uniqList(i).geneName,uni);
end

%% save the data

save([outPath species{1} '_' species{2} '_' 'GO_sorted.mat'],'uniqList')
% as table in .txt
T=struct2table(rmfield(uniqList,{'product'}));
