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

GOList = '/home/isabel/Dokumente/ExpEvol/Ngo/dictionaries/Nmeningitidis_MC58/Nmen_GO_editted.txt';
species = 'Nmen';

outPath = '/home/isabel/Dokumente/ExpEvol/Ngo/dictionaries/Nmeningitidis_MC58/';

fid = fopen(GOList); go = textscan(fid,'%s %s %s %s %s','delimiter','\t');
fclose(fid);

%%
% go through all lines of go and crate a uniq list in the form 
% geneName \t GOID1,GOID2,GOID3...

uniqList = struct('geneName',{},'GO',[],'product',[]);
searchFields = fieldnames(uniqList); 

m=0; 
for i=1:numel(go{1})
    % the first condition always gives the first blast hit
    if i==1 || ~strcmp(go{2}{i},go{2}{i-1})
        m = m+1; k = 1;
        uniqList(m).geneName = string(go{2}{i});
        uniqList(m).GO(k) = string(go{3}{i});
        uniqList(m).product(k) = string(go{4}{i});
    elseif strcmp(go{2}{i},go{2}{i-1})
        k = k+1;
        uniqList(m).GO(k) = string(go{3}{i});
        uniqList(m).product(k) = string(go{4}{i});
    end
    
end

%% save the data

save([outPath species 'GO_sorted' '.mat'],'uniqList')
% as table in .txt
T=struct2table(rmfield(uniqList,{'product'}));
writetable(T,[outPath species 'GO_sorted' '.txt'],'Delimiter','\t')


