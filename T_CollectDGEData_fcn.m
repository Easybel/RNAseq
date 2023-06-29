function [collectExps] = T_CollectDGEData_fcn(sampleName,inPath,conditions,GONames) 

%%%%%%%%%%%%%%%%%%%%%%%%% Function to collect DGE data %%%%%%%%%%%%%%%%%%%
%% What this function does %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT:
% This function takes the path to your output DGE data that you generate
% with R and it takes the names of the files and the according conditions.
% Then it reads in all the data into a big structure and as a ...
% OUTPUT:
% gives you:
%%% -- locustag_in_Ngo 
%%% -- all DGE expression info
%%% -- GOHit: is there GO info available for this gene in one of the
%%%    species?
%%% -- info2GO: the gene name or locustag that fit the GO information that
%%%    is available. 
%%% -- the genename and locustag (from where there is GO info) 
%%% -- product from Ngo 


%% here the script is run
for s=1:numel(sampleName)
    
    TranscriptData = inPath + sampleName(s) + ".tab";
    
    fid = fopen(TranscriptData);
    varNames = strsplit(fgetl(fid), '\t');
    fclose(fid);
    
    %transc data
    
    fid = fopen(TranscriptData); transcriptome = textscan(fid,'%s %f %f %f %f %f %f %s %s %s %s %s %s','delimiter','\t','HeaderLine',1);
    fclose(fid);
    
    transc.locustag_in_Ngo = transcriptome{1};
    
    for i=1:numel(varNames)
        a = varNames{i};
        fieldN = string(a);
        transc.(fieldN)=transcriptome{i+1};
    end
    
    % fix the problem with the last rows
    transcNames = fieldnames(transc);
    transc_range = (1:numel(transc.padj));
    for i=1:numel(transcNames)
        transc.(transcNames{i}) =  transc.(transcNames{i})(transc_range);
    end
    
    
    %% actually what we want is a list of all genes and their expression
    expData = struct('locustag_in_Ngo',num2cell(transc.locustag_in_Ngo),'baseMean',num2cell(transc.baseMean),'log2FoldChange',num2cell(transc.log2FoldChange),...
        'pvalue',num2cell(transc.pvalue),'padj',num2cell(transc.padj),'info2GO',[],'GOHit',[],'gene',[],'locustag',[],'product',num2cell(transc.product),...
        'locustag_in_Nmen',num2cell(transc.locustag_in_Nmen),'locustag_in_FA1090',num2cell(transc.locustag_in_FA1090));
    foi = {'Name','product'};
    
    % take the empty expData list and write the locustags of Ngo into it
    for i=1:numel(transc.locustag_in_Ngo)
        
        splitname1 = regexp(transc.locustag_in_Ngo{i},'-','split');
        splitname2 = regexp(splitname1{2},'"','split');
        expData(i).locustag_in_Ngo = string(splitname2{1});
        
    end
    
    % now which information may best be used to connect it to the GO ann
    % find the individually best match for each entry
    
    searchGOinfo = {'locustag_in_Nmen','GN_in_Nmen','locustag_in_FA1090','GN_in_FA1090'};
    
    hit_where = zeros(numel([expData.locustag_in_Ngo]),1);
    
    GOHits = strings(numel([expData.locustag_in_Ngo]),1);
    
    for i=1:numel(searchGOinfo)
        
        foi_search = cellfun(@(x) string(x),transc.(searchGOinfo{i}));
        
        [GO_yes,GO_where] = ismember(foi_search,GONames);
        
        err = sum((foi_search(GO_yes) == GONames(GO_where(GO_yes))') ==0);
        if err ~=0
            error("something is wrong here")
        end
        
        hit_where(:,i) = GO_where;
    end
    
    
    for i =1:numel(transc.locustag_in_Ngo)
        % initialize
        expData(i).info2GO  = "";
        expData(i).gene     = "";
        expData(i).locustag = "";
        
        expData(i).GOHit = sum([hit_where(i,:)] ~= 0) > 0;
        
        if hit_where(i,1) ~= 0
            expData(i).info2GO = GONames(hit_where(i,1));
        elseif hit_where(i,2) ~= 0
            expData(i).info2GO = GONames(hit_where(i,2));
        elseif hit_where(i,3) ~= 0
            expData(i).info2GO = GONames(hit_where(i,3));
        elseif hit_where(i,4) ~= 0
            expData(i).info2GO = GONames(hit_where(i,4));
        end
        
        if hit_where(i,1) ~= 0 || hit_where(i,2) ~= 0
            expData(i).gene = string(transc.GN_in_Nmen{i});
            expData(i).locustag = string(transc.locustag_in_Nmen{i});
        elseif hit_where(i,3) ~= 0 || hit_where(i,4) ~= 0
            expData(i).gene = string(transc.GN_in_FA1090{i});
            expData(i).locustag = string(transc.locustag_in_FA1090{i});
        end
        
        if expData(i).gene=="" && ~isempty(transc.GN_in_Ngo{i})
            expData(i).gene = string(transc.GN_in_Ngo{i});
        elseif expData(i).gene=="" && ~isempty(transc.GN_in_Nmen{i})
            expData(i).gene = string(transc.GN_in_Nmen{i});
         elseif expData(i).gene=="" && ~isempty(transc.GN_in_FA1090{i})
            expData(i).gene = string(transc.GN_in_FA1090{i});
        end
    end
    
    
    
    % write expData into another structure
    collectExps(s).condition = conditions(s);
    collectExps(s).expData = expData;

 
end
