%
%
%%%%%%%%%%%%%%%%   CNP2genes created by MonaIsa  %%%%%%%%%%%%%%%%%%%%%%%%
%
%  What this script does: It converts cluster information to genenames,
%  counts how often genes are hit and does multipleHitStatistics for them
%
%%%%%%%%%%%%%%
%  Input: - Cluster information: either CNPSummary, deldup or a txt file
%         - bed file; masterlist
%         - dataset: where in SNPSummary or deldup are the data of
%         interest?
%         - savepath: where are the data saved?
%       !!! You will want to exclude accessory genome genes (aka genes that are (partially) hit by acc. genome regions)
%         - Input: AccMM2Genes_Hits.mat
%
%        Additional input:
%             - you can load a txt file with a lsit of BSU names of genes that
%               interest you and do the whole analysis just for these genes, otherwise
%               it is done with all genes
%
%%%%%%%%%%%%%%%%
%  Output
%         HitOutput: struct of every hit on a gene for every replicate
%             -> struct('RepNo',[],'Cluster',[],'Type',[],'BSU',[],'Genename',[],'Frac',[],'IdentMl',[]);     
%         Type: did the hit occur .. in, over a complete, at the front or
%         tail .. of a gene
%
%          GeneOutput: struct of every gene and how it was hit for every replicate
%             -> struct('RepNo',[],'Cluster',[],'Type',[],'BSU',[],'Genename',[],'Frac',[],'IdentMl',[]);
%
%          RepSummary: summarizies information per replicate
%              -> struct('RepNo',[],'BSUHit',[],'GeneHit',[],'Replicate',[],'UniqBSU',[]);
%
%          MultiHitStat = struct('BSU',[],'Genename',[],'SumHit',[],'HitIndex',[]);
%
%
% things to keep in mind:
%    - the genes that belong to accessory genome have to be excluded from
%    further analysis. For this reason 

clear all; close all

%% Here, you have to make some inputs: 
% Define paths
masterlist = '/home/isabel/Dokumente/ExpEvol/dictionaries/BsubNC_000964wt/ml/BsubW23_2Bs166NCe_DP50_ml.txt';
recipbed = '/home/isabel/Dokumente/ExpEvol/Neisseria/dictionaries/MS11_ncbi_IR/MS11_Features_filtered.txt';
specGenes = '/home/isabel/Dokumente/ExpEvol/Compare_W_V/LeuOperon/leuOperon.txt';

% where are your replicates of interest in CNPSummary or deldup??  
% dataset= [36:40];
% 
% % Path/names/suffix of the clusterOutput:
% CNPname = '/home/isabel/Dokumente/ExpEvol/EvolExp_Ws_ns/DNA/W23_Compare_s_ns/Cluster/data/CNPSummary_Ws_cy10_20_15july2020_IR.mat';
% CNPSummary=load(CNPname);
% CNPSummary=CNPSummary.CNPSummary(dataset);

% /home/isabel/Dokumente/ExpEvol/EvolExp_Ws_ns/DNA/W23_Compare_s_ns/Cluster/CNPSummary_Wnscy10_20_DP50.mat
% /home/isabel/Dokumente/ExpEvol/EvolExp_Ws_ns/DNA/W23_Compare_s_ns/Cluster/CNPSummary_JeffW23_cy21_DP50.mat

%Path to the MultiHitStat output wth the acc and mm genes
% excAccmm = load('/home/isabel/Dokumente/ExpEvol/dictionaries/BsubNC_000964wt/acc/BsubW23_AccMM2Genes.mat');
% excAccmm = excAccmm.AccMM2Genes;
% exclude_thr = 1.1;

%Path/names/suffix of the deldup output:
% DDname = '/home/isabel/Dokumente/ExpEvol/EvolExp_Ws_ns/DNA/W23_Compare_s_ns/Del_Dup/Output_deldup_Wscy20_woRep10_15julyIR.mat';
% D=load(DDname);
% DD.del=D.deldup.del(dataset);
% DD.dup=D.deldup.dup(dataset);

% You can also input a txt file. It should have a column with start and end
% If you add a third column with a name, then each entry will count as its
% own replicate, if not, then all start and end will be taken together and
% a randomn sample name as rep will be given
txtdoc = '/home/isabel/Dokumente/ExpEvol/Neisseria/dictionaries/MS11_ncbi_IR/multimapper/Ms11_mm_regions_min40bp.txt';
fid = fopen(txtdoc);
imp = textscan(fid,'%f %f %f','delimiter','\t');
fclose(fid);
txt.S=imp{1}; txt.E=imp{2}; txt.L=imp{3};
clear imp
% 
% txtdoc2 = '/home/isabel/Dokumente/ExpEvol/EvolExp_Ws_ns/DNA/W23_Compare_s_ns/CNP2Genes/spbeta.txt';
% fid = fopen(txtdoc2);
% imp = textscan(fid,'%f %f %s','delimiter','\t');
% fclose(fid);
% txt1.S=imp{1}; txt1.E=imp{2}; txt1.N=imp{3};
% clear imp
% 
% txt.S = [txt.S;txt1.S];
% txt.E = [txt.E;txt1.E];
% txt.N = [txt.N;txt1.N];

% the plots and hotcolds are saved here
savepath = '/home/isabel/Dokumente/ExpEvol/Neisseria/dictionaries/MS11_ncbi_IR/'

%%%%%Load variables%%%%%%%%%%%%
recipsize = 2233640;
yes=1; no=0;

%% Here the data is loaded

% %Load recipient annonated bed file
fid = fopen(recipbed); bed = textscan(fid,'%s %f %f %s %s','delimiter',' ');
fclose(fid);
gene168.GN=bed{4}; gene168.S=bed{2}; gene168.E=bed{3}; gene168.BSU=bed{5}; 
gene168.L=gene168.E-gene168.S+1;
clear bed
%Load recipient/donor specific master list
fid = fopen(masterlist); imp = textscan(fid,'%f %s %s');
fclose(fid);
refchr.pos=imp{1}; refchr.atcg=imp{3}; 
clear imp
clear prompt str fid ans imp

%Do you want to look for specific genes? They can be provided as a list
quest = 'Do you want to search for specific genes that you provide in specGenes?? [yes/no] ';
spec = input(quest,'s');
if strcmp(spec,'yes')
% load special genes
fid = fopen(specGenes);
imp = textscan(fid,'%s','delimiter','\t');
fclose(fid);
specGen.N=imp{1};
clear imp fid

% and look for the indices of the special genes
[Lia,Locb] = ismember(horzcat({gene168.BSU{:}}),{specGen.N{:}});
speIdx = find(Lia);
gene168.BSU = {gene168.BSU{speIdx}}';
gene168.GN = {gene168.GN{speIdx}}';
gene168.S = gene168.S(speIdx);
gene168.E = gene168.E(speIdx);
gene168.L = gene168.L(speIdx);
end

%User input: Use which cluster information?
prompt = 'Do you want to search the genes in C, Adist, denovo, del or dup or txt? : ';
which = input(prompt,'s');

if strcmp(which,'C')
    C={CNPSummary(:).C}';
    for i=1:size(C,1)
        if ~isempty(C{i,1})
            Cluster{i,1}(1,:) = C{i,1}(1,:); Cluster{i,1}(2,:) = C{i,1}(2,:);
        else
            Cluster{i,1}(1,:) = 0; Cluster{i,1}(2,:) = 0;
        end
        rep = {CNPSummary(:).Samples};
        ORI = {CNPSummary(:).ORI_Crossing};
    end
    
elseif strcmp(which,'Adist')
    Adist={CNPSummary(:).Adist}';
    for i=1:size(Adist,1)
        if ~isempty(Adist{i,1})
            Cluster{i,1}(1,:) = Adist{i,1}(:,2)'; Cluster{i,1}(2,:) = Adist{i,1}(:,3)';
        else
            Cluster{i,1}(1,:) = 0; Cluster{i,1}(2,:) = 0;
        end
        rep = {CNPSummary(:).Samples};
        ORI = {CNPSummary(:).ORI_Crossing};
    end
    
elseif strcmp(which,'denovo')
    denovo={CNPSummary(:).denovo}';
    for i=1:size(denovo,1)
        if ~isempty(denovo{i,1})
            Cluster{i,1}(1,:) = denovo{i,1}(1,:)'; Cluster{i,1}(2,:) = denovo{i,1}(1,:)';
        else
            Cluster{i,1}(1,:) = 0; Cluster{i,1}(2,:) = 0;
        end
        rep = {CNPSummary(:).Samples};
        ORI = {CNPSummary(:).ORI_Crossing};
    end
    
elseif strcmp(which,'del')
    delstart={DD.del(:).start}';
    deledge={DD.del(:).edge}';
    for i=1:size(DD.del,2)
        if ~isempty(delstart{i,1})
            Cluster{i,1}(1,:) = delstart{i,1}(:)'; Cluster{i,1}(2,:) = deledge{i,1}(:)';
        else
            Cluster{i,1}(1,:) = 0; Cluster{i,1}(2,:) = 0;
        end 
        rep = {DD.del(:).sample};
    end
elseif strcmp(which,'dup')
    dupstart={DD.dup(:).start}';
    dupedge={DD.dup(:).edge}';
    for i=1:size(DD.dup,2)
        if ~isempty(dupstart{i,1})
            Cluster{i,1}(1,:) = dupstart{i,1}(:)'; Cluster{i,1}(2,:) = dupedge{i,1}(:)';
        else
            Cluster{i,1}(1,:) = 0; Cluster{i,1}(2,:) = 0;
        end 
        rep = {DD.dup(:).sample};
    end
elseif strcmp(which,'txt')
    txtstart=txt.S;
    txtedge=txt.E;
    f = fieldnames(txt); ff = {f{:}};
    if sum(ismember(ff,'N'))>0
        for i=1:numel(txtstart)
            Cluster{i,1}(1,:) = txtstart(i); Cluster{i,1}(2,:) = txtedge(i);
        end
        rep = {txt.N{:}};      
    else
        Cluster{1,1}(1,:) = txtstart; Cluster{1,1}(2,:) = txtedge;
        rep = {'Platzhalter'}
    end
end
clear f ff 

% in the case that we are looking at denovos, deletions and duplications, acc. genes
% are not excluded!
if strcmp(which,'del') || strcmp(which,'dup') || strcmp(which,'denovos')
    exclude_thr = 1.1;
end
%% HitOutput is generated 
%Here for every replicate i all cluster c are checked for hits with genes and all different hits are written as entry into HitOutput.
%   -> This means: genes can appear multiple times per replicate -> this is
%      fixed in GeneOutput
HitOutput = struct('RepNo',[],'Cluster',[],'Type',[],'BSU',[],'Genename',[],'Frac',[],'Sample',[]);
m = 0; 
for i=1:size(Cluster,1) % loops over all samples
    % fix the problem with ORI first
    d = 1;
    while d <= size(Cluster{i},2) % loops over all cluster in ith sample
        % if we have a cluster where start > edge, we assume ORI crossing and
        % split up this entry in 2: start(1):recipsize & 1:end(1)
        if Cluster{i}(1,d) > Cluster{i}(2,d)
            first = [1 Cluster{i}(2,1)]; last = [Cluster{i}(1,1) recipsize];
            Cluster{i}(:,1) = first; Cluster{i}(:,end+1) = last;
            disp('A cluster was found that crossed the ORI. It is cut in 2 pieces!');
        end
        d = d+1;
    end
    
    for c=1:size(Cluster{i},2) % loops over all cluster in ith sample
        ptail = find( ((Cluster{i}(1,c) > gene168.S) & (Cluster{i}(1,c) < gene168.E) & (Cluster{i}(2,c) >= gene168.E)) == 1);
        %ptMask = [ptMask ptail];
        for pt=1:numel(ptail)
            m = m+1;
            HitOutput(m).Type = 'tail';
            HitOutput(m).RepNo = i; HitOutput(m).Sample = rep{i};   HitOutput(m).Cluster = c;
            HitOutput(m).Genename = gene168.GN(ptail(pt)); HitOutput(m).BSU = gene168.BSU(ptail(pt));
            HitOutput(m).Frac = (gene168.E(ptail(pt)) - Cluster{i}(1,c) + 1)/gene168.L(ptail(pt));
        end
        complete = find( ((Cluster{i}(1,c) <= gene168.S) & (Cluster{i}(2,c) >= gene168.E)) == 1 );
        for co=1:numel(complete)
            m = m+1;
            HitOutput(m).Type = 'complete';
            HitOutput(m).RepNo = i; HitOutput(m).Sample = rep{i};    HitOutput(m).Cluster = c;
            HitOutput(m).Genename = gene168.GN(complete(co)); HitOutput(m).BSU = gene168.BSU(complete(co));
            HitOutput(m).Frac = 1;
        end
        pin = find(  ((Cluster{i}(1,c) > gene168.S) & (Cluster{i}(2,c) < gene168.E))==1  );
        for pi=1:numel(pin)
            m = m+1;
            HitOutput(m).Type = 'in';
            HitOutput(m).RepNo = i; HitOutput(m).Sample = rep{i};    HitOutput(m).Cluster = c;
            HitOutput(m).Genename = gene168.GN(pin(pi)); HitOutput(m).BSU = gene168.BSU(pin(pi));
            HitOutput(m).Frac = (Cluster{i}(2,c) - Cluster{i}(1,c) + 1)/gene168.L(pin(pi));
        end
        pfront = find( ((Cluster{i}(1,c) <= gene168.S) & (Cluster{i}(2,c) > gene168.S) & (Cluster{i}(2,c) < gene168.E))==1);
        for pf=1:numel(pfront)
            m= m+1;
            HitOutput(m).Type = 'front';
            HitOutput(m).RepNo = i; HitOutput(m).Sample = rep{i};    HitOutput(m).Cluster = c;
            HitOutput(m).Genename = gene168.GN(pfront(pf)); HitOutput(m).BSU = gene168.BSU(pfront(pf));
            HitOutput(m).Frac = (Cluster{i}(2,c) - gene168.S(pfront(pf)) + 1)/gene168.L(pfront(pf));
        end
        clear complete pin pfront ptail c co d pf pi pt 
    end
end

%% GeneOutput and RepSummary are generated 
GeneOutput = struct('RepNo',[],'Cluster',[],'Type',[],'BSU',[],'Genename',[],'Frac',[],'Sample',[]);
RepSummary = struct('RepNo',[],'BSUHit',[],'GeneHit',[],'Sample',[]);


l = 0;
for i=1:numel(rep)
    idx = find([HitOutput.RepNo] == i);
    b = {HitOutput(idx).Genename}';    BSUcollectRep = vertcat(b{:});
    for j=1:numel(idx)
        match = find(strcmp(BSUcollectRep{j},BSUcollectRep));
        matchD = match(match~=j);
        if length(match)==1
            l = l + 1;
            GeneOutput(l) = HitOutput(idx(j));
        elseif matchD>j
            l = l + 1;
            GeneOutput(l) = HitOutput(idx(j));
            %add the other gene infos to this entry
            GeneOutput(l).Frac = sum([HitOutput(idx(match)).Frac]);
            GeneOutput(l).Cluster = [HitOutput(idx(match)).Cluster];
            GeneOutput(l).Type = {HitOutput(idx(match)).Type};
        end
    end
    
    %Now for each Replicate I collect the summary
    idxUniq = find([GeneOutput.RepNo] == i);
    bb = {GeneOutput(idxUniq).BSU}';           RepSummary(i).BSUHit = vertcat(bb{:});
    bbb = {GeneOutput(idxUniq).Genename}';     RepSummary(i).GeneHit = vertcat(bbb{:});
    RepSummary(i).RepNo = i;               RepSummary(i).Sample = rep{i};
    RepSummary(i).GeneHitNo = length(idxUniq);   
end

%% Multiple HitStatistics: How often is each Gene hit over all replicates??
% Prepare gene list so that it excludes the genes that are excAccmm and 
%that are above the threshold
acm = {{excAccmm(:).BSU}'}; accmm.BSU = vertcat(acm{:});

mask = [excAccmm(:).FracMean]>=exclude_thr;
ex = {{excAccmm(mask).BSU}'}; exc.BSU = vertcat(ex{:});

mlSNPsinGene = 0;
MultiHitStat = struct('BSU',[],'Genename',[],'SumHit',[],'HitIndex',[],'FracMean',[],'mlGeneIdent',[],'Samples',[],'accxluded',[],'FracList',[]);
c = {GeneOutput(:).BSU}';          BSUcollect = vertcat(c{:});
Samplecollect = {GeneOutput(:).Sample}';   
n = 0;
% here, i goes through all the entries in the bed file
for i=1:numel(gene168.BSU)
    
    MultiHitStat(i).BSU = gene168.BSU{i};
    MultiHitStat(i).Genename = gene168.GN{i};
    
    %add the gene Identity of the i-th gene based on the masterlist to
    %MultiHitStat- if the i-th gene also appears on the excAnnMM list, then
    %calculate it a bit differently
    idx_acc = find(strcmp(gene168.BSU{i},accmm.BSU));
    if ~isempty(idx_acc)
        % here a very small number is added just to not to devide by 0
        l = gene168.L(i)*(1-(excAccmm(idx_acc).FracMean));
        snps = sum(ismember(refchr.pos,gene168.S(i):gene168.E(i)));
        MultiHitStat(i).mlGeneIdent =  1-(snps/l);
        if snps>l; snps = l; MultiHitStat(i).mlGeneIdent =  1-(snps/l); end
        if l == 0; l=0.0000001; MultiHitStat(i).mlGeneIdent =  0; end        
    else
        MultiHitStat(i).mlGeneIdent =  1-(sum(ismember(refchr.pos,gene168.S(i):gene168.E(i)))/gene168.L(i));
    end
    
% if the gene also appears in exc.BSU, mark it as excluded in the
% accxluded field
    exc_idx = find(strcmp(gene168.BSU{i},exc.BSU));
% at which index idx did the i-th gene fit the geneoutput list?
    idx = find(strcmp(gene168.BSU{i},BSUcollect));
    if ~isempty(exc_idx)
        MultiHitStat(i).accxluded = 1;
        MultiHitStat(i).SumHit = 0;
        MultiHitStat(i).FracMean = 0;
        
% if the gene was hit at least once, calculate ..
    elseif isempty(exc_idx) && ~isempty(idx)
        MultiHitStat(i).accxluded = 0;
        MultiHitStat(i).SumHit = numel(idx);
        MultiHitStat(i).HitIndex = idx;
        MultiHitStat(i).FracList = [GeneOutput(idx).Frac];
        MultiHitStat(i).FracMean = mean([GeneOutput(idx).Frac]);
        MultiHitStat(i).Samples = {Samplecollect{idx}};
        
    elseif isempty(exc_idx) && isempty(idx)
        MultiHitStat(i).accxluded = 0;
        MultiHitStat(i).SumHit = 0;
        MultiHitStat(i).FracMean = 0;
    end
    
    clear idx
end

%% Create a list of the most interesting hotspots
% save data
saveit = input('Do you want to save the HotColdGenes data? (yes/no) ');
if isempty(saveit) || saveit == 1
    cutit = input('Which minimal number of MultiHits interests you? Set cutoff: ');
    cutoff=cutit; c = 0; Liste = struct('BSU',[],'Genename',[],'SumHit',[],'HitIndex',[],'FracMean',[],'mlGeneIdent',[],'Samples',[],'accxluded',[],'FracList',[]);
    for i=1:numel(gene168.BSU)
        if MultiHitStat(i).SumHit>=cutoff
            c = c+1;
            Liste(c) = MultiHitStat(i);
        end
    end
    
    if strcmp(which,'denovo')
    Liste = rmfield(Liste,{'HitIndex','mlGeneIdent','FracMean','accxluded','FracList'});
    else
    Liste = rmfield(Liste,{'HitIndex','mlGeneIdent'});
    end
    
    T=struct2table(Liste);
    writetable(T,[savepath 'den2Genes_Wscy20.txt'],'Delimiter',' ')
else
    disp('Ok, HotColdGenes are not saved.')
end

%% Hot and cold plot
% be aware that the genes that appear here as acc. with ident=0 are the ones that 
% are excluded with the threshold

rows = 6; totalgenes = numel(gene168.L);
ident = [MultiHitStat(:).mlGeneIdent]; sumhit = [MultiHitStat(:).SumHit];
left_color = [0.1 0.1 0.1];
right_color = 'b';

figure(1)
subplot(6,1,1)
title('Hot and Cold plot')
hold on
for i=1:rows
    subplot(rows,1,i)
    region = [floor(1+(i-1)*(totalgenes/rows)):ceil((totalgenes/rows)*i)];
    yyaxis left;
    plot(region,ident(region),'Color',[left_color 0.3])
    ylim([0.8 1])
    hold on
    yyaxis right
    bar(region,sumhit(region),'BarWidth',1,'FaceColor',right_color)
    xlim([1+(i-1)*(totalgenes/rows) (totalgenes/rows)*i ])
    ylim([0 max(sumhit)+1])
    ax = gca;
    ax.YAxis(1).Color = left_color;
ax.YAxis(2).Color = right_color;
    
end
subplot(6,1,ceil(rows/2))
yyaxis left; ylabel('Identity')
yyaxis right; ylabel('# of hits')
subplot(6,1,rows)
xlabel('Gene number');

%% MultiHitStatistic Plot - raw data
if strcmp(which,'txt') ~= 1
    % in this plot the acc. genes that you defined with exclude_thr are
    % excluded    
    sumhitmask_woacc = [MultiHitStat(:).accxluded]==0;
    sumhit_woacc = sumhit(sumhitmask_woacc);
    
    figure(2);
    h1=histogram(sumhit,'BinWidth',1); hold on
    h1.BinEdges = [h1.BinEdges max(h1.BinEdges)+1] - h1.BinWidth/2;
    h2=histogram(sumhit_woacc,'BinWidth',1,'EdgeColor','r','FaceColor','none');
    h2.BinEdges = [h2.BinEdges max(h2.BinEdges)+1] - h2.BinWidth/2;
    set(gca, 'YScale', 'log')
    title('Multiple Hit Statistic')
    xlabel('Gene Hits'); ylabel('# of hits')
    xlim([-1 numel(h1.Values)+2]); ylim([1 max(h1.Values)+500])
    legend([h1 h2],{'Ws cy20','Ws cy20 wo AccGenes'})
end

allBSU = {gene168.BSU}'; allBSUU = vertcat(allBSU{:});
for i=1:numel([MultiHitStat(:).FracMean])
    idx=find(strcmp(string(MultiHitStat(i).BSU),allBSUU));
    MultiHitStat(i).GeneNum = idx;
end
