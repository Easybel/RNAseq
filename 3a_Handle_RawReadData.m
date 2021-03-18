%% Daten vom 2august
%all data that shows 4hcellsin4hSs
clear all; close all
red = [1 0.698 0.698];
blue = [0.698 0.698 1];
gray = [0.8 0.8 0.8]
clr = {'k','c','b','g','y','r','m',blue,gray};

path = '/home/isabel/Dokumente/P5_ExpEvol_Ngo/RNA/2021Feb_Run2/';

name = 'NgoG4_Run2_rRNA.xlsx'
%sheet = 'DiffSpecies_inMedium';
sheet = 'Conclusion';
x1Range = 'A1:N19';
[num,txt,raw] = xlsread(strcat(path,name),sheet,x1Range);
clear path; clear name; clear sheet; clear x1Range; 

% read in data
data = struct();
cat = txt(1,1:end);
data.(cat{1}) = txt(2:end,1);
data.(cat{2}) = txt(2:end,2);

for i=1:length(cat)-2
    data.(cat{i+2}) = num(:,i);
end

% plot

% is the starting RNA concentration connected to the number of reads?
figure(1)
scatter(data.RNA_concetration,data.fastqRaw_Reads,'Displayname','fastqReads','MarkerFaceColor',clr{1},'MarkerEdgeColor',clr{1})
hold on
scatter(data.RNA_concetration,data.uniqMap_reads,'Displayname','uniquely mapped reads','MarkerFaceColor',clr{2},'MarkerEdgeColor',clr{2},'Marker','d')
scatter(data.RNA_concetration,data.MultiMap_reads,'Displayname','multi mapped reads','MarkerFaceColor',clr{3},'MarkerEdgeColor',clr{3},'Marker','*')
legend()
xlabel('RNA concentration in µg/ml'); ylabel('reads')

% 
figure(2)
scatter(data.fastqRaw_Reads,data.uniqMap_reads./data.fastqRaw_Reads,'Displayname','uniquely mapped reads','MarkerFaceColor',clr{2},'MarkerEdgeColor',clr{2},'Marker','d')
hold on 
scatter(data.fastqRaw_Reads,data.MultiMap_reads./data.fastqRaw_Reads,'Displayname','multi mapped reads','MarkerFaceColor',clr{3},'MarkerEdgeColor',clr{3},'Marker','p')
scatter(data.fastqRaw_Reads,(data.uniqMap_reads+data.MultiMap_reads)./data.fastqRaw_Reads,'Displayname','sum (uniquely + multi mapped reads)','MarkerFaceColor','r','MarkerEdgeColor','r','Marker','p')
legend()
xlabel('# raw reads'); ylabel('reads/# fastq raw reads')
ylim([0 2.2])

% 
figure(3)
scatter(data.uniqMap_reads,data.assCDS_featureC,'Displayname','assigned mRNA counts','MarkerFaceColor',clr{1},'MarkerEdgeColor',clr{1})
hold on
scatter(data.uniqMap_reads,data.unassMM_Exon_featureC,'Displayname','rRNA proxy -- multimapper','MarkerFaceColor',clr{2},'MarkerEdgeColor',clr{2},'Marker','d')
legend()
xlabel('mRNA reads'); ylabel('reads')

% other influences

 %
figure(11)
boxplot(data.fastqRaw_Reads,data.Day)
xlabel('day of treatment'); ylabel('raw reads')

figure(12)
boxplot(data.unassMM_Exon_featureC,data.Day)
xlabel('day of treatment'); ylabel('proxy for rRNA')
figure(121)
boxplot(data.assCDS_featureC,data.Day)
xlabel('day of treatment'); ylabel('mRNA counts')

%
figure(13)
boxplot(data.fastqRaw_Reads,data.Condition)
xlabel('treatment'); ylabel('raw reads')


figure(14)
boxplot(data.unassMM_Exon_featureC,data.Condition)
xlabel('treatment'); ylabel('proxy for rRNA')


%
figure(15)
boxplot(data.fastqRaw_Reads,data.Time)
xlabel('duration of treatment'); ylabel('raw reads')


figure(16)
boxplot(data.unassMM_Exon_featureC,data.Time)
xlabel('duration of treatment'); ylabel('proxy for rRNA')
