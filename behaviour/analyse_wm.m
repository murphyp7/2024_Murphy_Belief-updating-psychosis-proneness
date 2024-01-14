clear, close all

% Path/subject stuff
wmdatpath = '/mnt/homes/home028/gmonov/SCZ/WM_analysis/clean_WM_data/';

addpath '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/Gen_fun';
addpath(genpath('/mnt/homes/home024/pmurphy/Surprise_scz/'))

allsubj = {'t010';'t011';'t012';'t013';'t014';'t015';'t016';'t017';'t018';'t019';...
           't020';'t024';'t025';'t027';'t028';'t029';...
           't030';'t031';'t032';'t033';'t034';'t035';'t036';'t037';'t038';'t039';...
           't040';'t041';'t043';'t044';'t045';'t047';'t048';'t049';...
           't050';'t051';'t052';'t053';'t054';'t055';'t056';'t057';'t058';'t059';...
           't060';'t061';'t062';'t063';'t064';'t065';'t066';'t067';'t068';'t069';...
           't070';'t071';'t072';'t073';'t074';'t075';'t076';'t078';'t079';...
           't080';'t081';'t082';'t083';'t084';'t085';'t086';'t087';'t088';...
           't090';'t091';'t092';'t093';'t094';'t095';'t096';'t097';'t098';'t099';...
           't100';'t101';'t102';'t103';'t104';'t105'};
% Bad participants that are currently excluded: t022, t023, t042, t046
% Following participants with no CAPE scores excluded: t077, t089
% Participants with no WM: t021, t026
       
nsubj = length(allsubj);

modeltype = 'H_noise';

% Loop through subjects
for subj = 1:nsubj
    % Compute WM accuracy
    load([wmdatpath,allsubj{subj},'clean_allbehav.mat'])
    WMacc(subj,1) = sum(allbehav(:,6))./size(allbehav,1).*100;
end


% Figure/analysis settings
jitbnd = [-0.3 0.3];
scatsize = 9;
axlw=0.5;
colL = [0.5 0.5 0.5];
colH = [0 0 0];

nperm=10000;


% Pull CAPE scores
[CAPE,names] = get_CAPE_new(allsubj);
CAPE(:,strcmp(names,'P')) = (CAPE(:,strcmp(names,'P'))+20)./20;  % rescaling scores so that choice scale is 1-4 (rather than coded 0-3) and total score is average over all items
CAPE(:,strcmp(names,'N')) = (CAPE(:,strcmp(names,'N'))+14)./14;
CAPE(:,strcmp(names,'D')) = (CAPE(:,strcmp(names,'D'))+8)./8;


% Plot relationship of WM accuracy to scores on each CAPE
fig_w = 13; % figure width
fig_h = 8.5; % figure height

hf =  findobj('type','figure');
cfig = length(hf)+1;
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),
sp=1;

for sc = 1:length(names)
    cuts = prctile(CAPE(:,sc),[20 80]);
    gL = find(CAPE(:,sc)<=cuts(1));
    gH = find(CAPE(:,sc)>=cuts(2));
    

    subplot(2,3,sp), hold on
    rnds = dist_scatters(WMacc(gL),0.1);  % get x jitters
    S=scatter(ones(size(WMacc(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),WMacc(gL),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
    plot(jitbnd+1,[mean(WMacc(gL)) mean(WMacc(gL))],'LineWidth',2,'Color',colL)

    rnds = dist_scatters(WMacc(gH),0.1);  % get x jitters
    S=scatter(ones(size(WMacc(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),WMacc(gH),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
    plot(jitbnd+2,[mean(WMacc(gH)) mean(WMacc(gH))],'LineWidth',2,'Color',colH)
    
    p = permutationTest(WMacc(gL), WMacc(gH), nperm);
    text(1.1,95,['p=',num2str(round(p,3))],'FontSize',7)
    
    set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{sc}],['High ',names{sc}]})
    ylabel('WM accuracy (%)'), xlim([0.5 2.5])
    
    
    subplot(2,3,sp+length(names)), hold on
    [rho,p]=corr(WMacc,CAPE(:,sc),'type','spearman','rows','pairwise');
    scatter(CAPE(:,sc),WMacc,scatsize), ylabel('WM accuracy (%)'), xlabel(names{sc})
    title(['rho=',num2str(round(rho,3)),', P=',num2str(round(p,3))])
    set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off')
    sp=sp+1;
end
