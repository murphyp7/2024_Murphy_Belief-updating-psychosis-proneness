clear, close all

% Path/subj stuff
modelpath = '/mnt/homes/home024/pmurphy/Surprise_scz/modelling/';
datpath = '/mnt/homes/home028/gmonov/SCZ/Data/decision_making/';

addpath '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/Gen_fun';
addpath(genpath('/mnt/homes/home024/pmurphy/Surprise_scz/'))

savepath = '/mnt/homes/home024/pmurphy/Surprise_scz/analysis/';

allsubj = {'t010';'t011';'t012';'t013';'t014';'t015';'t016';'t017';'t018';'t019';...
           't020';'t021';'t024';'t025';'t026';'t027';'t028';'t029';...
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
% N.B. 't048' uses extreme last-sample weighting
% N.B. 't071' has previous diagnosis od drug-induced psychosis
       
nsubj = length(allsubj);

modeltype = 'H_noise';  % model fit variant
CPPsmps = 4:10;   % sample positions over which to average CPP kernels
UNCsmps = 3:8;   % sample positions over which to average uncertainty kernels

% Loop through subjects
GA_LLR=[]; GAsurprise=[]; GAunc=[]; GAchoice=[]; GALLRcount=[]; GAacc=[]; GAacc_sw=[]; GAacc_i=[]; GAacc_i_sw=[]; GAacc_fit=[]; GAacc_perf=[]; GAacc_last=[]; GAntrials=[]; GA_pmfun={}; GA_pmfun_sum={}; GA_pmfun_end={};
for subj = 1:nsubj
    LLR_full=[]; LLRsum_full=[]; LLRend_full=[]; LPR_full=[]; LPR_full_i=[]; surprise_full=[]; surprise_full_i=[]; psi_full=[]; psi_full_i=[]; choices_full=[]; choices_all=[]; acc_full=[]; fdist_full=[];
    LPR_4acc_full=[]; LPR_fit_4acc_full=[]; LLRsum_4acc_full=[]; LLRend_4acc_full=[]; nswitch=[];
    
    % Load data for current subject
    fprintf('Pulling behavioural data for subject %s...\n',allsubj{subj})
    
    fsess = dir([datpath,allsubj{subj},filesep]);
    LLRin=[]; choices=[]; nblock=0;
    for s = 1:length(fsess)-2
        fblock = dir([datpath,allsubj{subj},filesep,fsess(s+2).name,filesep,'Behaviour',filesep,'*.mat']);
        fsmp = dir([datpath,allsubj{subj},filesep,fsess(s+2).name,filesep,'Sample_seqs',filesep,'*.mat']);
        for b = 1:length(fblock)
            load([fblock(b).folder,filesep,fblock(b).name])
            load([fsmp(b).folder,filesep,fsmp(b).name])
            
            % Converting sample and choice values to appropriate signs for choice regressions (silly inconsistency in signs of these variables in task code makes this necessary)
            stimIn = round(stimIn.*-1);
            choices = Behav(:,2)-1;
            acc = Behav(:,3);
            
            if subj==1 && s==1 && b==1, maxsamps = size(stimIn,2); end
            if s==1 && b==1, GA_H_i(subj,1) = gen.H; end
            
            % Convert stimulus values to true LLRs & 
            LLRin = log(normpdf(stimIn,gen.mu(2)*-1,gen.sigma(2))./normpdf(stimIn,gen.mu(1)*-1,gen.sigma(1)));
            
            % Calculate LPRs & sample-wise surprise using model fits
            pIn = cat(3,normpdf(stimIn,gen.mu(2)*-1,gen.sigma(2)),normpdf(stimIn,gen.mu(1)*-1,gen.sigma(1)));
            [LPRout_i,surprise_i,psi_i] = accGlaze_fast(LLRin,gen.H,0,'pCP',pIn);  % ideal observer
            
            if strcmp(modeltype,'H_noise')
                load([modelpath,modeltype,filesep,'fits',filesep,allsubj{subj},'_fit.mat'])
                GA_H(subj,1) = pm_fit(1);
                GA_noise(subj,1) = pm_fit(2);
                
                [LPRout,surprise,psi] = accGlaze_fast(LLRin,GA_H(subj),0,'pCP',pIn);
                
            elseif strcmp(modeltype,'H_noise_IU')
                load([modelpath,modeltype,filesep,'fits',filesep,allsubj{subj},'_fit.mat'])
                GA_H(subj,1) = pm_fit(1);
                GA_noise(subj,1) = pm_fit(2);
                GA_IU(subj,1) = pm_fit(3);
                
                [LPRout,surprise,psi] = accGlaze_InconUp_fast(LLRin,GA_H(subj),GA_IU(subj),0,'pCP',pIn);
            end
            surprise = log(surprise./(1-surprise));  % logit-transform CPP
            
            % Isolating useable full-sequence trials
            ts=[]; tsCE=[];
            for t = 1:length(choices)
                if sum(isnan(stimIn(t,:)))==0 && choices(t)<2, ts(end+1) = t; end
                if choices(t)<2, tsCE(end+1) = t; end
            end
            
            % Collating useable single trials
            LLR_full = [LLR_full; LLRin(ts,:)];
            LLRsum_full = [LLRsum_full; nansum(LLRin(ts,:),2)];
            LLRend_full = [LLRend_full; LLRin(ts,end)];
            LPR_full = [LPR_full; LPRout(ts,end)];
            LPR_full_i = [LPR_full_i; LPRout_i(ts,end)];
            surprise_full = [surprise_full; surprise(ts,:)];
            surprise_full_i = [surprise_full_i; surprise_i(ts,:)];
            psi_full = [psi_full; psi(ts,:)];
            psi_full_i = [psi_full_i; psi_i(ts,:)];
            choices_full = [choices_full; choices(ts)];
            
            fdist_full = [fdist_full; Behav(tsCE,1)-1];
            acc_full = [acc_full; acc(tsCE)];
            choices_all = [choices_all; choices(tsCE)];
            accB(subj,s,b) = mean(acc(tsCE));
            nswitch = [nswitch; nansum(pswitch(tsCE,:),2)];
            
            fdistc=[]; LPRc=[];
            for t = 1:length(choices)
                if choices(t)<2
                    nsmps = length(find(~isnan(LLRin(t,:))));
                    LPR_4acc_full(end+1,1) = LPRout_i(t,nsmps);
                    LPR_fit_4acc_full(end+1,1) = LPRout(t,nsmps);
                    LLRsum_4acc_full(end+1,1) = sum(LLRin(t,1:nsmps));
                    LLRend_4acc_full(end+1,1) = LLRin(t,nsmps);
                    
                    fdistc(end+1,1) = Behav(t,1)-1;
                    LPRc(end+1,1) = LPRout_i(t,nsmps);
                end
            end
            nblock = nblock+1; clear stimIn Behav gen
            
            % Calculate block-wise normative accuracy (for debugging)
            fdistc(fdistc==0) = -1;
            cacc = sign(LPRc);
            accBi(subj,s,b) = length(find(cacc-fdistc==0))./length(cacc)*100;
        end
    end
    GA_nblock(subj,1) = nblock;
    
    % Store subj-specific s.d. of LPRs
    GA_lpr_sd(subj,1) = std(LPR_full_i);
    
    % Running choice regressions with lasso regression
    [rB,STATS] = lassoglm(nanzscore([LLR_full LLR_full(:,2:end).*nanzscore(surprise_full(:,2:end)) LLR_full(:,2:end).*zscore(-abs(psi_full(:,2:end)))]),[choices_full],'binomial','Lambda',0.002,'CV',10);
    rGA_ssB_regU(subj,1:maxsamps) = rB(1:maxsamps);
    rGA_ssB_surpU(subj,:) = rB(maxsamps+1:end-length(2:maxsamps));
    rGA_ssB_uncertU(subj,:) = rB(end-length(2:maxsamps)+1:end);
    
    % Running fitted model choice regressions with lasso regression
    CPm = 0.5+0.5.*erf(LPR_full./GA_noise(subj,1));
    
    [rB,STATS] = lassoglm(nanzscore([LLR_full LLR_full(:,2:end).*nanzscore(surprise_full(:,2:end)) LLR_full(:,2:end).*zscore(-abs(psi_full(:,2:end)))]),[CPm],'binomial','Lambda',0.002,'CV',10);
    rGA_ssB_regUm(subj,1:maxsamps) = rB(1:maxsamps);
    rGA_ssB_surpUm(subj,:) = rB(maxsamps+1:end-length(2:maxsamps));
    rGA_ssB_uncertUm(subj,:) = rB(end-length(2:maxsamps)+1:end);
    
    % Compute second minus first half of regular & surprise kernels
    rGA_regU_halfdiff(subj,1) = mean(rGA_ssB_regU(subj,(maxsamps/2)+1:end))-mean(rGA_ssB_regU(subj,1:(maxsamps/2)));
    rGA_surpU_halfdiff(subj,1) = mean(rGA_ssB_surpU(subj,(maxsamps/2)+1:end))-mean(rGA_ssB_surpU(subj,1:(maxsamps/2)-1));
    
    % Compute time-averages of kernels
    rGA_ssB_regU_av(subj,1) = mean(rGA_ssB_regU(subj,CPPsmps));
    rGA_ssB_surpU_av(subj,1) = mean(rGA_ssB_surpU(subj,CPPsmps-1));
    rGA_ssB_uncU_av(subj,1) = mean(rGA_ssB_uncertU(subj,UNCsmps-1));
    
    rGA_ssB_regU_av_early(subj,1) = mean(rGA_ssB_regU(subj,1:min(CPPsmps)-1));
    rGA_ssB_surpU_av_early(subj,1) = mean(rGA_ssB_surpU(subj,1:min(CPPsmps)-2));
    
    % Calculate accuracy of idealized strategies
    fdistc = fdist_full; fdistc(fdistc==0) = -1;
    
    cacc = sign(LPR_4acc_full);
    accNorm(subj,1) = length(find(cacc-fdistc==0))./length(cacc)*100;
    consistNorm(subj,1) = length(find(cacc-choices_all==0))./length(cacc)*100;
    
    cp = 0.5+0.5.*erf(LPR_fit_4acc_full./GA_noise(subj,1));
    cacc = cp; cacc(find(fdistc==-1)) = 1-cacc(find(fdistc==-1));
    accNormFit(subj,1) = mean(cacc)*100;
    
    cacc = sign(LLRend_4acc_full);
    accLast(subj,1) = length(find(cacc-fdistc==0))./length(cacc)*100;
    
    cacc = sign(LLRsum_4acc_full);
    accSum(subj,1) = length(find(cacc-fdistc==0))./length(cacc)*100;                        
    
    % Computing per-subject accuracy (objective, relative to true gen. state)
    GAacc = [GAacc; nansum(acc_full)/length(find(~isnan(acc_full)))*100];
    GAacc_sw = [GAacc_sw; nansum(acc_full(nswitch>0))/length(find(~isnan(acc_full(nswitch>0))))*100];
    % Computing per-subject choice consistency relative to noiseless ideal observer
    cacc = sign(LPR_4acc_full);
    cacc(cacc==-1) = 0;
    GAacc_i = [GAacc_i; length(find(cacc(~isnan(choices_all))-choices_all(~isnan(choices_all))==0))/length(cacc(~isnan(choices_all)))];
    GAacc_i_sw = [GAacc_i_sw; length(find(cacc(~isnan(choices_all) & nswitch>0)-choices_all(~isnan(choices_all) & nswitch>0)==0))/length(cacc(~isnan(choices_all) & nswitch>0))];
    % Computing per-subject choice consistency relative to normative fit
    cacc = sign(LPR_fit_4acc_full);
    cacc(cacc==-1) = 0;
    GAacc_fit = [GAacc_fit; length(find(cacc(~isnan(choices_all))-choices_all(~isnan(choices_all))==0))/length(cacc(~isnan(choices_all)))];
    % Computing per-subject choice consistency relative to perfect accumulator
    cacc = sign(LLRsum_4acc_full);
    cacc(cacc==-1) = 0;
    GAacc_perf = [GAacc_perf; length(find(cacc(~isnan(choices_all))-choices_all(~isnan(choices_all))==0))/length(cacc(~isnan(choices_all)))];
    % Computing per-subject choice consistency relative to last sample heuristic
    cacc = sign(LLRend_4acc_full);
    cacc(cacc==-1) = 0;
    GAacc_last = [GAacc_last; length(find(cacc(~isnan(choices_all))-choices_all(~isnan(choices_all))==0))/length(cacc(~isnan(choices_all)))];
end


%%% Plot sorted single-subject accuracies
[GAacc,i] = sort(GAacc);
acc_norm_binary = accNormFit(i);
acc_fsmp_binary = accLast(i);
acc_perfacc_binary = accSum(i);
acc_ideal_nn_binary = accNorm(i);

fs = 7.5;
lw = 1.0;
axlw = 1.0;
dotsize = 8;

fig_w = 18; fig_h = 5.2;

h =  findobj('type','figure');
cfig = length(h)+1;
f=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','painters'), hold on

barw = 0.8; ylims = [0.6 0.821].*100;

f=fill([1-barw/2 1+barw/2 1+barw/2 1-barw/2],[ylims(1)+0.003 ylims(1)+0.003 mean(acc_fsmp_binary) mean(acc_fsmp_binary)],[1 0 1]); set(f,'EdgeColor',[1 1 1]);
plot([1 1],[mean(acc_fsmp_binary)-(std(acc_fsmp_binary)/sqrt(length(acc_fsmp_binary))) mean(acc_fsmp_binary)+(std(acc_fsmp_binary)/sqrt(length(acc_fsmp_binary)))],'k')

f=fill([2-barw/2 2+barw/2 2+barw/2 2-barw/2],[ylims(1)+0.003 ylims(1)+0.003 mean(acc_perfacc_binary) mean(acc_perfacc_binary)],[0 1 0]); set(f,'EdgeColor',[1 1 1]);
plot([2 2],[mean(acc_perfacc_binary)-(std(acc_perfacc_binary)/sqrt(length(acc_perfacc_binary))) mean(acc_perfacc_binary)+(std(acc_perfacc_binary)/sqrt(length(acc_perfacc_binary)))],'k')

for s = 3:length(acc_norm_binary)+2
    f=fill([s-barw/2 s+barw/2 s+barw/2 s-barw/2],[ylims(1)+0.003 ylims(1)+0.003 GAacc(s-2) GAacc(s-2)],[0.75 0.75 0.75]); set(f,'EdgeColor',[1 1 1]);
    S=scatter(s,acc_norm_binary(s-2),dotsize); set(S,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
end

f=fill([s+1-barw/2 s+1+barw/2 s+1+barw/2 s+1-barw/2],[ylims(1)+0.003 ylims(1)+0.003 mean(acc_ideal_nn_binary) mean(acc_ideal_nn_binary)],[0 0 1]); set(f,'EdgeColor',[1 1 1]);
plot([s+1 s+1],[mean(acc_ideal_nn_binary)-(std(acc_ideal_nn_binary)/sqrt(length(acc_ideal_nn_binary))) mean(acc_ideal_nn_binary)+(std(acc_ideal_nn_binary)/sqrt(length(acc_ideal_nn_binary)))],'k')

ylabel('Accuracy (%)','FontSize',10)
xlim([0.3 s+1.7]), ylim(ylims)
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','XTick',[1:s+1],'XTickLabel',[],'YTick',[0.5:0.1:0.9].*100)



% Plotting single-subject kernels
nperm = 10000;
fprintf('\nRunning cluster-based tests on kernels...')
[sig_ts1,sig_ts_uc1,cp1,~] = cluster_permWS_fast(cat(3,rGA_ssB_regU,zeros(size(rGA_ssB_regU))),nperm,0.01,0.05);
[sig_ts2,sig_ts_uc2,cp2,~] = cluster_permWS_fast(cat(3,rGA_ssB_surpU,zeros(size(rGA_ssB_surpU))),nperm,0.01,0.05);
[sig_ts3,sig_ts_uc3,cp3,~] = cluster_permWS_fast(cat(3,rGA_ssB_uncertU,zeros(size(rGA_ssB_uncertU))),nperm,0.01,0.05);

fs=8; axlw=0.5;

fig_w = 16; fig_h = 4.2;

h =  findobj('type','figure');
cfig = length(h)+1;
f=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

subplot(1,3,1), hold on
plot(1:maxsamps,rGA_ssB_regU,'Color',[0.6 0.6 0.6],'LineWidth',0.5)
gaM = mean(rGA_ssB_regUm,1); seM = std(rGA_ssB_regUm,[],1)/sqrt(size(rGA_ssB_regUm,1));
f=fill([1:maxsamps maxsamps:-1:1],[gaM+seM fliplr(gaM-seM)],[1 0 0]); set(f,'EdgeColor','none')
shadedErrorBar(1:maxsamps,mean(rGA_ssB_regU,1),std(rGA_ssB_regU,[],1)./sqrt(size(rGA_ssB_regU,1)),{'Color',[0 0 0],'LineWidth',1.25},0)
plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.75 0.75 0.75],'LineWidth',0.75)
plot(1:maxsamps,sig_ts1.*-0.35,'k','LineWidth',3)
ylabel('Weight on choice (a.u.)','FontSize',10)
xlim([0.5 maxsamps+0.5]), ylim([-0.5 5.4])
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','XTick',[2 6 10],'YTick',[0 2 4 6])

subplot(1,3,2), hold on
plot(2:maxsamps,rGA_ssB_surpU,'Color',[0.6 0.6 0.6],'LineWidth',0.5)
gaM = mean(rGA_ssB_surpUm,1); seM = std(rGA_ssB_surpUm,[],1)/sqrt(size(rGA_ssB_surpUm,1));
f=fill([2:maxsamps maxsamps:-1:2],[gaM+seM fliplr(gaM-seM)],[1 0 0]); set(f,'EdgeColor','none')
shadedErrorBar(2:maxsamps,mean(rGA_ssB_surpU,1),std(rGA_ssB_surpU,[],1)./sqrt(size(rGA_ssB_surpU,1)),{'Color',[0 0 0],'LineWidth',1.25},0)
plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.75 0.75 0.75],'LineWidth',0.75)
plot(2:maxsamps,sig_ts2.*-0.7,'k','LineWidth',3)
xlabel('Sample position','FontSize',10), 
xlim([0.5 maxsamps+0.5]), ylim([-0.75 1.25])
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','XTick',[2 6 10],'YTick',[-0.7 0 0.7 1.4])

subplot(1,3,3), hold on
plot(2:maxsamps,rGA_ssB_uncertU,'Color',[0.6 0.6 0.6],'LineWidth',0.5)
gaM = mean(rGA_ssB_uncertUm,1); seM = std(rGA_ssB_uncertUm,[],1)/sqrt(size(rGA_ssB_uncertUm,1));
f=fill([2:maxsamps maxsamps:-1:2],[gaM+seM fliplr(gaM-seM)],[1 0 0]); set(f,'EdgeColor','none')
shadedErrorBar(2:maxsamps,mean(rGA_ssB_uncertU,1),std(rGA_ssB_uncertU,[],1)./sqrt(size(rGA_ssB_uncertU,1)),{'Color',[0 0 0],'LineWidth',1.25},0)
plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.75 0.75 0.75],'LineWidth',0.75)
plot(2:maxsamps,sig_ts3.*-0.7,'k','LineWidth',3)
xlabel('Sample position','FontSize',10), 
xlim([0.5 maxsamps+0.5]), ylim([-0.75 1.25])
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','XTick',[2 6 10],'YTick',[-0.7 0 0.7 1.4])


% Plotting single-subject kernels WITH MODEL PREDICTIONS
fs = 8;
lw = 1.0;
dotsize = 12;

fig_w = 16; fig_h = 4;

h =  findobj('type','figure');
cfig = length(h)+1;
f=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

gaO = mean(rGA_ssB_regU,1); seO = std(rGA_ssB_regU,[],1)/sqrt(size(rGA_ssB_regU,1));
gaM = mean(rGA_ssB_regUm,1); seM = std(rGA_ssB_regUm,[],1)/sqrt(size(rGA_ssB_regUm,1));
ylims = [-0.2 2.7];
subplot(1,3,1), hold on
plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.7 0.7 0.7],'LineWidth',lw)
f=fill([1:maxsamps maxsamps:-1:1],[gaM+seM fliplr(gaM-seM)],[1 0 0]); set(f,'EdgeColor','none')
for p = 1:length(gaO)
    line([p p],[gaO(p)+seO(p) gaO(p)-seO(p)],'Color',[0 0 0],'LineWidth',lw)
end
plot(1:maxsamps,gaO,'Color',[0 0 0],'LineWidth',lw)
ylabel('Weight on choice (a.u.)','FontSize',10)
xlim([-0.5 maxsamps+0.5]), ylim(ylims)
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','XTick',[2:2:12],'YTick',[0 1.2 2.4])

gaO = mean(rGA_ssB_surpU,1); seO = std(rGA_ssB_surpU,[],1)/sqrt(size(rGA_ssB_surpU,1));
gaM = mean(rGA_ssB_surpUm,1); seM = std(rGA_ssB_surpUm,[],1)/sqrt(size(rGA_ssB_surpUm,1));
ylims = [-0.2 0.43];
subplot(1,3,2), hold on
plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.7 0.7 0.7],'LineWidth',lw)
f=fill([2:maxsamps maxsamps:-1:2],[gaM+seM fliplr(gaM-seM)],[1 0 0]); set(f,'EdgeColor','none')
for p = 1:length(gaO)
    line([p+1 p+1],[gaO(p)+seO(p) gaO(p)-seO(p)],'Color',[0 0 0],'LineWidth',lw)
end
plot(2:maxsamps,gaO,'Color',[0 0 0],'LineWidth',lw)
xlabel('Sample position','FontSize',10)
xlim([0.5 maxsamps+0.5]), ylim(ylims)
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','XTick',[2:2:12],'YTick',[-0.2 0 0.2 0.4])

gaO = mean(rGA_ssB_uncertU,1); seO = std(rGA_ssB_uncertU,[],1)/sqrt(size(rGA_ssB_uncertU,1));
gaM = mean(rGA_ssB_uncertUm,1); seM = std(rGA_ssB_uncertUm,[],1)/sqrt(size(rGA_ssB_uncertUm,1));
subplot(1,3,3), hold on
plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.7 0.7 0.7],'LineWidth',lw)
f=fill([2:maxsamps maxsamps:-1:2],[gaM+seM fliplr(gaM-seM)],[1 0 0]); set(f,'EdgeColor','none')
for p = 1:length(gaO)
    line([p+1 p+1],[gaO(p)+seO(p) gaO(p)-seO(p)],'Color',[0 0 0],'LineWidth',lw)
end
plot(2:maxsamps,gaO,'Color',[0 0 0],'LineWidth',lw)
xlim([0.5 maxsamps+0.5]), ylim(ylims)
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','XTick',[2:2:12],'YTick',[-0.2 0 0.2 0.4])


% Plotting accuracies relative to idealized strategies
figure,

jitbnd = [-0.3 0.3];
scatsize = 9;

hold on

rnds = dist_scatters(accLast,0.02);  % get x jitters
S=scatter(ones(size(accLast,1),1) + diff(jitbnd).*rnds + jitbnd(1),accLast,scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colsF(1,:))
plot(jitbnd+1,[median(accLast) median(accLast)],'Color',colsF(1,:),'LineStyle','-','LineWidth',3)

rnds = dist_scatters(accSum,0.02);  % get x jitters
S=scatter(ones(size(accSum,1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),accSum,scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colsF(2,:))
plot(jitbnd+2,[median(accSum) median(accSum)],'Color',colsF(2,:),'LineStyle','-','LineWidth',3)

rnds = dist_scatters(GAacc,0.02);  % get x jitters
S=scatter(ones(size(GAacc,1),1).*3 + diff(jitbnd).*rnds + jitbnd(1),GAacc,scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colsF(3,:))
plot(jitbnd+3,[median(GAacc) median(GAacc)],'Color',colsF(3,:),'LineStyle','-','LineWidth',3)

rnds = dist_scatters(accNorm,0.02);  % get x jitters
S=scatter(ones(size(accNorm,1),1).*4 + diff(jitbnd).*rnds + jitbnd(1),accNorm,scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colsF(4,:))
plot(jitbnd+4,[median(accNorm) median(accNorm)],'Color',colsF(4,:),'LineStyle','-','LineWidth',3)

set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','ydir','normal','XTick',[1:4],'XTickLabel',{'Last Smp.','Perfect','Data','Ideal'})
ylabel('Accuracy'), xlim([0.5 4.5])



% Get CAPE scores
[CAPE,names] = get_CAPE_new(allsubj);
CAPE(:,strcmp(names,'P')) = (CAPE(:,strcmp(names,'P'))+20)./20;  % rescaling scores so that choice scale is 1-4 (rather than coded 0-3) and total score is average (rather than sum) over all items
CAPE(:,strcmp(names,'N')) = (CAPE(:,strcmp(names,'N'))+14)./14;
CAPE(:,strcmp(names,'D')) = (CAPE(:,strcmp(names,'D'))+8)./8;


% Plot CAPE score histograms with quintile cutoffs
nbins = 20;
hist_col = [0.4 0.4 0.4];

fig_w = 16; fig_h = 4;

h =  findobj('type','figure');
cfig = length(h)+1;
f=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

subplot(1,3,1), hold on
h1 = histogram(CAPE(:,strcmp(names,'P')),nbins);
cutoffs = prctile(CAPE(:,strcmp(names,'P')),[20 80]);
plot(repmat(cutoffs,2,1),[0 0 ;10 10],'LineStyle','--','Color','k')
text(0.9,10.25,num2str(round(mean(CAPE(CAPE(:,strcmp(names,'P'))<cutoffs(1),strcmp(names,'P'))),2)),'FontName','Arial','FontSize',fs)
text(1.86,10.25,num2str(round(mean(CAPE(CAPE(:,strcmp(names,'P'))>=cutoffs(2),strcmp(names,'P'))),2)),'FontName','Arial','FontSize',fs)
set(h1,'FaceColor',[1 1 1])
set(h1,'EdgeColor',hist_col)
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off')
ylabel('No. subjects'), xlabel('P-score'), ylim([0 11])

subplot(1,3,2), hold on
h1 = histogram(CAPE(:,strcmp(names,'N')),nbins);
cutoffs = prctile(CAPE(:,strcmp(names,'N')),[20 80]);
plot(repmat(cutoffs,2,1),[0 0 ;18 18],'LineStyle','--','Color','k')
text(0.9,18.5,num2str(round(mean(CAPE(CAPE(:,strcmp(names,'N'))<cutoffs(1),strcmp(names,'N'))),2)),'FontName','Arial','FontSize',fs)
text(2.35,18.5,num2str(round(mean(CAPE(CAPE(:,strcmp(names,'N'))>=cutoffs(2),strcmp(names,'N'))),2)),'FontName','Arial','FontSize',fs)
set(h1,'FaceColor',[1 1 1])
set(h1,'EdgeColor',hist_col)
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off')
ylabel('No. subjects'), xlabel('N-score'), ylim([0 19.4])

subplot(1,3,3), hold on
h1 = histogram(CAPE(:,strcmp(names,'D')),nbins);
cutoffs = prctile(CAPE(:,strcmp(names,'D')),[20 80]);
plot(repmat(cutoffs,2,1),[0 0 ;24 24],'LineStyle','--','Color','k')
text(0.98,25,num2str(round(mean(CAPE(CAPE(:,strcmp(names,'D'))<cutoffs(1),strcmp(names,'D'))),2)),'FontName','Arial','FontSize',fs)
text(2.4,25,num2str(round(mean(CAPE(CAPE(:,strcmp(names,'D'))>=cutoffs(2),strcmp(names,'D'))),2)),'FontName','Arial','FontSize',fs)
set(h1,'FaceColor',[1 1 1])
set(h1,'EdgeColor',hist_col)
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off')
ylabel('No. subjects'), xlabel('D-score'), ylim([0 26.5])


% Median split by CAPE scores and compare kernels
nperm=10000; clustalpha=0.05; alpha=0.05;
colL = [0.5 0.5 0.5];
colH = [0 0 0];
colI = [0 0 0.8];

transbars = 0;

figure
for sc = 1:length(names)

    cuts = prctile(CAPE(:,sc),[20 80]);
    gL = find(CAPE(:,sc)<=cuts(1));
    gH = find(CAPE(:,sc)>=cuts(2));

    % plot LLR kernels
    subplot(length(names),3,1+3*(sc-1)), hold on
    sig_ts=[]; for t=1:maxsamps, sig_ts(t)=permutationTest(rGA_ssB_regU(gL,t), rGA_ssB_regU(gH,t), nperm); end
    shadedErrorBar(1:maxsamps,mean(rGA_ssB_regU(gL,:),1),std_err(rGA_ssB_regU(gL,:),1),{'Color',colL},transbars)
    shadedErrorBar(1:maxsamps,mean(rGA_ssB_regU(gH,:),1),std_err(rGA_ssB_regU(gH,:),1),{'Color',colH},transbars)
    s2=scatter(find(sig_ts<=0.05),ones(1,length(find(sig_ts<0.05))).*0,20); set(s2,'MarkerEdgeColor',[0.8 0 0])
    plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.75 0.75 0.75],'LineWidth',0.75)
    xlim([0.5 maxsamps+0.5])
    xlabel('Sample position'), ylabel('beta (a.u.)'), set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[2 6 10])
    if sc==1, title('LLR'), end
    
    % plot surprise kernels
    subplot(length(names),3,2+3*(sc-1)), hold on
    sig_ts=[]; for t=1:maxsamps-1, sig_ts(t)=permutationTest(rGA_ssB_surpU(gL,t), rGA_ssB_surpU(gH,t), nperm); end
    shadedErrorBar(2:maxsamps,mean(rGA_ssB_surpU(gL,:),1),std_err(rGA_ssB_surpU(gL,:),1),{'Color',colL},transbars)
    shadedErrorBar(2:maxsamps,mean(rGA_ssB_surpU(gH,:),1),std_err(rGA_ssB_surpU(gH,:),1),{'Color',colH},transbars)
    s2=scatter(find(sig_ts<=0.05)+1,ones(1,length(find(sig_ts<0.05))).*0,20); set(s2,'MarkerEdgeColor',[0.8 0 0])
    plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.75 0.75 0.75],'LineWidth',0.75)
    xlim([0.5 maxsamps+0.5])
    xlabel('Sample position'), ylabel('beta (a.u.)'), set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[2 6 10])
    if sc==1, title('CPP'), end
    
    % plot uncertainty kernels
    subplot(length(names),3,3+3*(sc-1)), hold on
    sig_ts=[]; for t=1:maxsamps-1, sig_ts(t)=permutationTest(rGA_ssB_uncertU(gL,t), rGA_ssB_uncertU(gH,t), nperm); end
    shadedErrorBar(2:maxsamps,mean(rGA_ssB_uncertU(gL,:),1),std_err(rGA_ssB_uncertU(gL,:),1),{'Color',colL},transbars)
    shadedErrorBar(2:maxsamps,mean(rGA_ssB_uncertU(gH,:),1),std_err(rGA_ssB_uncertU(gH,:),1),{'Color',colH},transbars)
    s2=scatter(find(sig_ts<=0.05)+1,ones(1,length(find(sig_ts<0.05))).*0,20); set(s2,'MarkerEdgeColor',[0.8 0 0])
    plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.75 0.75 0.75],'LineWidth',0.75)
    xlim([0.5 maxsamps+0.5])
    xlabel('Sample position'), ylabel('beta (a.u.)'), set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[2 6 10])
    if sc==1, title('Uncertainty'), end
end


% plot CPP + UNC kernels
cuts = prctile(CAPE(:,strcmp(names,'P')),[20 80]);
gL = find(CAPE(:,strcmp(names,'P'))<=cuts(1));
gH = find(CAPE(:,strcmp(names,'P'))>=cuts(2));

corrtype = 'spearman';
    
figure,
subplot(2,5,1), hold on
sig_ts=[]; for t=1:maxsamps-1, sig_ts(t)=permutationTest(rGA_ssB_surpU(gL,t), rGA_ssB_surpU(gH,t), nperm); end
f1=fill([CPPsmps(1) CPPsmps(1) CPPsmps(end) CPPsmps(end)],[-0.24 0.4 0.4 -0.24],[0.9 0.9 0.9]); set(f1,'EdgeColor',[1 1 1])
shadedErrorBar(2:maxsamps,mean(rGA_ssB_surpU(gL,:),1),std_err(rGA_ssB_surpU(gL,:),1),{'Color',colL},transbars)
shadedErrorBar(2:maxsamps,mean(rGA_ssB_surpU(gH,:),1),std_err(rGA_ssB_surpU(gH,:),1),{'Color',colH},transbars)
s2=scatter(find(sig_ts<=0.05)+1,ones(1,length(find(sig_ts<0.05))).*0,20); set(s2,'MarkerEdgeColor',[0.8 0 0])
plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.75 0.75 0.75],'LineWidth',0.75)
xlim([0.5 maxsamps+0.5]), ylim([-0.24 0.4])
xlabel('Sample position'), ylabel('Beta weight (a.u.)'), set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[2 6 10])
title('CPP modulation')

subplot(2,5,2), hold on
sig_ts=[]; for t=1:maxsamps-1, sig_ts(t)=permutationTest(rGA_ssB_uncertU(gL,t), rGA_ssB_uncertU(gH,t), nperm); end
f1=fill([UNCsmps(1) UNCsmps(1) UNCsmps(end) UNCsmps(end)],[-0.24 0.4 0.4 -0.24],[0.9 0.9 0.9]); set(f1,'EdgeColor',[1 1 1])
shadedErrorBar(2:maxsamps,mean(rGA_ssB_uncertU(gL,:),1),std_err(rGA_ssB_uncertU(gL,:),1),{'Color',colL},transbars)
shadedErrorBar(2:maxsamps,mean(rGA_ssB_uncertU(gH,:),1),std_err(rGA_ssB_uncertU(gH,:),1),{'Color',colH},transbars)
s2=scatter(find(sig_ts<=0.05)+1,ones(1,length(find(sig_ts<0.05))).*0,20); set(s2,'MarkerEdgeColor',[0.8 0 0])
plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.75 0.75 0.75],'LineWidth',0.75)
xlim([0.5 maxsamps+0.5]), ylim([-0.24 0.4])
xlabel('Sample position'), ylabel('Beta weight (a.u.)'), set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[2 6 10])
title('UNC modulation')

% avg CPP modulation
subplot(2,5,3), hold on
rnds = dist_scatters(rGA_ssB_surpU_av(gL),0.1);  % get x jitters
S=scatter(ones(size(rGA_ssB_surpU_av(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),rGA_ssB_surpU_av(gL),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
plot(jitbnd+1,[mean(rGA_ssB_surpU_av(gL)) mean(rGA_ssB_surpU_av(gL))],'LineWidth',2,'Color',colL)

rnds = dist_scatters(rGA_ssB_surpU_av(gH),0.1);  % get x jitters
S=scatter(ones(size(rGA_ssB_surpU_av(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),rGA_ssB_surpU_av(gH),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
plot(jitbnd+2,[mean(rGA_ssB_surpU_av(gH)) mean(rGA_ssB_surpU_av(gH))],'LineWidth',2,'Color',colH)

p = permutationTest(rGA_ssB_surpU_av(gL), rGA_ssB_surpU_av(gH), nperm);
text(1.1,0.46,['p=',num2str(round(p,3))],'FontSize',7)

set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
ylabel('Mean beta weight (a.u.)'), xlim([0.5 2.5])
title('CPP kernel magnitude')

subplot(2,5,8), hold on
[rho,p]=corr(CAPE(:,strcmp(names,'P')),rGA_ssB_surpU_av,'type',corrtype,'rows','pairwise');
scatter(CAPE(:,strcmp(names,'P')),rGA_ssB_surpU_av,scatsize), ylabel('CPP modulation'), xlabel(names{1})
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

% avg UNC modulation
subplot(2,5,4), hold on
rnds = dist_scatters(rGA_ssB_uncU_av(gL),0.1);  % get x jitters
S=scatter(ones(size(rGA_ssB_uncU_av(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),rGA_ssB_uncU_av(gL),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
plot(jitbnd+1,[mean(rGA_ssB_uncU_av(gL)) mean(rGA_ssB_uncU_av(gL))],'LineWidth',2,'Color',colL)

rnds = dist_scatters(rGA_ssB_uncU_av(gH),0.1);  % get x jitters
S=scatter(ones(size(rGA_ssB_uncU_av(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),rGA_ssB_uncU_av(gH),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
plot(jitbnd+2,[mean(rGA_ssB_uncU_av(gH)) mean(rGA_ssB_uncU_av(gH))],'LineWidth',2,'Color',colH)

p = permutationTest(rGA_ssB_uncU_av(gL), rGA_ssB_uncU_av(gH), nperm);
text(1.1,0.08,['p=',num2str(round(p,3))],'FontSize',7)

set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
ylabel('Mean beta weight (a.u.)'), xlim([0.5 2.5])
title('UNC kernel magnitude')

subplot(2,5,9), hold on
[rho,p]=corr(CAPE(:,strcmp(names,'P')),rGA_ssB_uncU_av,'type',corrtype,'rows','pairwise');
scatter(CAPE(:,strcmp(names,'P')),rGA_ssB_uncU_av,scatsize), ylabel('UNC modulation'), xlabel(names{1})
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

% plot average magnitude of CPP+UNC modulations
subplot(2,5,5), hold on
rnds = dist_scatters(rGA_ssB_surpU_av(gL)+rGA_ssB_uncU_av(gL),0.1);  % get x jitters
S=scatter(ones(size(rGA_ssB_surpU_av(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),rGA_ssB_surpU_av(gL)+rGA_ssB_uncU_av(gL),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
plot(jitbnd+1,[mean(rGA_ssB_surpU_av(gL)+rGA_ssB_uncU_av(gL)) mean(rGA_ssB_surpU_av(gL)+rGA_ssB_uncU_av(gL))],'LineWidth',2,'Color',colL)

rnds = dist_scatters(rGA_ssB_surpU_av(gH)+rGA_ssB_uncU_av(gH),0.1);  % get x jitters
S=scatter(ones(size(rGA_ssB_surpU_av(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),rGA_ssB_surpU_av(gH)+rGA_ssB_uncU_av(gH),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
plot(jitbnd+2,[mean(rGA_ssB_surpU_av(gH)+rGA_ssB_uncU_av(gH)) mean(rGA_ssB_surpU_av(gH)+rGA_ssB_uncU_av(gH))],'LineWidth',2,'Color',colH)

p = permutationTest(rGA_ssB_surpU_av(gL)+rGA_ssB_uncU_av(gL), rGA_ssB_surpU_av(gH)+rGA_ssB_uncU_av(gH), nperm);
text(1.1,0.5,['p=',num2str(round(p,3))],'FontSize',7)

set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
ylabel('Integrated beta weight (a.u.)'), xlim([0.5 2.5])
title('CPP+UNC kernel magnitude')

subplot(2,5,10), hold on
[rho,p]=corr(CAPE(:,strcmp(names,'P')),rGA_ssB_surpU_av+rGA_ssB_uncU_av,'type',corrtype,'rows','pairwise');
scatter(CAPE(:,strcmp(names,'P')),rGA_ssB_surpU_av+rGA_ssB_uncU_av,scatsize), ylabel('CPP+UNC modulation'), xlabel(names{1})
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')



%%% Iterate through different data splits (by P-scale only) and plot key metrics
splits = [49.999 50.001;...
          33.33 66.66;...
          25 75;...
          20 80];
      
figure,
for i = 1:size(splits,1)
    cuts = prctile(CAPE(:,strcmp(names,'P')),[splits(i,1) splits(i,2)]);
    gL = find(CAPE(:,strcmp(names,'P'))<=cuts(1));
    gH = find(CAPE(:,strcmp(names,'P'))>=cuts(2));
    
    % plot P-score histogram with cutoffs
    subplot(size(splits,1),9,1+9*(i-1)), hold on
    h1 = histogram(CAPE(:,strcmp(names,'P')),nbins);
    plot(repmat(cuts,2,1),[0 0 ;10 10],'LineStyle','--','Color','k')
    text(0.9,10.25,num2str(round(mean(CAPE(CAPE(:,strcmp(names,'P'))<cuts(1),strcmp(names,'P'))),2)),'FontName','Arial','FontSize',fs)
    text(1.86,10.25,num2str(round(mean(CAPE(CAPE(:,strcmp(names,'P'))>=cuts(2),strcmp(names,'P'))),2)),'FontName','Arial','FontSize',fs)
    set(h1,'FaceColor',[1 1 1])
    set(h1,'EdgeColor',hist_col)
    set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off')
    ylabel('No. subjects'), xlabel('P-score'), ylim([0 11])
    
    % plot LLR kernels
    subplot(size(splits,1),9,2+9*(i-1)), hold on
    sig_ts=[]; for t=1:maxsamps, sig_ts(t)=permutationTest(rGA_ssB_regU(gL,t), rGA_ssB_regU(gH,t), nperm); end
    shadedErrorBar(1:maxsamps,mean(rGA_ssB_regU(gL,:),1),std_err(rGA_ssB_regU(gL,:),1),{'Color',colL},transbars)
    shadedErrorBar(1:maxsamps,mean(rGA_ssB_regU(gH,:),1),std_err(rGA_ssB_regU(gH,:),1),{'Color',colH},transbars)
    s2=scatter(find(sig_ts<=0.05),ones(1,length(find(sig_ts<0.05))).*0,20); set(s2,'MarkerEdgeColor',[0.8 0 0])
    plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.75 0.75 0.75],'LineWidth',0.75)
    xlim([0.5 maxsamps+0.5])
    xlabel('Sample position'), ylabel('Beta weight (a.u.)'), set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[2 6 10])
    if i==1, title('Evidence weighting'), end
    
    % plot LLR kernel 2nd minus 1st halves
    subplot(size(splits,1),9,3+9*(i-1)), hold on
    rnds = dist_scatters(rGA_regU_halfdiff(gL),0.1);  % get x jitters
    S=scatter(ones(size(rGA_regU_halfdiff(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),rGA_regU_halfdiff(gL),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
    plot(jitbnd+1,[mean(rGA_regU_halfdiff(gL)) mean(rGA_regU_halfdiff(gL))],'LineWidth',2,'Color',colL)

    rnds = dist_scatters(rGA_regU_halfdiff(gH),0.1);  % get x jitters
    S=scatter(ones(size(rGA_regU_halfdiff(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),rGA_regU_halfdiff(gH),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
    plot(jitbnd+2,[mean(rGA_regU_halfdiff(gH)) mean(rGA_regU_halfdiff(gH))],'LineWidth',2,'Color',colH)
    
    p = permutationTest(rGA_regU_halfdiff(gL), rGA_regU_halfdiff(gH), nperm);
    text(1.1,1.8,['p=',num2str(round(p,3))],'FontSize',7)
    
    set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
    ylabel('delta beta weight (a.u.)'), xlim([0.5 2.5])
    if i==1, title('Kernel half difference'), end
    
    % plot surprise kernels
    subplot(size(splits,1),9,4+9*(i-1)), hold on
    sig_ts=[]; for t=1:maxsamps-1, sig_ts(t)=permutationTest(rGA_ssB_surpU(gL,t), rGA_ssB_surpU(gH,t), nperm); end
    shadedErrorBar(2:maxsamps,mean(rGA_ssB_surpU(gL,:),1),std_err(rGA_ssB_surpU(gL,:),1),{'Color',colL},transbars)
    shadedErrorBar(2:maxsamps,mean(rGA_ssB_surpU(gH,:),1),std_err(rGA_ssB_surpU(gH,:),1),{'Color',colH},transbars)
    s2=scatter(find(sig_ts<=0.05)+1,ones(1,length(find(sig_ts<0.05))).*0,20); set(s2,'MarkerEdgeColor',[0.8 0 0])
    plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.75 0.75 0.75],'LineWidth',0.75)
    xlim([0.5 maxsamps+0.5])
    xlabel('Sample position'), ylabel('Beta weight (a.u.)'), set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[2 6 10])
    if i==1, title('CPP modulation'), end
        
    % plot average magnitude of CPP+UNC kernels
    subplot(size(splits,1),9,5+9*(i-1)), hold on
    rnds = dist_scatters(rGA_ssB_surpU_av(gL)+rGA_ssB_uncU_av(gL),0.1);  % get x jitters
    S=scatter(ones(size(rGA_ssB_surpU_av(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),rGA_ssB_surpU_av(gL)+rGA_ssB_uncU_av(gL),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
    plot(jitbnd+1,[mean(rGA_ssB_surpU_av(gL)+rGA_ssB_uncU_av(gL)) mean(rGA_ssB_surpU_av(gL)+rGA_ssB_uncU_av(gL))],'LineWidth',2,'Color',colL)

    rnds = dist_scatters(rGA_ssB_surpU_av(gH)+rGA_ssB_uncU_av(gH),0.1);  % get x jitters
    S=scatter(ones(size(rGA_ssB_surpU_av(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),rGA_ssB_surpU_av(gH)+rGA_ssB_uncU_av(gH),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
    plot(jitbnd+2,[mean(rGA_ssB_surpU_av(gH)+rGA_ssB_uncU_av(gH)) mean(rGA_ssB_surpU_av(gH)+rGA_ssB_uncU_av(gH))],'LineWidth',2,'Color',colH)
    
    p = permutationTest(rGA_ssB_surpU_av(gL)+rGA_ssB_uncU_av(gL), rGA_ssB_surpU_av(gH)+rGA_ssB_uncU_av(gH), nperm);
    text(1.1,0.5,['p=',num2str(round(p,3))],'FontSize',7)
    
    set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
    ylabel('Integrated beta weight (a.u.)'), xlim([0.5 2.5])
    if i==1, title('CPP+UNC kernel magnitude, sample positions 5-10'), end
    
    % plot choice accuracy
    subplot(size(splits,1),9,6+9*(i-1)), hold on
    rnds = dist_scatters(GAacc_sw(gL,1)./100,0.2);  % get x jitters
    S=scatter(ones(size(GAacc_sw(gL,1)./100,1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),GAacc_sw(gL,1)./100,scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
    plot(jitbnd+1,[mean(GAacc_sw(gL)) mean(GAacc_sw(gL))]./100,'LineWidth',2,'Color',colL)

    rnds = dist_scatters(GAacc_sw(gH,1)./100,0.2);  % get x jitters
    S=scatter(ones(size(GAacc_sw(gH,1)./100,1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),GAacc_sw(gH,1)./100,scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
    plot(jitbnd+2,[mean(GAacc_sw(gH)) mean(GAacc_sw(gH))]./100,'LineWidth',2,'Color',colH)
    
    p = permutationTest(GAacc_sw(gL,1)./100, GAacc_sw(gH,1)./100, nperm);
    text(1.1,0.78,['p=',num2str(round(p,3))],'FontSize',7)
    
    set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
    ylabel('Proportion correct'), xlim([0.5 2.5])
    if i==1, title('Accuracy, CP trials only'), end
    
    % plot choice constistency with ideal observer
    subplot(size(splits,1),9,7+9*(i-1)), hold on
    rnds = dist_scatters(GAacc_i_sw(gL,1),0.2);  % get x jitters
    S=scatter(ones(size(GAacc_i_sw(gL,1),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),GAacc_i_sw(gL,1),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
    plot(jitbnd+1,[mean(GAacc_i_sw(gL)) mean(GAacc_i_sw(gL))],'LineWidth',2,'Color',colL)

    rnds = dist_scatters(GAacc_i_sw(gH,1),0.2);  % get x jitters
    S=scatter(ones(size(GAacc_i_sw(gH,1),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),GAacc_i_sw(gH,1),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
    plot(jitbnd+2,[mean(GAacc_i_sw(gH)) mean(GAacc_i_sw(gH))],'LineWidth',2,'Color',colH)
    
    p = permutationTest(GAacc_i_sw(gL,1), GAacc_i_sw(gH,1), nperm);
    text(1.1,0.88,['p=',num2str(round(p,3))],'FontSize',7)
    
    set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
    ylabel('Proportion consistent'), xlim([0.5 2.5])
    if i==1, title('Consistency with ideal observer, CP trials only'), end
    
    % plot H
    subplot(size(splits,1),9,8+9*(i-1)), hold on
    rnds = dist_scatters(GA_H(gL,1),0.2);  % get x jitters
    S=scatter(ones(size(GA_H(gL,1),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),GA_H(gL,1),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
    plot(jitbnd+1,[mean(GA_H(gL)) mean(GA_H(gL))],'LineWidth',2,'Color',colL)

    rnds = dist_scatters(GA_H(gH,1),0.2);  % get x jitters
    S=scatter(ones(size(GA_H(gH,1),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),GA_H(gH,1),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
    plot(jitbnd+2,[mean(GA_H(gH)) mean(GA_H(gH))],'LineWidth',2,'Color',colH)
    
    p = permutationTest(GA_H(gL,1), GA_H(gH,1), nperm);
    text(1.1,0.29,['p=',num2str(round(p,3))],'FontSize',7)
    
    set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
    ylabel('H'), xlim([0.5 2.5])
    if i==1, title('H'), end
    
    % plot noise
    subplot(size(splits,1),9,9+9*(i-1)), hold on
    rnds = dist_scatters(GA_noise(gL,1),0.2);  % get x jitters
    S=scatter(ones(size(GA_noise(gL,1),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),GA_noise(gL,1),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
    plot(jitbnd+1,[mean(GA_noise(gL)) mean(GA_noise(gL))],'LineWidth',2,'Color',colL)

    rnds = dist_scatters(GA_noise(gH,1),0.2);  % get x jitters
    S=scatter(ones(size(GA_noise(gH,1),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),GA_noise(gH,1),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
    plot(jitbnd+2,[mean(GA_noise(gH)) mean(GA_noise(gH))],'LineWidth',2,'Color',colH)
    
    p = permutationTest(GA_noise(gL,1), GA_noise(gH,1), nperm);
    text(1.1,3.8,['p=',num2str(round(p,3))],'FontSize',7)
    
    set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
    ylabel('noise'), xlim([0.5 2.5])
    if i==1, title('noise'), end
end



for i = 1:size(splits,1)
    cuts = prctile(CAPE(:,strcmp(names,'P')),[splits(i,1) splits(i,2)]);
    gL = find(CAPE(:,strcmp(names,'P'))<=cuts(1));
    gH = find(CAPE(:,strcmp(names,'P'))>=cuts(2));
    
    % kernel difference
    kdiffBIN(i,1:3) = [mean(rGA_regU_halfdiff(gL)) mean(rGA_regU_halfdiff(gH)) permutationTest(rGA_regU_halfdiff(gL), rGA_regU_halfdiff(gH), nperm)];
    % surprise modulation
    surpmodBIN(i,1:3) = [mean(rGA_ssB_surpU_av(gL)) mean(rGA_ssB_surpU_av(gH)) permutationTest(rGA_ssB_surpU_av(gL), rGA_ssB_surpU_av(gH), nperm)];
    % combined modulation
    combmodBIN(i,1:3) = [mean(rGA_ssB_surpU_av(gL)+rGA_ssB_uncU_av(gL)) mean(rGA_ssB_surpU_av(gH)+rGA_ssB_uncU_av(gH)) permutationTest(rGA_ssB_surpU_av(gL)+rGA_ssB_uncU_av(gL), rGA_ssB_surpU_av(gH)+rGA_ssB_uncU_av(gH), nperm)];
    % choice accuracy
    accswBIN(i,1:3) = [mean(GAacc_sw(gL,1)./100) mean(GAacc_sw(gH,1)./100) permutationTest(GAacc_sw(gL,1)./100, GAacc_sw(gH,1)./100, nperm)];
    % choice consistency
    acciswBIN(i,1:3) = [mean(GAacc_i_sw(gL,1)) mean(GAacc_i_sw(gH,1)) permutationTest(GAacc_i_sw(gL,1), GAacc_i_sw(gH,1), nperm)];
    % H
    hBIN(i,1:3) = [mean(GA_H(gL,1)) mean(GA_H(gH,1)) permutationTest(GA_H(gL,1), GA_H(gH,1), nperm)];
    % noise
    noiseBIN(i,1:3) = [mean(GA_noise(gL,1)) mean(GA_noise(gH,1)) permutationTest(GA_noise(gL,1), GA_noise(gH,1), nperm)];
end

markersize = 23;
colSig = [1 0 0];
colTrend = [1 0.6 0.6];

figure,
subplot(6,1,1), hold on
plot(1:4,kdiffBIN(:,1),'Color',colL)
plot(1:4,kdiffBIN(:,2),'Color',colH)
for i = 1:size(splits,1)
    if kdiffBIN(i,3)<=0.05
        S = scatter([i i],kdiffBIN(i,1:2),markersize); set(S,'MarkerEdgeColor',colSig,'MarkerFaceColor',colSig)
    elseif kdiffBIN(i,3)<=0.1
        S = scatter([i i],kdiffBIN(i,1:2),markersize); set(S,'MarkerEdgeColor',colTrend,'MarkerFaceColor',colTrend)
    else
        S = scatter([i],kdiffBIN(i,1),markersize); set(S,'MarkerEdgeColor',colL,'MarkerFaceColor',colL)
        S = scatter([i],kdiffBIN(i,2),markersize); set(S,'MarkerEdgeColor',colH,'MarkerFaceColor',colH)
    end
end
set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:4],'XTickLabel',{'Median split','Tertile','Quartile','Quintile'},'XTickLabelRotation',0)
ylabel('delta beta weight (a.u.)'), xlim([0.5 4.5])
title('Kernel half difference')

subplot(6,1,2), hold on
plot(1:4,combmodBIN(:,1),'Color',colL)
plot(1:4,combmodBIN(:,2),'Color',colH)
for i = 1:size(splits,1)
    if combmodBIN(i,3)<=0.05
        S = scatter([i i],combmodBIN(i,1:2),markersize); set(S,'MarkerEdgeColor',colSig,'MarkerFaceColor',colSig)
    elseif combmodBIN(i,3)<=0.1
        S = scatter([i i],combmodBIN(i,1:2),markersize); set(S,'MarkerEdgeColor',colTrend,'MarkerFaceColor',colTrend)
    else
        S = scatter([i],combmodBIN(i,1),markersize); set(S,'MarkerEdgeColor',colL,'MarkerFaceColor',colL)
        S = scatter([i],combmodBIN(i,2),markersize); set(S,'MarkerEdgeColor',colH,'MarkerFaceColor',colH)
    end
end
set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:4],'XTickLabel',{'Median split','Tertile','Quartile','Quintile'},'XTickLabelRotation',0)
ylabel('Integrated beta weight (a.u.)'), xlim([0.5 4.5])
title('CPP+UNC kernel magnitude, sample positions 5-10')

subplot(6,1,3), hold on
plot(1:4,accswBIN(:,1),'Color',colL)
plot(1:4,accswBIN(:,2),'Color',colH)
for i = 1:size(splits,1)
    if accswBIN(i,3)<=0.05
        S = scatter([i i],accswBIN(i,1:2),markersize); set(S,'MarkerEdgeColor',colSig,'MarkerFaceColor',colSig)
    elseif accswBIN(i,3)<=0.1
        S = scatter([i i],accswBIN(i,1:2),markersize); set(S,'MarkerEdgeColor',colTrend,'MarkerFaceColor',colTrend)
    else
        S = scatter([i],accswBIN(i,1),markersize); set(S,'MarkerEdgeColor',colL,'MarkerFaceColor',colL)
        S = scatter([i],accswBIN(i,2),markersize); set(S,'MarkerEdgeColor',colH,'MarkerFaceColor',colH)
    end
end
set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:4],'XTickLabel',{'Median split','Tertile','Quartile','Quintile'},'XTickLabelRotation',0)
ylabel('Proportion correct'), xlim([0.5 4.5])
title('Accuracy, CP trials only')

subplot(6,1,4), hold on
plot(1:4,acciswBIN(:,1),'Color',colL)
plot(1:4,acciswBIN(:,2),'Color',colH)
for i = 1:size(splits,1)
    if acciswBIN(i,3)<=0.05
        S = scatter([i i],acciswBIN(i,1:2),markersize); set(S,'MarkerEdgeColor',colSig,'MarkerFaceColor',colSig)
    elseif acciswBIN(i,3)<=0.1
        S = scatter([i i],acciswBIN(i,1:2),markersize); set(S,'MarkerEdgeColor',colTrend,'MarkerFaceColor',colTrend)
    else
        S = scatter([i],acciswBIN(i,1),markersize); set(S,'MarkerEdgeColor',colL,'MarkerFaceColor',colL)
        S = scatter([i],acciswBIN(i,2),markersize); set(S,'MarkerEdgeColor',colH,'MarkerFaceColor',colH)
    end
end
set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:4],'XTickLabel',{'Median split','Tertile','Quartile','Quintile'},'XTickLabelRotation',0)
ylabel('Proportion consistent'), xlim([0.5 4.5])
title('Consistency with ideal observer, CP trials only')

subplot(6,1,5), hold on
plot(1:4,hBIN(:,1),'Color',colL)
plot(1:4,hBIN(:,2),'Color',colH)
for i = 1:size(splits,1)
    if hBIN(i,3)<=0.05
        S = scatter([i i],hBIN(i,1:2),markersize); set(S,'MarkerEdgeColor',colSig,'MarkerFaceColor',colSig)
    elseif hBIN(i,3)<=0.1
        S = scatter([i i],hBIN(i,1:2),markersize); set(S,'MarkerEdgeColor',colTrend,'MarkerFaceColor',colTrend)
    else
        S = scatter([i],hBIN(i,1),markersize); set(S,'MarkerEdgeColor',colL,'MarkerFaceColor',colL)
        S = scatter([i],hBIN(i,2),markersize); set(S,'MarkerEdgeColor',colH,'MarkerFaceColor',colH)
    end
end
set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:4],'XTickLabel',{'Median split','Tertile','Quartile','Quintile'},'XTickLabelRotation',0)
ylabel('H'), xlim([0.5 4.5])
title('Subjective hazard rate')

subplot(6,1,6), hold on
plot(1:4,noiseBIN(:,1),'Color',colL)
plot(1:4,noiseBIN(:,2),'Color',colH)
for i = 1:size(splits,1)
    if noiseBIN(i,3)<=0.05
        S = scatter([i i],noiseBIN(i,1:2),markersize); set(S,'MarkerEdgeColor',colSig,'MarkerFaceColor',colSig)
    elseif noiseBIN(i,3)<=0.1
        S = scatter([i i],noiseBIN(i,1:2),markersize); set(S,'MarkerEdgeColor',colTrend,'MarkerFaceColor',colTrend)
    else
        S = scatter([i],noiseBIN(i,1),markersize); set(S,'MarkerEdgeColor',colL,'MarkerFaceColor',colL)
        S = scatter([i],noiseBIN(i,2),markersize); set(S,'MarkerEdgeColor',colH,'MarkerFaceColor',colH)
    end
end
set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:4],'XTickLabel',{'Median split','Tertile','Quartile','Quintile'},'XTickLabelRotation',0)
ylabel('Noise'), xlim([0.5 4.5])
title('Decision noise')



% various plots relating P-scores to behaviour and model parameters, both quintile split and continuous correlations
corrtype = 'spearman';

figure

splits = [20 80]; i=1;

cuts = prctile(CAPE(:,strcmp(names,'P')),[splits(i,1) splits(i,2)]);
gL = find(CAPE(:,strcmp(names,'P'))<=cuts(1));
gH = find(CAPE(:,strcmp(names,'P'))>=cuts(2));

% accuracy
subplot(2,6,1), hold on
rnds = dist_scatters(GAacc(gL,1),0.2);  % get x jitters
S=scatter(ones(size(GAacc(gL,1),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),GAacc(gL,1),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
plot(jitbnd+1,[mean(GAacc(gL)) mean(GAacc(gL))],'LineWidth',2,'Color',colL)

rnds = dist_scatters(GAacc(gH,1),0.2);  % get x jitters
S=scatter(ones(size(GAacc(gH,1),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),GAacc(gH,1),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
plot(jitbnd+2,[mean(GAacc(gH)) mean(GAacc(gH))],'LineWidth',2,'Color',colH)

p = permutationTest(GAacc(gL,1), GAacc(gH,1), nperm);
text(1.1,79,['p=',num2str(round(p,4))],'FontSize',7)

set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
ylabel('H'), xlim([0.5 2.5])
if i==1, title('Accuracy'), end

subplot(2,6,7), hold on
[rho,p]=corr(CAPE(:,strcmp(names,'P')),GAacc,'type',corrtype,'rows','pairwise');
scatter(CAPE(:,strcmp(names,'P')),GAacc,scatsize), ylabel('Accuracy'), xlabel(names{1})
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

% accuracy, CP trials only
subplot(2,6,2), hold on
rnds = dist_scatters(GAacc_sw(gL,1),0.2);  % get x jitters
S=scatter(ones(size(GAacc_sw(gL,1),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),GAacc_sw(gL,1),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
plot(jitbnd+1,[mean(GAacc_sw(gL)) mean(GAacc_sw(gL))],'LineWidth',2,'Color',colL)

rnds = dist_scatters(GAacc_sw(gH,1),0.2);  % get x jitters
S=scatter(ones(size(GAacc_sw(gH,1),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),GAacc_sw(gH,1),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
plot(jitbnd+2,[mean(GAacc_sw(gH)) mean(GAacc_sw(gH))],'LineWidth',2,'Color',colH)

p = permutationTest(GAacc_sw(gL,1), GAacc_sw(gH,1), nperm);
text(1.1,74,['p=',num2str(round(p,4))],'FontSize',7)

set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
ylabel('H'), xlim([0.5 2.5])
if i==1, title('Accuracy, CP-only'), end

subplot(2,6,8), hold on
[rho,p]=corr(CAPE(:,strcmp(names,'P')),GAacc_sw,'type',corrtype,'rows','pairwise');
scatter(CAPE(:,strcmp(names,'P')),GAacc_sw,scatsize), ylabel('Accuracy, CP-only'), xlabel(names{1})
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

% accuracy relative to ideal
subplot(2,6,3), hold on
rnds = dist_scatters(GAacc_i(gL,1),0.2);  % get x jitters
S=scatter(ones(size(GAacc_i(gL,1),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),GAacc_i(gL,1),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
plot(jitbnd+1,[mean(GAacc_i(gL)) mean(GAacc_i(gL))],'LineWidth',2,'Color',colL)

rnds = dist_scatters(GAacc_i(gH,1),0.2);  % get x jitters
S=scatter(ones(size(GAacc_i(gH,1),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),GAacc_i(gH,1),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
plot(jitbnd+2,[mean(GAacc_i(gH)) mean(GAacc_i(gH))],'LineWidth',2,'Color',colH)

p = permutationTest(GAacc_i(gL,1), GAacc_i(gH,1), nperm);
text(1.1,0.91,['p=',num2str(round(p,4))],'FontSize',7)

set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
ylabel('H'), xlim([0.5 2.5])
if i==1, title('Accuracy rel. to ideal'), end

subplot(2,6,9), hold on
[rho,p]=corr(CAPE(:,strcmp(names,'P')),GAacc_i,'type',corrtype,'rows','pairwise');
scatter(CAPE(:,strcmp(names,'P')),GAacc_i,scatsize), ylabel('Accuracy rel. to ideal'), xlabel(names{1})
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

% accuracy relative to ideal, CP trials only
subplot(2,6,4), hold on
rnds = dist_scatters(GAacc_i_sw(gL,1),0.2);  % get x jitters
S=scatter(ones(size(GAacc_i_sw(gL,1),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),GAacc_i_sw(gL,1),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
plot(jitbnd+1,[mean(GAacc_i_sw(gL)) mean(GAacc_i_sw(gL))],'LineWidth',2,'Color',colL)

rnds = dist_scatters(GAacc_i_sw(gH,1),0.2);  % get x jitters
S=scatter(ones(size(GAacc_i_sw(gH,1),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),GAacc_i_sw(gH,1),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
plot(jitbnd+2,[mean(GAacc_i_sw(gH)) mean(GAacc_i_sw(gH))],'LineWidth',2,'Color',colH)

p = permutationTest(GAacc_i_sw(gL,1), GAacc_i_sw(gH,1), nperm);
text(1.1,0.89,['p=',num2str(round(p,4))],'FontSize',7)

set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
ylabel('H'), xlim([0.5 2.5])
if i==1, title('Accuracy rel. to ideal, CP-only'), end

subplot(2,6,10), hold on
[rho,p]=corr(CAPE(:,strcmp(names,'P')),GAacc_i_sw,'type',corrtype,'rows','pairwise');
scatter(CAPE(:,strcmp(names,'P')),GAacc_i_sw,scatsize), ylabel('Accuracy rel. to ideal, CP-only'), xlabel(names{1})
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

% H
subplot(2,6,5), hold on
rnds = dist_scatters(GA_H(gL,1),0.2);  % get x jitters
S=scatter(ones(size(GA_H(gL,1),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),GA_H(gL,1),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
plot(jitbnd+1,[mean(GA_H(gL)) mean(GA_H(gL))],'LineWidth',2,'Color',colL)

rnds = dist_scatters(GA_H(gH,1),0.2);  % get x jitters
S=scatter(ones(size(GA_H(gH,1),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),GA_H(gH,1),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
plot(jitbnd+2,[mean(GA_H(gH)) mean(GA_H(gH))],'LineWidth',2,'Color',colH)

p = permutationTest(GA_H(gL,1), GA_H(gH,1), nperm);
text(1.1,0.29,['p=',num2str(round(p,4))],'FontSize',7)

set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
ylabel('H'), xlim([0.5 2.5])
if i==1, title('H'), end

subplot(2,6,11), hold on
[rho,p]=corr(CAPE(:,strcmp(names,'P')),GA_H,'type',corrtype,'rows','pairwise');
scatter(CAPE(:,strcmp(names,'P')),GA_H,scatsize), ylabel('H'), xlabel(names{1})
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

% plot noise
subplot(2,6,6), hold on
rnds = dist_scatters(GA_noise(gL,1),0.2);  % get x jitters
S=scatter(ones(size(GA_noise(gL,1),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),GA_noise(gL,1),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
plot(jitbnd+1,[mean(GA_noise(gL)) mean(GA_noise(gL))],'LineWidth',2,'Color',colL)

rnds = dist_scatters(GA_noise(gH,1),0.2);  % get x jitters
S=scatter(ones(size(GA_noise(gH,1),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),GA_noise(gH,1),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
plot(jitbnd+2,[mean(GA_noise(gH)) mean(GA_noise(gH))],'LineWidth',2,'Color',colH)

p = permutationTest(GA_noise(gL,1), GA_noise(gH,1), nperm);
text(1.1,3.8,['p=',num2str(round(p,4))],'FontSize',7)

set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
ylabel('noise'), xlim([0.5 2.5])
if i==1, title('noise'), end

subplot(2,6,12), hold on
[rho,p]=corr(CAPE(:,strcmp(names,'P')),GA_noise,'type',corrtype,'rows','pairwise');
scatter(CAPE(:,strcmp(names,'P')),GA_noise,scatsize), ylabel('Noise'), xlabel(names{1})
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')




% various plots relating P-scores (both quintile split and continuous correlations) to kernel measures
figure,

% plot LLR kernels
subplot(3,7,1), hold on
sig_ts=[]; for t=1:maxsamps, sig_ts(t)=permutationTest(rGA_ssB_regU(gL,t), rGA_ssB_regU(gH,t), nperm); end
shadedErrorBar(1:maxsamps,mean(rGA_ssB_regU(gL,:),1),std_err(rGA_ssB_regU(gL,:),1),{'Color',colL},transbars)
shadedErrorBar(1:maxsamps,mean(rGA_ssB_regU(gH,:),1),std_err(rGA_ssB_regU(gH,:),1),{'Color',colH},transbars)
s2=scatter(find(sig_ts<=0.05),ones(1,length(find(sig_ts<0.05))).*0,20); set(s2,'MarkerEdgeColor',[0.8 0 0])
plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.75 0.75 0.75],'LineWidth',0.75)
xlim([0.5 maxsamps+0.5])
xlabel('Sample position'), ylabel('beta (a.u.)'), set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off','XTick',[2 6 10])
title('LLR kernel')

% plot LLR kernel avg (LATE)
subplot(3,7,2), hold on

rnds = dist_scatters(rGA_ssB_regU_av(gL),0.1);  % get x jitters
S=scatter(ones(size(rGA_ssB_regU_av(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),rGA_ssB_regU_av(gL),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
plot(jitbnd+1,[mean(rGA_ssB_regU_av(gL)) mean(rGA_ssB_regU_av(gL))],'LineWidth',2,'Color',colL)

rnds = dist_scatters(rGA_ssB_regU_av(gH),0.1);  % get x jitters
S=scatter(ones(size(rGA_ssB_regU_av(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),rGA_ssB_regU_av(gH),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
plot(jitbnd+2,[mean(rGA_ssB_regU_av(gH)) mean(rGA_ssB_regU_av(gH))],'LineWidth',2,'Color',colH)

p = permutationTest(rGA_ssB_regU_av(gL), rGA_ssB_regU_av(gH), nperm);
text(1.1,1.05,['p=',num2str(round(p,4))],'FontSize',7)

set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
ylabel('beta (a.u.)'), xlim([0.5 2.5])
title('Kernel avg.')

subplot(3,7,3), hold on
[rho,p]=corr(CAPE(:,strcmp(names,'P')),rGA_ssB_regU_av,'type',corrtype,'rows','pairwise');
scatter(CAPE(:,strcmp(names,'P')),rGA_ssB_regU_av,scatsize), ylabel('beta (a.u.)'), xlabel(names{1})
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))]), xlim([0.9 2.4])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

% plot LLR kernel avg (EARLY)
subplot(3,7,4), hold on

rnds = dist_scatters(rGA_ssB_regU_av_early(gL),0.1);  % get x jitters
S=scatter(ones(size(rGA_ssB_regU_av_early(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),rGA_ssB_regU_av_early(gL),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
plot(jitbnd+1,[mean(rGA_ssB_regU_av_early(gL)) mean(rGA_ssB_regU_av_early(gL))],'LineWidth',2,'Color',colL)

rnds = dist_scatters(rGA_ssB_regU_av_early(gH),0.1);  % get x jitters
S=scatter(ones(size(rGA_ssB_regU_av_early(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),rGA_ssB_regU_av_early(gH),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
plot(jitbnd+2,[mean(rGA_ssB_regU_av_early(gH)) mean(rGA_ssB_regU_av_early(gH))],'LineWidth',2,'Color',colH)

p = permutationTest(rGA_ssB_regU_av_early(gL), rGA_ssB_regU_av_early(gH), nperm);
text(1.1,1.05,['p=',num2str(round(p,4))],'FontSize',7)

set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
ylabel('beta (a.u.)'), xlim([0.5 2.5])
title('Kernel avg. (early smps)')

subplot(3,7,5), hold on
[rho,p]=corr(CAPE(:,strcmp(names,'P')),rGA_ssB_regU_av_early,'type',corrtype,'rows','pairwise');
scatter(CAPE(:,strcmp(names,'P')),rGA_ssB_regU_av_early,scatsize), ylabel('beta (a.u.)'), xlabel(names{1})
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))]), xlim([0.9 2.4])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

% plot LLR kernel 2nd minus 1st halves
subplot(3,7,6), hold on

rnds = dist_scatters(rGA_regU_halfdiff(gL),0.1);  % get x jitters
S=scatter(ones(size(rGA_regU_halfdiff(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),rGA_regU_halfdiff(gL),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
plot(jitbnd+1,[mean(rGA_regU_halfdiff(gL)) mean(rGA_regU_halfdiff(gL))],'LineWidth',2,'Color',colL)

rnds = dist_scatters(rGA_regU_halfdiff(gH),0.1);  % get x jitters
S=scatter(ones(size(rGA_regU_halfdiff(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),rGA_regU_halfdiff(gH),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
plot(jitbnd+2,[mean(rGA_regU_halfdiff(gH)) mean(rGA_regU_halfdiff(gH))],'LineWidth',2,'Color',colH)

p = permutationTest(rGA_regU_halfdiff(gL), rGA_regU_halfdiff(gH), nperm);
text(1.1,1.8,['p=',num2str(round(p,4))],'FontSize',7)

set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
ylabel('delta beta (a.u.)'), xlim([0.5 2.5])
title('Kernel half diff.')

subplot(3,7,7), hold on
[rho,p]=corr(rGA_regU_halfdiff,CAPE(:,strcmp(names,'P')),'type',corrtype,'rows','pairwise');
scatter(CAPE(:,strcmp(names,'P')),rGA_regU_halfdiff,scatsize), ylabel('delta beta (a.u.)'), xlabel(names{1})
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))]), xlim([0.9 2.4])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

% plot LLR*CPP kernels
subplot(3,7,8), hold on
sig_ts=[]; for t=1:maxsamps-1, sig_ts(t)=permutationTest(rGA_ssB_surpU(gL,t), rGA_ssB_surpU(gH,t), nperm); end
shadedErrorBar(2:maxsamps,mean(rGA_ssB_surpU(gL,:),1),std_err(rGA_ssB_surpU(gL,:),1),{'Color',colL},transbars)
shadedErrorBar(2:maxsamps,mean(rGA_ssB_surpU(gH,:),1),std_err(rGA_ssB_surpU(gH,:),1),{'Color',colH},transbars)
s2=scatter(find(sig_ts<=0.05)+1,ones(1,length(find(sig_ts<0.05))).*0,20); set(s2,'MarkerEdgeColor',[0.8 0 0])
plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.75 0.75 0.75],'LineWidth',0.75)
xlim([0.5 maxsamps+0.5])
xlabel('Sample position'), ylabel('beta (a.u.)'), set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off','XTick',[2 6 10])
if i==1, title('LLR*CPP kernel'), end

% plot LLR*CPP kernel avg (LATE)
subplot(3,7,9), hold on
rnds = dist_scatters(rGA_ssB_surpU_av(gL)+rGA_ssB_uncU_av(gL),0.1);  % get x jitters
S=scatter(ones(size(rGA_ssB_surpU_av(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),rGA_ssB_surpU_av(gL)+rGA_ssB_uncU_av(gL),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
plot(jitbnd+1,[mean(rGA_ssB_surpU_av(gL)+rGA_ssB_uncU_av(gL)) mean(rGA_ssB_surpU_av(gL)+rGA_ssB_uncU_av(gL))],'LineWidth',2,'Color',colL)

rnds = dist_scatters(rGA_ssB_surpU_av(gH)+rGA_ssB_uncU_av(gH),0.1);  % get x jitters
S=scatter(ones(size(rGA_ssB_surpU_av(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),rGA_ssB_surpU_av(gH)+rGA_ssB_uncU_av(gH),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
plot(jitbnd+2,[mean(rGA_ssB_surpU_av(gH)+rGA_ssB_uncU_av(gH)) mean(rGA_ssB_surpU_av(gH)+rGA_ssB_uncU_av(gH))],'LineWidth',2,'Color',colH)

p = permutationTest(rGA_ssB_surpU_av(gL)+rGA_ssB_uncU_av(gL), rGA_ssB_surpU_av(gH)+rGA_ssB_uncU_av(gH), nperm);
text(1.1,0.3,['p=',num2str(round(p,4))],'FontSize',7)

set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
ylabel('Integrated beta weight (a.u.)'), xlim([0.5 2.5])

subplot(3,7,10), hold on
[rho,p]=corr(rGA_ssB_surpU_av+rGA_ssB_uncU_av,CAPE(:,strcmp(names,'P')),'type',corrtype,'rows','pairwise');
scatter(CAPE(:,strcmp(names,'P')),rGA_ssB_surpU_av+rGA_ssB_uncU_av,scatsize), ylabel('Integrated beta (a.u.)'), xlabel(names{1})
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))]), xlim([0.9 2.4])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

% plot LLR*CPP kernel avg (EARLY)
subplot(3,7,11), hold on

rnds = dist_scatters(rGA_ssB_surpU_av_early(gL),0.1);  % get x jitters
S=scatter(ones(size(rGA_ssB_surpU_av_early(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),rGA_ssB_surpU_av_early(gL),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
plot(jitbnd+1,[mean(rGA_ssB_surpU_av_early(gL)) mean(rGA_ssB_surpU_av_early(gL))],'LineWidth',2,'Color',colL)

rnds = dist_scatters(rGA_ssB_surpU_av_early(gH),0.1);  % get x jitters
S=scatter(ones(size(rGA_ssB_surpU_av_early(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),rGA_ssB_surpU_av_early(gH),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
plot(jitbnd+2,[mean(rGA_ssB_surpU_av_early(gH)) mean(rGA_ssB_surpU_av_early(gH))],'LineWidth',2,'Color',colH)

p = permutationTest(rGA_ssB_surpU_av_early(gL), rGA_ssB_surpU_av_early(gH), nperm);
text(1.1,0.3,['p=',num2str(round(p,4))],'FontSize',7)

set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
ylabel('beta (a.u.)'), xlim([0.5 2.5])

subplot(3,7,12), hold on
[rho,p]=corr(rGA_ssB_surpU_av_early,CAPE(:,strcmp(names,'P')),'type',corrtype,'rows','pairwise');
scatter(CAPE(:,strcmp(names,'P')),rGA_ssB_surpU_av_early,scatsize), ylabel('beta (a.u.)'), xlabel(names{1})
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))]), xlim([0.9 2.4])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

% plot LLR*CPP kernel 2nd minus 1st halves
subplot(3,7,13), hold on

rnds = dist_scatters(rGA_surpU_halfdiff(gL),0.1);  % get x jitters
S=scatter(ones(size(rGA_surpU_halfdiff(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),rGA_surpU_halfdiff(gL),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
plot(jitbnd+1,[mean(rGA_surpU_halfdiff(gL)) mean(rGA_surpU_halfdiff(gL))],'LineWidth',2,'Color',colL)

rnds = dist_scatters(rGA_surpU_halfdiff(gH),0.1);  % get x jitters
S=scatter(ones(size(rGA_surpU_halfdiff(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),rGA_surpU_halfdiff(gH),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
plot(jitbnd+2,[mean(rGA_surpU_halfdiff(gH)) mean(rGA_surpU_halfdiff(gH))],'LineWidth',2,'Color',colH)

p = permutationTest(rGA_surpU_halfdiff(gL), rGA_surpU_halfdiff(gH), nperm);
text(1.1,0.3,['p=',num2str(round(p,4))],'FontSize',7)

set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
ylabel('beta (a.u.)'), xlim([0.5 2.5])

subplot(3,7,14), hold on
[rho,p]=corr(rGA_surpU_halfdiff,CAPE(:,strcmp(names,'P')),'type',corrtype,'rows','pairwise');
scatter(CAPE(:,strcmp(names,'P')),rGA_surpU_halfdiff,scatsize), ylabel('beta (a.u.)'), xlabel(names{1})
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))]), xlim([0.9 2.4])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')



% save some variables for correlation with pupil
allsubj_B = allsubj;   % renaming as will conflict with variable names in pupil scripts otherwise
CAPE_B = CAPE;
save([savepath,'kernel_output.mat'],'allsubj_B','CAPE_B','names','CPPsmps','rGA_ssB_regU','rGA_ssB_surpU','rGA_ssB_uncertU','rGA_regU_halfdiff','rGA_ssB_surpU_av','rGA_ssB_uncU_av')

