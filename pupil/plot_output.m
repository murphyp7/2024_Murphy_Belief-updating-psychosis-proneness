clear, close all

loadpath = '/mnt/homes/home024/pmurphy/Surprise_scz/pupil/5.output/';
behavpath = '/mnt/homes/home028/gmonov/SCZ/Data/decision_making/';
modelpath = '/mnt/homes/home024/pmurphy/Surprise_scz/modelling/';
behavoutpath = '/mnt/homes/home024/pmurphy/Surprise_scz/analysis/';

addpath '/mnt/homes/home024/pmurphy/Surprise_scz/pupil/'
addpath '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/Gen_fun';
addpath(genpath('/mnt/homes/home024/pmurphy/Surprise_scz/'))

model_type = 'H_noise';

subj_exclude = {'t022', 't023', 't042', 't046', 't071', 't077', 't089'};  % subjects to exclude based on task performance/lack of CAPE scores

% Analysis stuff
first_deriv = 1;
ntrial_cutoff = 400;  % minimum # trials for participant to be included

maxsamps = 10;

fullwin = [-2 6.2];  % window for plotting full dilation response aligned to trial onset
sampwin = [0 1.4];  % window for plotting dilation response to indiviudal samples (default = [-0.15 1.5])
fbwin = [-0.5 3.5];    % window for plotting dilation response to feeback

% get list of all files
if first_deriv
    files = dir([loadpath,'*_d1_',model_type,'.mat']);
else files = dir([loadpath,'*_raw_',model_type,'.mat']);
end
for f = 1:length(files), allsubj{f}=files(f).name(1:4); end
allsubj = unique(allsubj);

% Loop through subjects
for subj = 1:length(files)
    
    % load files
    load([loadpath,files(subj).name]);
    
    % concatenate output variables over subjects
    GAdil_full(subj,:) = dil_full_av;
    GAdil_samp(subj,:,:) = dil_samp_av;
    GAdil_fbC(subj,:) = dil_fbC_av;
    GAdil_fbE(subj,:) = dil_fbEav;
    
    GAntrials(subj,1) = ntrials;
    GApart(subj,1) = p_artifact;
    
    GAdil_pCP_B(subj,:,:) = dil_pCP_B;
    GAdil_pCP_Bext(subj,:,:,:) = dil_pCP_Bext;
    
    GAdil_PPI_B(subj,:,:) = dil_PPI_B;
    GAdil_PPI_Bext(subj,:,:) = dil_PPI_Bext;
    GAdil_PPI_B_lasso(subj,:,:) = dil_PPI_B_lasso;
    GAdil_PPI_Bext_lasso(subj,:,:) = dil_PPI_Bext_lasso;
    
    GA_H(subj,1) = mfit.H;
    GA_noise(subj,1) = mfit.noise;
    if strcmp(modeltype,'H_noise_IU')
        GA_IU(subj,1) = mfit.IU;
    end
end

% Load output from kernel analysis script
load([behavoutpath,'kernel_output.mat'])  % 'allsubj_B','CAPE_B','names','CPPsmps','rGA_ssB_surpU','rGA_regU_halfdiff','rGA_ssB_surpU_av'

% Get useable participants
include = find(GAntrials>ntrial_cutoff & ~ismember(allsubj',subj_exclude'));

ext_subj = find(zscore(squeeze(std(mean(GAdil_PPI_Bext,2),[],3)))<5);  % for unregularized choice regressions, excluding those with very extreme fits
include2 = intersect(include,ext_subj);


% Linear projection of single-subject pupil responses onto group-avg
if first_deriv, dil_mean_win = [0.2 1]; else dil_mean_win = [0.5 5.5]; end
dil_proj_win = [0.8 5.4]; dil_ts = find(fulltimes>=dil_proj_win(1) & fulltimes<=dil_proj_win(2));
reg_proj_win = sampwin; reg_ts = find(samptimes>=reg_proj_win(1) & samptimes<=reg_proj_win(2));

for subj = 1:length(files)
    dil_proj(subj,1) = dot(GAdil_full(subj,dil_ts),mean(GAdil_full(include,dil_ts),1))/(norm(mean(GAdil_full(include,dil_ts),1))^2);
    for r = 1:size(GAdil_pCP_Bext,4)
        dil_pCP_Bext_proj(subj,r) = dot(squeeze(mean(GAdil_pCP_Bext(subj,:,reg_ts,r),2)),squeeze(mean(mean(GAdil_pCP_Bext(include,:,reg_ts,r),2),1)))/(norm(squeeze(mean(mean(GAdil_pCP_Bext(include,:,reg_ts,r),2),1)))^2);
    end
    PPI_proj(subj,1) = dot(mean(GAdil_full(subj,reg_ts,:),3),squeeze(mean(mean(GAdil_full(include,reg_ts),3),1)))/(norm(squeeze(mean(mean(GAdil_full(include,reg_ts),3),1)))^2);
end
dil_mean = mean(GAdil_full(:,fulltimes>=dil_mean_win(1) & fulltimes<=dil_mean_win(2)),2);

% Avg/peak measures derived from sample-related responses
CPPsmps = 5:10;   % sample positions over which to average CPP kernels

cpp_win = [0.4 1];
unc_win = [0.2 0.45];
llr_win = [0.35 0.45];

for subj = 1:length(files)
    dil_pCP_Bext_av(subj,1) = squeeze(mean(mean(GAdil_pCP_Bext(subj,CPPsmps-1,samptimes>=cpp_win(1) & samptimes<=cpp_win(2),1),2),3));
    dil_pCP_Bext_av(subj,2) = squeeze(mean(mean(GAdil_pCP_Bext(subj,CPPsmps-1,samptimes>=unc_win(1) & samptimes<=unc_win(2),2),2),3));
    dil_pCP_Bext_av(subj,3) = squeeze(mean(mean(GAdil_pCP_Bext(subj,CPPsmps-1,samptimes>=llr_win(1) & samptimes<=llr_win(2),3),2),3));
    
    dil_pCP_Bext_peak(subj,1) = max(squeeze(mean(GAdil_pCP_Bext(subj,CPPsmps-1,samptimes>=cpp_win(1) & samptimes<=cpp_win(2),1),2)));
    dil_pCP_Bext_peak(subj,2) = max(squeeze(mean(GAdil_pCP_Bext(subj,CPPsmps-1,samptimes>=unc_win(1) & samptimes<=unc_win(2),2),2)));
    dil_pCP_Bext_peak(subj,3) = max(squeeze(mean(GAdil_pCP_Bext(subj,CPPsmps-1,samptimes>=llr_win(1) & samptimes<=llr_win(2),3),2)));
end


% Run cluster-based permutation tests
if first_deriv
    plotwin = [0 1.2];
else plotwin = [0 1.4];
end
plotts = find(samptimes>=plotwin(1) & samptimes<=plotwin(2));
plottimes = samptimes(samptimes>=plotwin(1) & samptimes<=plotwin(2));

nperm = 10000;
fprintf('\nRunning cluster-based tests...')
% [sig_ts1,sig_ts_uc1,cp1,~] = cluster_permWS_fast(cat(3,squeeze(mean(GAdil_pCP_B(include,:,plotts),2)),zeros(length(include),length(plotts))),nperm,0.05,0.05);
% 
% [sig_ts2,sig_ts_uc2,cp2,~] = cluster_permWS_fast(cat(3,squeeze(mean(GAdil_pCP_Bext(include,:,plotts,1),2)),zeros(length(include),length(plotts))),nperm,0.01,0.05);
% [sig_ts3,sig_ts_uc3,cp3,~] = cluster_permWS_fast(cat(3,squeeze(mean(GAdil_pCP_Bext(include,:,plotts,2),2)),zeros(length(include),length(plotts))),nperm,0.01,0.05);
% [sig_ts4,sig_ts_uc4,cp4,~] = cluster_permWS_fast(cat(3,squeeze(mean(GAdil_pCP_Bext(include,:,plotts,3),2)),zeros(length(include),length(plotts))),nperm,0.01,0.05);
% 
% [sig_ts5,sig_ts_uc5,cp5,~] = cluster_permWS_fast(cat(3,squeeze(mean(GAdil_PPI_Bext(include2,:,plotts),2)),zeros(length(include2),length(plotts))),nperm,0.05,0.05);
% [sig_ts6,sig_ts_uc6,cp6,~] = cluster_permWS_fast(cat(3,squeeze(mean(GAdil_PPI_Bext_lasso(include,:,plotts),2)),...
%                                                  zeros(length(include),length(plotts))),nperm,0.05,0.05);
% 
% [sig_ts7,sig_ts_uc7,cp7,~] = cluster_permWS_fast(cat(3,GAdil_full(include,:),zeros(length(include),size(GAdil_full,2))),nperm,0.05,0.05);


% Plotting histogram of trial counts
figure,
subplot(1,2,1), hold on
hist(GAntrials(include),15), xlabel('# trials'), ylabel('# subjects')
subplot(1,2,2), hold on
hist(GApart(include),15), xlabel('fraction of rejected trials'), ylabel('# subjects')


% Plotting average trial-evoked and sample-evoked pupil responses
smarks = repmat(0.4:0.4:0.4*10,2,1);
sampcols = [linspace(0.8,0,size(GAdil_samp,2))' linspace(0.8,0,size(GAdil_samp,2))' linspace(0.8,0,size(GAdil_samp,2))'];

figure,
subplot(1,2,1), hold on,
lx = line(fullwin,[0 0]); set(lx,'LineWidth',0.5,'LineStyle','--','Color',[0 0 0])
% for subj = include
%     plot(fulltimes,GAdil_full(subj,:),'Color',[0.7 0.7 0.7],'LineWidth',0.75)
% end
se = std(GAdil_full(include,:),[],1)./sqrt(length(include));
shadedErrorBar(fulltimes,mean(GAdil_full(include,:),1),se,{'Color',[0.3 0.1 1],'LineWidth',1.5},0);
plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
l = line([smarks],repmat(get(gca, 'ylim')',1,size(smarks,2))); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0.6 0.6 0.6])
% plot(fulltimes,sig_ts7.*-0.1,'k','LineWidth',3)
xlim(fullwin), xlabel('Time relative to trial onset (s)'), ylabel('Pupil dilation (zs^-^1)')
set(gca,'TickDir','out','box','off')
% 
% subplot(1,2,2), hold on,
% lx = line(sampwin,[0 0]); set(lx,'LineWidth',0.5,'LineStyle','--','Color',[0 0 0])
% for samp = 1:size(GAdil_samp,2)
%     plot(samptimes,squeeze(mean(GAdil_samp(include,samp,:),1)),'Color',sampcols(samp,:));
% end
% se = squeeze(std(mean(GAdil_samp(include,:,:),2),[],1))./sqrt(length(include));
% shadedErrorBar(samptimes,squeeze(mean(mean(GAdil_samp(include,:,:),2),1)),se,{'Color',[0.3 0.1 1],'LineWidth',1.5},0);
% plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
% l = line(repmat(0.4:0.4:2.0,2,1),repmat(get(gca, 'ylim')',1,length(0.4:0.4:2.0))); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0.6 0.6 0.6])
% xlim(sampwin), xlabel('Time relative to sample onset (s)'), ylabel('Pupil dilation (zs^-^1)')
% set(gca,'TickDir','out','box','off')


% Plotting CPP/uncertainty encoding
if first_deriv
    avwin = [0.4 0.8];
    avwinU = [0.1 0.5];
    avwinL = [0.3 0.55];
else
    avwin = [0.6 1.4];
    avwinU = [0.3 0.6];
    avwinL = [1.0 1.4];
end

% figure,
% subplot(2,2,1), hold on,
% for samp = 1:size(GAdil_pCP_B,2)
%     p1=plot(plottimes,squeeze(mean(mean(GAdil_pCP_B(include,samp,plotts),2),1)),'Color',sampcols(samp,:),'LineWidth',0.5);
% end
% seS = std(squeeze(mean(GAdil_pCP_B(include,:,plotts),2)),[],1)./sqrt(length(include));
% shadedErrorBar(plottimes,squeeze(mean(mean(GAdil_pCP_B(include,:,plotts),2),1)),seS,{'Color',[1 0 0],'LineWidth',1.5},0);
% plot(plottimes,sig_ts1.*-0.015,'r','LineWidth',3)
% l = line(plotwin,[0 0]); set(l,'LineWidth',0.5,'LineStyle','--')
% xlim([0 plotwin(2)]), ylabel('beta (a.u.'), xlabel('Time relative to sample onset (s)')
% 
% av_dil_surp_B = squeeze(mean(mean(GAdil_pCP_B(include,:,samptimes>=avwin(1) & samptimes<=avwin(2)),3),1));
% seS = std(squeeze(mean(GAdil_pCP_B(include,:,samptimes>=avwin(1) & samptimes<=avwin(2),1),3)),[],1)./sqrt(length(include));
% subplot(2,2,2), hold on,
% shadedErrorBar(2:maxsamps,av_dil_surp_B,seS,{'Color',[1 0 0],'LineStyle','--'},0)
% plot(2:maxsamps,av_dil_surp_B,'Color',[1 0 0],'LineStyle','--');
% l = line([0 maxsamps],[0 0]); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0 0 0])
% for samp = 2:maxsamps
%     S1=scatter(samp,av_dil_surp_B(samp-1)); set(S1,'MarkerEdgeColor',sampcols(samp-1,:));
% end
% xlim([1 maxsamps+0.5]), ylabel('beta (a.u.'), xlabel('Sample position')% , ylim([-0.005 0.055])
% set(gca,'TickDir','out')
% 
% subplot(2,2,3), hold on,
% for samp = 1:size(GAdil_pCP_Bext,2)
%     p1=plot(plottimes,squeeze(mean(mean(GAdil_pCP_Bext(include,:,plotts,1),2),1)),'Color',[1 0 0],'LineWidth',1.5);
%     p2=plot(plottimes,squeeze(mean(mean(GAdil_pCP_Bext(include,:,plotts,2),2),1)),'Color',[0 0 1],'LineWidth',1.5);
%     p3=plot(plottimes,squeeze(mean(mean(GAdil_pCP_Bext(include,:,plotts,3),2),1)),'Color',[0 1 0],'LineWidth',1.5);
% end
% seS = std(squeeze(mean(GAdil_pCP_Bext(include,:,plotts,1),2)),[],1)./sqrt(length(include));
% seL = std(squeeze(mean(GAdil_pCP_Bext(include,:,plotts,2),2)),[],1)./sqrt(length(include));
% seG = std(squeeze(mean(GAdil_pCP_Bext(include,:,plotts,3),2)),[],1)./sqrt(length(include));
% shadedErrorBar(plottimes,squeeze(mean(mean(GAdil_pCP_Bext(include,:,plotts,1),2),1)),seS,{'Color',[1 0 0],'LineWidth',1.5},0);
% shadedErrorBar(plottimes,squeeze(mean(mean(GAdil_pCP_Bext(include,:,plotts,2),2),1)),seL,{'Color',[0 0 1],'LineWidth',1.5},0);
% shadedErrorBar(plottimes,squeeze(mean(mean(GAdil_pCP_Bext(include,:,plotts,3),2),1)),seG,{'Color',[0 1 0],'LineWidth',1.5},0);
% plot(plottimes,sig_ts2.*-0.014,'r','LineWidth',3)
% plot(plottimes,sig_ts3.*-0.016,'b','LineWidth',3)
% plot(plottimes,sig_ts4.*-0.018,'g','LineWidth',3)
% l = line(plotwin,[0 0]); set(l,'LineWidth',0.5,'LineStyle','--')
% xlim([0 plotwin(2)]), ylabel('beta (a.u.)'), xlabel('Time relative to sample onset (s)')
% legend([p1 p2],{'CPP','-|psi|'})
% 
% av_dil_surp_B = squeeze(mean(mean(GAdil_pCP_Bext(include,:,samptimes>=avwin(1) & samptimes<=avwin(2),1),3),1));
% av_dil_deltaL_B = squeeze(mean(mean(GAdil_pCP_Bext(include,:,samptimes>=avwinU(1) & samptimes<=avwinU(2),2),3),1));
% av_dil_absLLR_B = squeeze(mean(mean(GAdil_pCP_Bext(include,:,samptimes>=avwinL(1) & samptimes<=avwinL(2),3),3),1));
% seS = std(squeeze(mean(GAdil_pCP_Bext(include,:,samptimes>=avwin(1) & samptimes<=avwin(2),1),3)),[],1)./sqrt(length(include));
% seL = std(squeeze(mean(GAdil_pCP_Bext(include,:,samptimes>=avwinU(1) & samptimes<=avwinU(2),2),3)),[],1)./sqrt(length(include));
% seG = std(squeeze(mean(GAdil_pCP_Bext(include,:,samptimes>=avwinL(1) & samptimes<=avwinL(2),3),3)),[],1)./sqrt(length(include));
% subplot(2,2,4), hold on,
% shadedErrorBar(2:maxsamps,av_dil_absLLR_B,seG,{'Color',[0 1 0],'LineStyle','--'},0)
% shadedErrorBar(2:maxsamps,av_dil_deltaL_B,seL,{'Color',[0 0 1],'LineStyle','--'},0)
% shadedErrorBar(2:maxsamps,av_dil_surp_B,seS,{'Color',[1 0 0],'LineStyle','--'},0)
% plot(2:maxsamps,av_dil_absLLR_B,'Color',[0 1 0],'LineStyle','--');
% plot(2:maxsamps,av_dil_deltaL_B,'Color',[0 0 1],'LineStyle','--');
% plot(2:maxsamps,av_dil_surp_B,'Color',[1 0 0],'LineStyle','--');
% l = line([0 maxsamps],[0 0]); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0 0 0])
% for samp = 2:maxsamps
%     S2=scatter(samp,av_dil_absLLR_B(samp-1)); set(S2,'MarkerEdgeColor',sampcols(samp-1,:));
%     S=scatter(samp,av_dil_deltaL_B(samp-1)); set(S,'MarkerEdgeColor',sampcols(samp-1,:));
%     S1=scatter(samp,av_dil_surp_B(samp-1)); set(S1,'MarkerEdgeColor',sampcols(samp-1,:));
% end
% xlim([1 maxsamps+0.5]), ylabel('Encoding'), xlabel('Sample position')%, ylim([-0.005 0.055])
% set(gca,'TickDir','out')


% % Plot PPI
% if first_deriv
%     avwin = [0.3 0.45];
% else avwin = [0.9 1.4];
% end
% 
% figure,
% subplot(2,2,1), hold on,
% for samp = 1:size(GAdil_PPI_Bext,2)
%     p1=plot(plottimes,squeeze(mean(mean(GAdil_PPI_Bext(include2,samp,plotts),2),1)),'Color',sampcols(samp,:),'LineWidth',0.5);
% end
% seS = std(squeeze(mean(GAdil_PPI_Bext(include2,:,plotts),2)),[],1)./sqrt(length(include2));
% shadedErrorBar(plottimes,squeeze(mean(mean(GAdil_PPI_Bext(include2,:,plotts),2),1)),seS,{'Color',[1 0 0],'LineWidth',1.5},0);
% plot(plottimes,sig_ts5.*-0.045,'k','LineWidth',3)
% l = line(plotwin,[0 0]); set(l,'LineWidth',0.5,'LineStyle','--')
% xlim([0 plotwin(2)]), ylabel('beta (a.u.'), xlabel('Time relative to sample onset (s)')
% 
% av_dil_surp_B = squeeze(mean(mean(GAdil_PPI_Bext(include2,:,samptimes>=avwin(1) & samptimes<=avwin(2)),3),1));
% seS = std(squeeze(mean(GAdil_PPI_Bext(include2,:,samptimes>=avwin(1) & samptimes<=avwin(2),1),3)),[],1)./sqrt(length(include2));
% subplot(2,2,2), hold on,
% shadedErrorBar(2:maxsamps,av_dil_surp_B,seS,{'Color',[1 0 0],'LineStyle','--'},0)
% plot(2:maxsamps,av_dil_surp_B,'Color',[1 0 0],'LineStyle','--');
% l = line([0 maxsamps],[0 0]); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0 0 0])
% for samp = 2:maxsamps
%     S1=scatter(samp,av_dil_surp_B(samp-1)); set(S1,'MarkerEdgeColor',sampcols(samp-1,:));
% end
% xlim([1 maxsamps+0.5]), ylabel('beta (a.u.'), xlabel('Sample position')% , ylim([-0.005 0.055])
% set(gca,'TickDir','out')
% 
% subplot(2,2,3), hold on,
% for samp = 1:size(GAdil_PPI_Bext_lasso,2)
%     p1=plot(plottimes,squeeze(mean(mean(GAdil_PPI_Bext_lasso(include,samp,plotts),2),1)),'Color',sampcols(samp,:),'LineWidth',0.5);
% end
% seS = std(squeeze(mean(GAdil_PPI_Bext_lasso(include,:,plotts),2)),[],1)./sqrt(length(include));
% shadedErrorBar(plottimes,squeeze(mean(mean(GAdil_PPI_Bext_lasso(include,:,plotts),2),1)),seS,{'Color',[1 0 0],'LineWidth',1.5},0);
% plot(plottimes,sig_ts6.*-0.045,'k','LineWidth',3)
% l = line(plotwin,[0 0]); set(l,'LineWidth',0.5,'LineStyle','--')
% xlim([0 plotwin(2)]), ylabel('beta (a.u.'), xlabel('Time relative to sample onset (s)')
% 
% av_dil_surp_B = squeeze(mean(mean(GAdil_PPI_Bext_lasso(include,:,samptimes>=avwin(1) & samptimes<=avwin(2)),3),1));
% seS = std(squeeze(mean(GAdil_PPI_Bext_lasso(include,:,samptimes>=avwin(1) & samptimes<=avwin(2),1),3)),[],1)./sqrt(length(include));
% subplot(2,2,4), hold on,
% shadedErrorBar(2:maxsamps,av_dil_surp_B,seS,{'Color',[1 0 0],'LineStyle','--'},0)
% plot(2:maxsamps,av_dil_surp_B,'Color',[1 0 0],'LineStyle','--');
% l = line([0 maxsamps],[0 0]); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0 0 0])
% for samp = 2:maxsamps
%     S1=scatter(samp,av_dil_surp_B(samp-1)); set(S1,'MarkerEdgeColor',sampcols(samp-1,:));
% end
% xlim([1 maxsamps+0.5]), ylabel('beta (a.u.'), xlabel('Sample position')% , ylim([-0.005 0.055])
% set(gca,'TickDir','out')




% Correlate pupil measures with kernel measures
fs = 8;
lw = 1.0;
axlw=0.5; 
scatsize = 9;
corrtype = 'spearman';

allsubj_P = allsubj(include)';
include_B = ismember(allsubj_B,allsubj_P);

figure,
subplot(4,2,1), hold on
[rho,p]=corr(dil_mean(include),rGA_regU_halfdiff(include_B),'type',corrtype,'rows','pairwise');
scatter(dil_mean(include),rGA_regU_halfdiff(include_B),scatsize), xlabel('Pupil dilation'), ylabel('Kernel half difference (a.u.)')
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])%, xlim([0.9 2.4])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

subplot(4,2,2), hold on
[rho,p]=corr(dil_mean(include),rGA_ssB_surpU_av(include_B)+rGA_ssB_uncU_av(include_B),'type',corrtype,'rows','pairwise');
scatter(dil_mean(include),rGA_ssB_surpU_av(include_B)+rGA_ssB_uncU_av(include_B),scatsize), xlabel('Pupil dilation'), ylabel('CPP + UNC kernel modulations (a.u.)')
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])%, xlim([0.9 2.4])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

subplot(4,2,3), hold on
[rho,p]=corr(dil_pCP_Bext_av(include,1),rGA_regU_halfdiff(include_B),'type',corrtype,'rows','pairwise');
scatter(dil_pCP_Bext_av(include,1),rGA_regU_halfdiff(include_B),scatsize), xlabel('CPP beta weight (a.u.)'), ylabel('Kernel half difference (a.u.)')
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])%, xlim([0.9 2.4])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

subplot(4,2,4), hold on
[rho,p]=corr(dil_pCP_Bext_av(include,1),rGA_ssB_surpU_av(include_B),'type',corrtype,'rows','pairwise');
scatter(dil_pCP_Bext_av(include,1),rGA_ssB_surpU_av(include_B),scatsize), xlabel('CPP beta weight (a.u.)'), ylabel('CPP kernel magnitude (a.u.)')
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])%, xlim([0.9 2.4])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

subplot(4,2,5), hold on
[rho,p]=corr(dil_pCP_Bext_av(include,2),rGA_regU_halfdiff(include_B),'type',corrtype,'rows','pairwise');
scatter(dil_pCP_Bext_av(include,2),rGA_regU_halfdiff(include_B),scatsize), xlabel('Unc beta weight (a.u.)'), ylabel('Kernel half difference (a.u.)')
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])%, xlim([0.9 2.4])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

subplot(4,2,6), hold on
[rho,p]=corr(dil_pCP_Bext_av(include,2),rGA_ssB_uncU_av(include_B),'type',corrtype,'rows','pairwise');
scatter(dil_pCP_Bext_av(include,2),rGA_ssB_uncU_av(include_B),scatsize), xlabel('Unc beta weight (a.u.)'), ylabel('UNC kernel magnitude (a.u.)')
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])%, xlim([0.9 2.4])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

subplot(4,2,7), hold on
[rho,p]=corr(dil_pCP_Bext_av(include,1)+dil_pCP_Bext_av(include,2),rGA_regU_halfdiff(include_B),'type',corrtype,'rows','pairwise');
scatter(dil_pCP_Bext_av(include,1)+dil_pCP_Bext_av(include,2),rGA_regU_halfdiff(include_B),scatsize), xlabel('CPP+Unc beta weight (a.u.)'), ylabel('Kernel half difference (a.u.)')
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])%, xlim([0.9 2.4])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

subplot(4,2,8), hold on
[rho,p]=corr(dil_pCP_Bext_av(include,1)+dil_pCP_Bext_av(include,2),rGA_ssB_surpU_av(include_B)+rGA_ssB_uncU_av(include_B),'type',corrtype,'rows','pairwise');
scatter(dil_pCP_Bext_av(include,1)+dil_pCP_Bext_av(include,2),rGA_ssB_surpU_av(include_B)+rGA_ssB_uncU_av(include_B),scatsize), xlabel('CPP+Unc beta weight (a.u.)'), ylabel('CPP + UNC kernel modulations (a.u.)')
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])%, xlim([0.9 2.4])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')


% % Correlate with CAPE scores
[CAPE,names] = get_CAPE_new(allsubj');
CAPE(:,strcmp(names,'P')) = (CAPE(:,strcmp(names,'P'))+20)./20;  % rescaling scores so that choice scale is 1-4 (rather than coded 0-3) and total score is average over all items
CAPE(:,strcmp(names,'N')) = (CAPE(:,strcmp(names,'N'))+14)./14;
CAPE(:,strcmp(names,'D')) = (CAPE(:,strcmp(names,'D'))+8)./8;


% 
% figure,
% subplot(3,3,1), hold on
% for t = 1:size(GAdil_full,2)
%     [r1(t),pt1(t)] = corr(GAdil_full(include,t),CAPE(include,1),'type','spearman','rows','complete');
% end
% c_ts1=nan(1,length(fulltimes)); c_ts1(pt1<0.05)=1;
% plot(fulltimes,r1,'Color',[0.3 0.1 1],'LineWidth',1.5)
% plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
% l = line([smarks],repmat(get(gca, 'ylim')',1,size(smarks,2))); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0.6 0.6 0.6])
% plot(fulltimes,c_ts1.*0,'Color',[0 0 0],'LineWidth',3)
% xlim(fullwin), xlabel('Time relative to trial onset (s)'), ylabel('\rho')
% set(gca,'TickDir','out','box','off'), title('Dilation response vs ',names{1})
% 
% subplot(3,3,2), hold on
% for t = 1:size(GAdil_full,2)
%     [r2(t),pt2(t)] = corr(GAdil_full(include,t),CAPE(include,2),'type','spearman','rows','complete');
% end
% c_ts2=nan(1,length(fulltimes)); c_ts2(pt2<0.05)=1;
% plot(fulltimes,r2,'Color',[0.3 0.1 1],'LineWidth',1.5)
% plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
% l = line([smarks],repmat(get(gca, 'ylim')',1,size(smarks,2))); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0.6 0.6 0.6])
% plot(fulltimes,c_ts2.*0,'Color',[0 0 0],'LineWidth',3)
% xlim(fullwin), xlabel('Time relative to trial onset (s)'), ylabel('\rho')
% set(gca,'TickDir','out','box','off'), title('Dilation response vs ',names{2})
% 
% subplot(3,3,3), hold on
% for t = 1:size(GAdil_full,2)
%     [r3(t),pt3(t)] = corr(GAdil_full(include,t),CAPE(include,3),'type','spearman','rows','complete');
% end
% c_ts3=nan(1,length(fulltimes)); c_ts3(pt3<0.05)=1;
% plot(fulltimes,r3,'Color',[0.3 0.1 1],'LineWidth',1.5)
% plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
% l = line([smarks],repmat(get(gca, 'ylim')',1,size(smarks,2))); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0.6 0.6 0.6])
% plot(fulltimes,c_ts3.*0,'Color',[0 0 0],'LineWidth',3)
% xlim(fullwin), xlabel('Time relative to trial onset (s)'), ylabel('\rho')
% set(gca,'TickDir','out','box','off'), title('Dilation response vs ',names{3})
% 
% 
% subplot(3,3,4), hold on
% for t = 1:size(GAdil_pCP_Bext,3)
%     [r4a(t),pt4a(t)] = corr(squeeze(mean(GAdil_pCP_Bext(include,:,t,1),2)),CAPE(include,1),'type','spearman','rows','complete');
%     [r4b(t),pt4b(t)] = corr(squeeze(mean(GAdil_pCP_Bext(include,:,t,2),2)),CAPE(include,1),'type','spearman','rows','complete');
%     [r4c(t),pt4c(t)] = corr(squeeze(mean(GAdil_pCP_Bext(include,:,t,3),2)),CAPE(include,1),'type','spearman','rows','complete');
% end
% c_ts4a=nan(1,length(samptimes)); c_ts4a(pt4a<0.05)=1;
% c_ts4b=nan(1,length(samptimes)); c_ts4b(pt4b<0.05)=1;
% c_ts4c=nan(1,length(samptimes)); c_ts4c(pt4c<0.05)=1;
% plot(samptimes,r4a,'Color',[1 0 0],'LineWidth',1.5)
% plot(samptimes,r4b,'Color',[0 0 1],'LineWidth',1.5)
% plot(samptimes,r4c,'Color',[0 1 0],'LineWidth',1.5)
% plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
% plot(samptimes,c_ts4a.*0.03,'Color',[1 0 0],'LineWidth',3)
% plot(samptimes,c_ts4b.*-0.03,'Color',[0 0 1],'LineWidth',3)
% plot(samptimes,c_ts4c.*-0,'Color',[0 1 0],'LineWidth',3)
% xlim(plotwin), xlabel('Time relative to sample onset (s)'), ylabel('\rho')
% set(gca,'TickDir','out','box','off'), title('Encoding vs ',names{1})
% 
% subplot(3,3,5), hold on
% for t = 1:size(GAdil_pCP_Bext,3)
%     [r4a(t),pt4a(t)] = corr(squeeze(mean(GAdil_pCP_Bext(include,:,t,1),2)),CAPE(include,2),'type','spearman','rows','complete');
%     [r4b(t),pt4b(t)] = corr(squeeze(mean(GAdil_pCP_Bext(include,:,t,2),2)),CAPE(include,2),'type','spearman','rows','complete');
%     [r4c(t),pt4c(t)] = corr(squeeze(mean(GAdil_pCP_Bext(include,:,t,3),2)),CAPE(include,2),'type','spearman','rows','complete');
% end
% c_ts4a=nan(1,length(samptimes)); c_ts4a(pt4a<0.05)=1;
% c_ts4b=nan(1,length(samptimes)); c_ts4b(pt4b<0.05)=1;
% c_ts4c=nan(1,length(samptimes)); c_ts4c(pt4c<0.05)=1;
% plot(samptimes,r4a,'Color',[1 0 0],'LineWidth',1.5)
% plot(samptimes,r4b,'Color',[0 0 1],'LineWidth',1.5)
% plot(samptimes,r4c,'Color',[0 1 0],'LineWidth',1.5)
% plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
% plot(samptimes,c_ts4a.*0.03,'Color',[1 0 0],'LineWidth',3)
% plot(samptimes,c_ts4b.*-0.03,'Color',[0 0 1],'LineWidth',3)
% plot(samptimes,c_ts4c.*-0,'Color',[0 1 0],'LineWidth',3)
% xlim(plotwin), xlabel('Time relative to sample onset (s)'), ylabel('\rho')
% set(gca,'TickDir','out','box','off'), title('Encoding vs ',names{2})
% 
% subplot(3,3,6), hold on
% for t = 1:size(GAdil_pCP_Bext,3)
%     [r4a(t),pt4a(t)] = corr(squeeze(mean(GAdil_pCP_Bext(include,:,t,1),2)),CAPE(include,3),'type','spearman','rows','complete');
%     [r4b(t),pt4b(t)] = corr(squeeze(mean(GAdil_pCP_Bext(include,:,t,2),2)),CAPE(include,3),'type','spearman','rows','complete');
%     [r4c(t),pt4c(t)] = corr(squeeze(mean(GAdil_pCP_Bext(include,:,t,3),2)),CAPE(include,3),'type','spearman','rows','complete');
% end
% c_ts4a=nan(1,length(samptimes)); c_ts4a(pt4a<0.05)=1;
% c_ts4b=nan(1,length(samptimes)); c_ts4b(pt4b<0.05)=1;
% c_ts4c=nan(1,length(samptimes)); c_ts4c(pt4c<0.05)=1;
% plot(samptimes,r4a,'Color',[1 0 0],'LineWidth',1.5)
% plot(samptimes,r4b,'Color',[0 0 1],'LineWidth',1.5)
% plot(samptimes,r4c,'Color',[0 1 0],'LineWidth',1.5)
% plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
% plot(samptimes,c_ts4a.*0.03,'Color',[1 0 0],'LineWidth',3)
% plot(samptimes,c_ts4b.*-0.03,'Color',[0 0 1],'LineWidth',3)
% plot(samptimes,c_ts4c.*-0,'Color',[0 1 0],'LineWidth',3)
% xlim(plotwin), xlabel('Time relative to sample onset (s)'), ylabel('\rho')
% set(gca,'TickDir','out','box','off'), title('Encoding vs ',names{3})
% 
% 
% subplot(3,3,7), hold on
% for t = 1:size(GAdil_PPI_Bext_lasso,3)
%     [r7(t),pt7(t)] = corr(squeeze(mean(GAdil_PPI_Bext_lasso(include,:,t),2)),CAPE(include,1),'type','spearman','rows','complete');
% end
% c_ts7=nan(1,length(samptimes)); c_ts7(pt7<0.05)=1;
% plot(samptimes,r7,'Color',[0.3 0.3 0.3],'LineWidth',1.5)
% plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
% plot(samptimes,c_ts7.*0,'Color',[0 0 0],'LineWidth',3)
% xlim(plotwin), xlabel('Time relative to sample onset (s)'), ylabel('\rho')
% set(gca,'TickDir','out','box','off'), title('PPI vs ',names{1})
% 
% subplot(3,3,8), hold on
% for t = 1:size(GAdil_PPI_Bext_lasso,3)
%     [r8(t),pt8(t)] = corr(squeeze(mean(GAdil_PPI_Bext_lasso(include,:,t),2)),CAPE(include,2),'type','spearman','rows','complete');
% end
% c_ts8=nan(1,length(samptimes)); c_ts8(pt8<0.05)=1;
% plot(samptimes,r8,'Color',[0.3 0.3 0.3],'LineWidth',1.5)
% plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
% plot(samptimes,c_ts8.*0,'Color',[0 0 0],'LineWidth',3)
% xlim(plotwin), xlabel('Time relative to sample onset (s)'), ylabel('\rho')
% set(gca,'TickDir','out','box','off'), title('PPI vs ',names{2})
% 
% subplot(3,3,9), hold on
% for t = 1:size(GAdil_PPI_Bext_lasso,3)
%     [r9(t),pt9(t)] = corr(squeeze(mean(GAdil_PPI_Bext_lasso(include,:,t),2)),CAPE(include,3),'type','spearman','rows','complete');
% end
% c_ts9=nan(1,length(samptimes)); c_ts9(pt9<0.05)=1;
% plot(samptimes,r9,'Color',[0.3 0.3 0.3],'LineWidth',1.5)
% plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
% plot(samptimes,c_ts9.*0,'Color',[0 0 0],'LineWidth',3)
% xlim(plotwin), xlabel('Time relative to sample onset (s)'), ylabel('\rho')
% set(gca,'TickDir','out','box','off'), title('PPI vs ',names{3})





% % Bin by CAPE score
% figure
% for sc = 1:length(names)
%     
%     fprintf('\nProcessing scale %d of %d for final plot...',sc,length(names))
%     
% %     if strcmp(names{sc},'P')
% %         % if P-scale, use defined cutoff
% %         gL = find(CAPE(include,sc)<8);
% %         gH = find(CAPE(include,sc)>=8);
% %     elseif strcmp(names{sc},'P_p_e_r')
% %         % if P-perceptual-subscale, compare zero to non-zero
% %         gL = find(CAPE(include,sc)==0);
% %         gH = find(CAPE(include,sc)>0);
% %     else
% %         % get groups via median split
% %         gL = find(CAPE(include,sc)<nanmedian(CAPE(include,sc)));
% %         gH = find(CAPE(include,sc)>nanmedian(CAPE(include,sc)));
% %     end
%     
%     % get groups via percentile split
%     if strcmp(names{sc},'P_p_e_r')
%         % if P-perceptual-subscale, compare zero to non-zero
%         gL = find(CAPE(include,sc)==0);
%         gH = find(CAPE(include,sc)>0);
%     else
%         cuts = prctile(CAPE(include,sc),[20 80]);
%         gL = find(CAPE(include,sc)<=cuts(1));
%         gH = find(CAPE(include,sc)>=cuts(2));
%     end
%     
%     subplot(3,5,sc), hold on
%     dat = GAdil_full(include,:);
%     p=[];
%     for t = 1:size(GAdil_full,2)
%         p(t)=permutationTest(dat(gL,t), dat(gH,t), nperm);
%     end
%     h=nan(size(p)); h(p<0.05)=1;
%     shadedErrorBar(fulltimes,mean(dat(gL,:),1),std(dat(gL,:),[],1)./sqrt(length(gL)),{'Color',[0.3 0.1 1],'LineWidth',1.5,'LineStyle','--'},0)
%     shadedErrorBar(fulltimes,mean(dat(gH,:),1),std(dat(gH,:),[],1)./sqrt(length(gH)),{'Color',[0.3 0.1 1],'LineWidth',1.5},0)
%     plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
%     l = line([smarks],repmat(get(gca, 'ylim')',1,size(smarks,2))); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0.6 0.6 0.6])
%     plot(fulltimes,h.*0,'Color',[0 0 0],'LineWidth',3)
%     xlim(fullwin), xlabel('Time relative to trial onset (s)'), ylabel('dilation (zs^-^1)')
%     set(gca,'TickDir','out','box','off'), title(['Dilation response vs ',names{sc}])
%     
%     subplot(3,5,sc+5), hold on
%     dat1 = squeeze(mean(GAdil_pCP_Bext(include,:,:,1),2));
%     dat2 = squeeze(mean(GAdil_pCP_Bext(include,:,:,2),2));
%     dat3 = squeeze(mean(GAdil_pCP_Bext(include,:,:,3),2));
%     p1=[]; p2=[]; p3=[];
%     for t = 1:size(GAdil_pCP_Bext,3)
%         p1(t)=permutationTest(dat1(gL,t), dat1(gH,t), nperm);
%         p2(t)=permutationTest(dat2(gL,t), dat2(gH,t), nperm);
%         p3(t)=permutationTest(dat3(gL,t), dat3(gH,t), nperm);
%     end
%     h1=nan(size(p1)); h1(p1<0.05)=1;
%     h2=nan(size(p2)); h2(p2<0.05)=1;
%     h3=nan(size(p3)); h3(p3<0.05)=1;
%     shadedErrorBar(samptimes,mean(dat1(gL,:),1),std(dat1(gL,:),[],1)./sqrt(length(gL)),{'Color',[1 0 0],'LineWidth',1.5,'LineStyle','--'},0)
%     shadedErrorBar(samptimes,mean(dat1(gH,:),1),std(dat1(gH,:),[],1)./sqrt(length(gH)),{'Color',[1 0 0],'LineWidth',1.5},0)
%     shadedErrorBar(samptimes,mean(dat2(gL,:),1),std(dat2(gL,:),[],1)./sqrt(length(gL)),{'Color',[0 0 1],'LineWidth',1.5,'LineStyle','--'},0)
%     shadedErrorBar(samptimes,mean(dat2(gH,:),1),std(dat2(gH,:),[],1)./sqrt(length(gH)),{'Color',[0 0 1],'LineWidth',1.5},0)
%     shadedErrorBar(samptimes,mean(dat3(gL,:),1),std(dat3(gL,:),[],1)./sqrt(length(gL)),{'Color',[0 1 0],'LineWidth',1.5,'LineStyle','--'},0)
%     shadedErrorBar(samptimes,mean(dat3(gH,:),1),std(dat3(gH,:),[],1)./sqrt(length(gH)),{'Color',[0 1 0],'LineWidth',1.5},0)
%     plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
%     plot(samptimes,h1.*0.003,'Color',[1 0 0],'LineWidth',3)
%     plot(samptimes,h2.*-0.003,'Color',[0 0 1],'LineWidth',3)
%     plot(samptimes,h3.*-0,'Color',[0 1 0],'LineWidth',3)
%     xlim(plotwin), xlabel('Time relative to sample onset (s)'), ylabel('beta (a.u.)')
%     set(gca,'TickDir','out','box','off'), title(['Encoding vs ',names{sc}])
%     
%     subplot(3,5,sc+10), hold on
%     dat = squeeze(mean(GAdil_PPI_Bext_lasso(include,:,:),2));
%     p=[];
%     for t = 1:size(GAdil_PPI_Bext_lasso,3)
%         p(t)=permutationTest(dat(gL,t), dat(gH,t), nperm);
%     end
%     h=nan(size(p)); h(p<0.05)=1;
%     shadedErrorBar(samptimes,mean(dat(gL,:),1),std(dat(gL,:),[],1)./sqrt(length(gL)),{'Color',[0 0 0],'LineWidth',1.5,'LineStyle','--'},0)
%     shadedErrorBar(samptimes,mean(dat(gH,:),1),std(dat(gH,:),[],1)./sqrt(length(gH)),{'Color',[0 0 0],'LineWidth',1.5},0)
%     plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
%     plot(samptimes,h.*0,'Color',[0 0 0],'LineWidth',3)
%     xlim(plotwin), xlabel('Time relative to sample onset (s)'), ylabel('beta (a.u.)')
%     set(gca,'TickDir','out','box','off'), title(['PPI vs ',names{sc}])
%     
% end



% temp FIGURE
corrtype = 'spearman';
colL = [0.5 0.5 0.5];
colH = [0 0 0];
axlw=0.5; fs = 7;
jitbnd = [-0.3 0.3];
scatsize = 9;

cuts = prctile(CAPE(include,strcmp(names,'P')),[20 80]);
gL = find(CAPE(include,strcmp(names,'P'))<=cuts(1));
gH = find(CAPE(include,strcmp(names,'P'))>=cuts(2));

figure

% trial-evoked response
% subplot(3,5,1), hold on
% dat = GAdil_full(include,:);
% p=[];
% for t = 1:size(GAdil_full,2)
%     p(t)=permutationTest(dat(gL,t), dat(gH,t), nperm);
% end
% h=nan(size(p)); h(p<0.05)=1;
% shadedErrorBar(fulltimes,mean(dat(gL,:),1),std(dat(gL,:),[],1)./sqrt(length(gL)),{'Color',[0 0 0],'LineWidth',1.5,'LineStyle','--'},0)
% shadedErrorBar(fulltimes,mean(dat(gH,:),1),std(dat(gH,:),[],1)./sqrt(length(gH)),{'Color',[0 0 0],'LineWidth',1.5},0)
% plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
% l = line([smarks],repmat(get(gca, 'ylim')',1,size(smarks,2))); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0.6 0.6 0.6])
% plot(fulltimes,h.*0,'Color',[0 0 0],'LineWidth',3)
% xlim(fullwin), xlabel('Time relative to trial onset (s)'), ylabel('dilation (zs^-^1)')
% set(gca,'FontName','Helvetica','FontSize',fs,'box','off','TickDir','out','LineWidth',axlw), title('Trial-evoked response')


subplot(3,5,6), hold on
% dat = dil_proj(include);
dat = dil_mean(include);

rnds = dist_scatters(dat(gL),0.1);  % get x jitters
S=scatter(ones(size(dat(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),dat(gL),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
plot(jitbnd+1,[mean(dat(gL)) mean(dat(gL))],'LineWidth',2,'Color',colL)

rnds = dist_scatters(dat(gH),0.1);  % get x jitters
S=scatter(ones(size(dat(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),dat(gH),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
plot(jitbnd+2,[mean(dat(gH)) mean(dat(gH))],'LineWidth',2,'Color',colH)

p = permutationTest(dat(gL), dat(gH), nperm);
text(1.1,2,['p=',num2str(round(p,4))],'FontSize',7)

set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
ylabel('Dilation lin. proj.'), xlim([0.5 2.5])


subplot(3,5,11), hold on
[rho,p]=corr(CAPE(include,strcmp(names,'P')),dat,'type',corrtype,'rows','pairwise');
scatter(CAPE(include,strcmp(names,'P')),dat,scatsize), ylabel('Dilation lin. proj.'), xlabel(names{1})
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')


% var encoding
cols = [1 0 0; 0 0 1; 0 1 0];
spnames = {'CPP','Unc','|LLR|'};
for sp = 2:4
%     subplot(3,5,sp), hold on
%     dat = squeeze(mean(GAdil_pCP_Bext(include,:,:,sp-1),2));
%     p1=[];
%     for t = 1:size(GAdil_pCP_Bext,3)
%         p1(t)=permutationTest(dat(gL,t), dat(gH,t), nperm);
%     end
%     h1=nan(size(p1)); h1(p1<0.05)=1;
%     shadedErrorBar(samptimes,mean(dat(gL,:),1),std(dat(gL,:),[],1)./sqrt(length(gL)),{'Color',cols(sp-1,:),'LineWidth',1.5,'LineStyle','--'},0)
%     shadedErrorBar(samptimes,mean(dat(gH,:),1),std(dat(gH,:),[],1)./sqrt(length(gH)),{'Color',cols(sp-1,:),'LineWidth',1.5},0)
%     plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
%     plot(samptimes,h1.*-0.003,'Color',[1 0 0],'LineWidth',3)
%     xlim(plotwin), xlabel('Time relative to sample onset (s)'), ylabel('beta (a.u.)')
%     set(gca,'FontName','Helvetica','FontSize',fs,'box','off','TickDir','out','LineWidth',axlw), title([spnames{sp-1},' encoding'])
    
    
    subplot(3,5,sp+5), hold on
    dat = dil_pCP_Bext_av(include,sp-1);
    
    rnds = dist_scatters(dat(gL),0.1);  % get x jitters
    S=scatter(ones(size(dat(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),dat(gL),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
    plot(jitbnd+1,[mean(dat(gL)) mean(dat(gL))],'LineWidth',2,'Color',cols(sp-1,:),'LineStyle','--')
    
    rnds = dist_scatters(dat(gH),0.1);  % get x jitters
    S=scatter(ones(size(dat(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),dat(gH),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
    plot(jitbnd+2,[mean(dat(gH)) mean(dat(gH))],'LineWidth',2,'Color',cols(sp-1,:))
    
    p = permutationTest(dat(gL), dat(gH), nperm);
    text(1.1,0,['p=',num2str(round(p,4))],'FontSize',7)
    
    set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
    ylabel('Encoding lin. proj. (a.u.)'), xlim([0.5 2.5])
    
    
    subplot(3,5,sp+10), hold on
    [rho,p]=corr(CAPE(include,strcmp(names,'P')),dat,'type',corrtype,'rows','pairwise');
    scatter(CAPE(include,strcmp(names,'P')),dat,scatsize,cols(sp-1,:)), ylabel('Encoding lin. proj. (a.u.)'), xlabel(names{1})
    title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])
    set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')
end 


% choice effect
% subplot(3,5,5), hold on
% dat = squeeze(mean(GAdil_PPI_Bext_lasso(include,:,:),2));
% p=[];
% for t = 1:size(GAdil_PPI_Bext_lasso,3)
%     p(t)=permutationTest(dat(gL,t), dat(gH,t), nperm);
% end
% h=nan(size(p)); h(p<0.05)=1;
% shadedErrorBar(samptimes,mean(dat(gL,:),1),std(dat(gL,:),[],1)./sqrt(length(gL)),{'Color',[0 0 0],'LineWidth',1.5,'LineStyle','--'},0)
% shadedErrorBar(samptimes,mean(dat(gH,:),1),std(dat(gH,:),[],1)./sqrt(length(gH)),{'Color',[0 0 0],'LineWidth',1.5},0)
% plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
% plot(samptimes,h.*0,'Color',[0 0 0],'LineWidth',3)
% xlim(plotwin), xlabel('Time relative to sample onset (s)'), ylabel('beta (a.u.)')
% set(gca,'FontName','Helvetica','FontSize',fs,'box','off','TickDir','out','LineWidth',axlw), title('Pupil choice effect')


subplot(3,5,10), hold on
dat = PPI_proj(include);
% dat = dil_pCP_Bext_av(include,1) + dil_pCP_Bext_av(include,2);

rnds = dist_scatters(dat(gL),0.1);  % get x jitters
S=scatter(ones(size(dat(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),dat(gL),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
plot(jitbnd+1,[mean(dat(gL)) mean(dat(gL))],'LineWidth',2,'Color',[0 0 0],'LineStyle','--')

rnds = dist_scatters(dat(gH),0.1);  % get x jitters
S=scatter(ones(size(dat(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),dat(gH),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
plot(jitbnd+2,[mean(dat(gH)) mean(dat(gH))],'LineWidth',2,'Color',[0 0 0])

p = permutationTest(dat(gL), dat(gH), nperm);
text(1.1,2,['p=',num2str(round(p,4))],'FontSize',7)

set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
ylabel('Choice eff. lin. proj. (a.u.)'), xlim([0.5 2.5])


subplot(3,5,15), hold on
[rho,p]=corr(CAPE(include,strcmp(names,'P')),dat,'type',corrtype,'rows','pairwise');
scatter(CAPE(include,strcmp(names,'P')),dat,scatsize,[0 0 0]), ylabel('Choice eff. lin. proj. (a.u.)'), xlabel(names{1})
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')



    


% Figure
fig_w = 22; % figure width
fig_h = 11; % figure height

hf =  findobj('type','figure');
cfig = length(hf)+1;
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),


% Plotting average trial-evoked response
smarks = repmat(0.4:0.4:0.4*10,2,1);
sampcols = [linspace(0.8,0,size(GAdil_samp,2))' linspace(0.8,0,size(GAdil_samp,2))' linspace(0.8,0,size(GAdil_samp,2))'];

subplot(3,5,1:2), hold on,
lx = line(fullwin,[0 0]); set(lx,'LineWidth',0.5,'LineStyle','--','Color',[0 0 0])
% for subj = include
%     plot(fulltimes,GAdil_full(subj,:),'Color',[0.7 0.7 0.7],'LineWidth',0.75)
% end
se = std(GAdil_full(include,:),[],1)./sqrt(length(include));
if first_deriv
    shadedErrorBar(fulltimes,mean(GAdil_full(include,:),1),se,{'Color',[0 0 0],'LineWidth',1.5,'LineStyle','--'},0);
else shadedErrorBar(fulltimes,mean(GAdil_full(include,:),1),se,{'Color',[0 0 0],'LineWidth',1.5},0);
end
plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5)
l = line([smarks],repmat(get(gca, 'ylim')',1,size(smarks,2))); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0.6 0.6 0.6])
% plot(fulltimes,sig_ts7.*-0.1,'k','LineWidth',3)
xlim([-0.5 5.5]), xlabel('Time relative to trial onset (s)'), ylabel('Pupil dilation (z / zs-1)')
set(gca,'TickDir','out','box','off','FontName','Helvetica','FontSize',fs,'LineWidth',axlw)

if first_deriv
    subplot(3,5,7), hold on
else subplot(3,5,6), hold on
end
% dat = dil_proj(include);
dat = dil_mean(include);

rnds = dist_scatters(dat(gL),0.1);  % get x jitters
S=scatter(ones(size(dat(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),dat(gL),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
plot(jitbnd+1,[mean(dat(gL)) mean(dat(gL))],'LineWidth',2,'Color',colL)

rnds = dist_scatters(dat(gH),0.1);  % get x jitters
S=scatter(ones(size(dat(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),dat(gH),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
plot(jitbnd+2,[mean(dat(gH)) mean(dat(gH))],'LineWidth',2,'Color',colH)

p = permutationTest(dat(gL), dat(gH), nperm);
text(1.1,1.2,['p=',num2str(round(p,4))],'FontSize',7)

set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
ylabel('Dilation lin. proj.'), xlim([0.5 2.5])
title([num2str(dil_mean_win(1)),'-',num2str(dil_mean_win(2)),' s'])


if first_deriv
    subplot(3,5,12), hold on
else subplot(3,5,11), hold on
end
[rho,p]=corr(CAPE(include,strcmp(names,'P')),dat,'type',corrtype,'rows','pairwise');
scatter(CAPE(include,strcmp(names,'P')),dat,scatsize,'k'), ylabel('Dilation lin. proj.'), xlabel(names{1})
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')


% Plotting sample-wise variable encoding
subplot(3,5,4:5), hold on,
seS = std(squeeze(mean(GAdil_pCP_Bext(include,:,plotts,1),2)),[],1)./sqrt(length(include));
seL = std(squeeze(mean(GAdil_pCP_Bext(include,:,plotts,2),2)),[],1)./sqrt(length(include));
seG = std(squeeze(mean(GAdil_pCP_Bext(include,:,plotts,3),2)),[],1)./sqrt(length(include));
shadedErrorBar(plottimes,squeeze(mean(mean(GAdil_pCP_Bext(include,:,plotts,1),2),1)),seS,{'Color',[1 0 0],'LineWidth',1.5},0);
shadedErrorBar(plottimes,squeeze(mean(mean(GAdil_pCP_Bext(include,:,plotts,2),2),1)),seL,{'Color',[0 0 1],'LineWidth',1.5},0);
shadedErrorBar(plottimes,squeeze(mean(mean(GAdil_pCP_Bext(include,:,plotts,3),2),1)),seG,{'Color',[0 1 0],'LineWidth',1.5},0);
% plot(plottimes,sig_ts2.*-0.014,'r','LineWidth',3)
% plot(plottimes,sig_ts3.*-0.016,'b','LineWidth',3)
% plot(plottimes,sig_ts4.*-0.018,'g','LineWidth',3)
l = line(plotwin,[0 0]); set(l,'LineWidth',0.5,'LineStyle','--')
xlim([0 1.0]), ylabel('beta (a.u.)'), xlabel('Time relative to sample onset (s)')
% legend([p1 p2],{'CPP','-|psi|'})
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')

cols = [1 0 0; 0 0 1; 0 1 0];
spnames = {'CPP','Unc'};
for sp = 2:3
%     subplot(3,5,sp), hold on
%     dat = squeeze(mean(GAdil_pCP_Bext(include,:,:,sp-1),2));
%     p1=[];
%     for t = 1:size(GAdil_pCP_Bext,3)
%         p1(t)=permutationTest(dat(gL,t), dat(gH,t), nperm);
%     end
%     h1=nan(size(p1)); h1(p1<0.05)=1;
%     shadedErrorBar(samptimes,mean(dat(gL,:),1),std(dat(gL,:),[],1)./sqrt(length(gL)),{'Color',cols(sp-1,:),'LineWidth',1.5,'LineStyle','--'},0)
%     shadedErrorBar(samptimes,mean(dat(gH,:),1),std(dat(gH,:),[],1)./sqrt(length(gH)),{'Color',cols(sp-1,:),'LineWidth',1.5},0)
%     plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
%     plot(samptimes,h1.*-0.003,'Color',[1 0 0],'LineWidth',3)
%     xlim(plotwin), xlabel('Time relative to sample onset (s)'), ylabel('beta (a.u.)')
%     set(gca,'FontName','Helvetica','FontSize',fs,'box','off','TickDir','out','LineWidth',axlw), title([spnames{sp-1},' encoding'])
    
    
    subplot(3,5,sp+6), hold on
    dat = dil_pCP_Bext_av(include,sp-1);
    
    rnds = dist_scatters(dat(gL),0.1);  % get x jitters
    S=scatter(ones(size(dat(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),dat(gL),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
    plot(jitbnd+1,[mean(dat(gL)) mean(dat(gL))],'LineWidth',2,'Color',cols(sp-1,:),'LineStyle','--')
    
    rnds = dist_scatters(dat(gH),0.1);  % get x jitters
    S=scatter(ones(size(dat(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),dat(gH),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
    plot(jitbnd+2,[mean(dat(gH)) mean(dat(gH))],'LineWidth',2,'Color',cols(sp-1,:))
    
    p = permutationTest(dat(gL), dat(gH), nperm);
    text(1.1,0,['p=',num2str(round(p,4))],'FontSize',7)
    
    set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
    ylabel('beta (a.u.)'), xlim([0.5 2.5])
    
    if sp==2
        title([num2str(cpp_win(1)),'-',num2str(cpp_win(2)),' s'])
    elseif sp==3
        title([num2str(unc_win(1)),'-',num2str(unc_win(2)),' s'])
    elseif sp==4
        title([num2str(llr_win(1)),'-',num2str(llr_win(2)),' s'])
    end
    
    subplot(3,5,sp+11), hold on
    [rho,p]=corr(CAPE(include,strcmp(names,'P')),dat,'type',corrtype,'rows','pairwise');
    scatter(CAPE(include,strcmp(names,'P')),dat,scatsize,cols(sp-1,:)), ylabel('beta (a.u.)'), xlabel(names{1})
    title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])
    set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')
end 

sp = sp+1;
subplot(3,5,sp+6), hold on
dat = dil_pCP_Bext_av(include,1) + dil_pCP_Bext_av(include,2);

rnds = dist_scatters(dat(gL),0.1);  % get x jitters
S=scatter(ones(size(dat(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),dat(gL),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
plot(jitbnd+1,[mean(dat(gL)) mean(dat(gL))],'LineWidth',2,'Color',cols(sp-1,:),'LineStyle','--')

rnds = dist_scatters(dat(gH),0.1);  % get x jitters
S=scatter(ones(size(dat(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),dat(gH),scatsize,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
plot(jitbnd+2,[mean(dat(gH)) mean(dat(gH))],'LineWidth',2,'Color',cols(sp-1,:))

p = permutationTest(dat(gL), dat(gH), nperm);
text(0.2,0,['p=',num2str(round(p,4))],'FontSize',7)

set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{1}],['High ',names{1}]})
ylabel('beta (a.u.)'), xlim([0.5 2.5])

if sp==2
    title([num2str(cpp_win(1)),'-',num2str(cpp_win(2)),' s'])
elseif sp==3
    title([num2str(unc_win(1)),'-',num2str(unc_win(2)),' s'])
elseif sp==4
    title('CPP + Unc')
end

subplot(3,5,sp+11), hold on
[rho,p]=corr(CAPE(include,strcmp(names,'P')),dat,'type',corrtype,'rows','pairwise');
scatter(CAPE(include,strcmp(names,'P')),dat,scatsize,cols(sp-1,:)), ylabel('beta (a.u.)'), xlabel(names{1})
title(['rho=',num2str(round(rho,3)),', p=',num2str(round(p,4))])
set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off')



% Multiple regression of P on pupil encs
m = regstats(CAPE(include,strcmp(names,'P')),[dil_pCP_Bext_av(include,1:2)],'linear');



% Kernel & pupil measures quintile-binned by N and D scores
transbars = 0;

fig_w = 22; % figure width
fig_h = 11; % figure height

h =  findobj('type','figure');
cfig = length(h)+1;
f=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

for sc = 2:3
        
    % get groups via percentile split
    cuts = prctile(CAPE(include,sc),[20 80]);
    gL = find(CAPE(include,sc)<=cuts(1));
    gH = find(CAPE(include,sc)>=cuts(2));
    
    % plot LLR kernels
    subplot(2,6,1+6*(sc-2)), hold on
    sig_ts=[]; for t=1:maxsamps, sig_ts(t)=permutationTest(rGA_ssB_regU(gL,t), rGA_ssB_regU(gH,t), nperm); end
    shadedErrorBar(1:maxsamps,mean(rGA_ssB_regU(gL,:),1),std_err(rGA_ssB_regU(gL,:),1),{'Color',colL},transbars)
    shadedErrorBar(1:maxsamps,mean(rGA_ssB_regU(gH,:),1),std_err(rGA_ssB_regU(gH,:),1),{'Color',colH},transbars)
    s2=scatter(find(sig_ts<=0.05),ones(1,length(find(sig_ts<0.05))).*0,20); set(s2,'MarkerEdgeColor',[0.8 0 0])
    plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.75 0.75 0.75],'LineWidth',0.75)
    xlim([0.5 maxsamps+0.5])
    xlabel('Sample position'), ylabel('Beta weight (a.u.)'), set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[2 6 10])
    if sc==1, title('Evidence weighting'), end
    
    % plot LLR kernel 2nd minus 1st halves
    subplot(2,6,2+6*(sc-2)), hold on
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
    
    set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{sc}],['High ',names{sc}]})
    ylabel('delta beta weight (a.u.)'), xlim([0.5 2.5])
    if sc==1, title('Kernel half difference'), end
    
    % plot surprise kernels
    subplot(2,6,3+6*(sc-2)), hold on
    sig_ts=[]; for t=1:maxsamps-1, sig_ts(t)=permutationTest(rGA_ssB_surpU(gL,t), rGA_ssB_surpU(gH,t), nperm); end
    shadedErrorBar(2:maxsamps,mean(rGA_ssB_surpU(gL,:),1),std_err(rGA_ssB_surpU(gL,:),1),{'Color',colL},transbars)
    shadedErrorBar(2:maxsamps,mean(rGA_ssB_surpU(gH,:),1),std_err(rGA_ssB_surpU(gH,:),1),{'Color',colH},transbars)
    s2=scatter(find(sig_ts<=0.05)+1,ones(1,length(find(sig_ts<0.05))).*0,20); set(s2,'MarkerEdgeColor',[0.8 0 0])
    plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.75 0.75 0.75],'LineWidth',0.75)
    xlim([0.5 maxsamps+0.5])
    xlabel('Sample position'), ylabel('Beta weight (a.u.)'), set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[2 6 10])
    if sc==1, title('CPP modulation'), end
        
    % plot average magnitude of CPP kernels
    subplot(2,6,4+6*(sc-2)), hold on
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
    
    set(gca,'FontName','Arial','LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{sc}],['High ',names{sc}]})
    ylabel('Mean beta weight (a.u.)'), xlim([0.5 2.5])
    if sc==1, title('CPP +UNC kernel modulations'), end
    
    % Dilation mean response
    subplot(2,6,5+6*(sc-2)), hold on
    dat = dil_mean(include);
    
    rnds = dist_scatters(dat(gL),0.1);  % get x jitters
    S=scatter(ones(size(dat(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),dat(gL),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
    plot(jitbnd+1,[mean(dat(gL)) mean(dat(gL))],'LineWidth',2,'Color',colL)
    
    rnds = dist_scatters(dat(gH),0.1);  % get x jitters
    S=scatter(ones(size(dat(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),dat(gH),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
    plot(jitbnd+2,[mean(dat(gH)) mean(dat(gH))],'LineWidth',2,'Color',colH)
    
    p = permutationTest(dat(gL), dat(gH), nperm);
    text(1.1,1.2,['p=',num2str(round(p,4))],'FontSize',7)
    
    set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{sc}],['High ',names{sc}]})
    ylabel('Mean response'), xlim([0.5 2.5]), ylim([-0.2 1.22])
    
    % Pupil CPP+Unc encoding
    subplot(2,6,6+6*(sc-2)), hold on
    dat = dil_pCP_Bext_av(include,1) + dil_pCP_Bext_av(include,2);
    
    rnds = dist_scatters(dat(gL),0.1);  % get x jitters
    S=scatter(ones(size(dat(gL),1),1).*1 + diff(jitbnd).*rnds + jitbnd(1),dat(gL),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colL)
    plot(jitbnd+1,[mean(dat(gL)) mean(dat(gL))],'LineWidth',2,'Color',cols(sp-1,:),'LineStyle','--')
    
    rnds = dist_scatters(dat(gH),0.1);  % get x jitters
    S=scatter(ones(size(dat(gH),1),1).*2 + diff(jitbnd).*rnds + jitbnd(1),dat(gH),scatsize,'k','o');
    set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colH)
    plot(jitbnd+2,[mean(dat(gH)) mean(dat(gH))],'LineWidth',2,'Color',cols(sp-1,:))
    
    p = permutationTest(dat(gL), dat(gH), nperm);
    text(1.1,0.08,['p=',num2str(round(p,4))],'FontSize',7)
    
    set(gca,'FontName','Helvetica','FontSize',fs,'LineWidth',axlw,'TickDir','out','box','off','XTick',[1:2],'XTickLabel',{['Low ',names{sc}],['High ',names{sc}]})
    ylabel('Beta weight (a.u.)'), xlim([0.5 2.5])
end
    

