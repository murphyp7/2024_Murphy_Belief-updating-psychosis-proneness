clear, clc
close all

% param settings
ntrials = [1376];  % full range: [1376, 10000]
Hs = [0.03 0.06 0.1];  % full range: [0.03 0.06 0.1]
noises = [0.001 2];  % full range: [0.001 2]
IUs = [1 1.4];    % full range: [1 1.4]

% path details
loadpath = '/mnt/homes/home024/pmurphy/Surprise_scz/param_recovery/sims/';

%%%%%%%%%%%%%
%%% NO IU %%%
%%%%%%%%%%%%%
p_noIU=[];
for t = 1:length(ntrials)
    for h = 1:length(Hs)
        for n = 1:length(noises)
            load([loadpath,'Param_rec_H',num2str(Hs(h)),'_noise',num2str(noises(n)),'_t',num2str(ntrials(t)),'.mat'],'pm_fit')
            p_noIU(t,h,n,:,:) = pm_fit;
        end
    end
end

Hcols = [1 0.5 0.5; 1 0 0; 0.5 0 0];
for t = 1:length(ntrials)
    figure,  % plot param[H,noise]*noise
    for n  = 1:length(noises)
        subplot(length(noises),2,n), hold on
        for h = 1:length(Hs)
            histogram(p_noIU(t,h,n,:,1),0:0.002:0.15,'DisplayStyle','stairs','EdgeColor',Hcols(h,:))
            plot([Hs(h) Hs(h)],[-2 20],'Color',Hcols(h,:),'LineStyle','--','LineWidth',2.5)
            plot([median(p_noIU(t,h,n,:,1)) median(p_noIU(t,h,n,:,1))],[-2 20],'Color',Hcols(h,:),'LineStyle','--','LineWidth',0.5)
        end
        xlabel('H'), ylabel('count'), title(['noise=',num2str(noises(n)),', ntrials=',num2str(ntrials(t))])
        subplot(length(noises),2,n+length(noises)), hold on
        for h = 1:length(Hs)
            histogram(p_noIU(t,h,n,:,2),40,'DisplayStyle','stairs','EdgeColor',Hcols(h,:))
            plot([noises(n) noises(n)],[-2 16],'Color',Hcols(h,:),'LineStyle','--','LineWidth',2.5)
            plot([median(p_noIU(t,h,n,:,2)) median(p_noIU(t,h,n,:,2))],[-2 16],'Color',Hcols(h,:),'LineStyle','--','LineWidth',0.5)
        end
        xlabel('noise'), ylabel('count')
    end
end


%%%%%%%%%%%%%
%%% WITH IU %%%
%%%%%%%%%%%%%
p_IU=[];
for t = 1:length(ntrials)
    for h = 1:length(Hs)
        for n = 1:length(noises)
            for i = 1:length(IUs)
                load([loadpath,'Param_rec_H',num2str(Hs(h)),'_noise',num2str(noises(n)),'_IU',num2str(IUs(i)),'_t',num2str(ntrials(t)),'.mat'],'pm_fit')
                p_IU(t,h,n,i,:,:) = pm_fit;
            end
        end
    end
end

for t = 1:length(ntrials)
    figure,  % plot param[H,noise,IU]*(noise+IU)
    for n  = 1:length(noises)
        for i = 1:length(IUs)
            subplot(3,length(noises)+length(IUs),n+length(noises)*(i-1)), hold on
            for h = 1:length(Hs)
                histogram(p_IU(t,h,n,i,:,1),0:0.002:0.15,'DisplayStyle','stairs','EdgeColor',Hcols(h,:))
                plot([Hs(h) Hs(h)],[-2 20],'Color',Hcols(h,:),'LineStyle','--','LineWidth',2.5)
                plot([median(p_IU(t,h,n,i,:,1)) median(p_IU(t,h,n,i,:,1))],[-2 20],'Color',Hcols(h,:),'LineStyle','--','LineWidth',0.5)
            end
            xlabel('H'), ylabel('count'), title(['noise=',num2str(noises(n)),', IU=',num2str(IUs(i)),', ntrials=',num2str(ntrials(t))])
            
            subplot(3,length(noises)+length(IUs),n+length(noises)*(i-1)+length(noises)+length(IUs)), hold on
            for h = 1:length(Hs)
                histogram(p_IU(t,h,n,i,:,2),40,'DisplayStyle','stairs','EdgeColor',Hcols(h,:))
                plot([noises(n) noises(n)],[-2 16],'Color',Hcols(h,:),'LineStyle','--','LineWidth',2.5)
                plot([median(p_IU(t,h,n,i,:,2)) median(p_IU(t,h,n,i,:,2))],[-2 16],'Color',Hcols(h,:),'LineStyle','--','LineWidth',0.5)
            end
            xlabel('noise'), ylabel('count')
            
            subplot(3,length(noises)+length(IUs),n+length(noises)*(i-1)+(length(noises)+length(IUs))*2), hold on
            for h = 1:length(Hs)
                histogram(p_IU(t,h,n,i,:,3),40,'DisplayStyle','stairs','EdgeColor',Hcols(h,:))
                plot([IUs(i) IUs(i)],[-2 16],'Color',Hcols(h,:),'LineStyle','--','LineWidth',2.5)
                plot([median(p_IU(t,h,n,i,:,3)) median(p_IU(t,h,n,i,:,3))],[-2 16],'Color',Hcols(h,:),'LineStyle','--','LineWidth',0.5)
            end
            xlabel('IU'), ylabel('count')
        end
    end
end

% correlate parameter estimates
for t = 1:length(ntrials)
    figure,
    for n  = 1:length(noises)
        for i = 1:length(IUs)
            for h = 1:length(Hs)
                subplot(length(Hs),length(noises)+length(IUs),n+length(noises)*(i-1)+(length(noises)+length(IUs))*(h-1)), hold on
                s=scatter(squeeze(p_IU(t,h,n,i,:,1)),squeeze(p_IU(t,h,n,i,:,3)),10); set(s,'MarkerEdgeColor',Hcols(h,:))
                l=lsline; set(l,'Color',[0.3 0.3 0.3],'LineWidth',1.5)
                s2=scatter(Hs(h),IUs(i),40,'k','+'); set(s2,'MarkerFaceColor',[0 0 0],'LineWidth',2.5)
                xlabel('H'), ylabel('IU'), title(['noise=',num2str(noises(n)),', ntrials=',num2str(ntrials(t))])
            end
        end
    end
end
