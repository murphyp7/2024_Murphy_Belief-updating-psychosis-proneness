% Analyze single-subject pupil data

function analyze_pupil(subj)

% Path stuff
loadpath = '/mnt/homes/home024/pmurphy/Surprise_scz/pupil/4.interpolated/';
savepath = '/mnt/homes/home024/pmurphy/Surprise_scz/pupil/5.output/';
behavpath = '/mnt/homes/home028/gmonov/SCZ/Data/decision_making/';
modelpath = '/mnt/homes/home024/pmurphy/Surprise_scz/modelling/';

addpath '/mnt/homes/home024/pmurphy/Surprise_scz/pupil/'
addpath '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/Gen_fun';
addpath(genpath('/mnt/homes/home024/pmurphy/Surprise_scz/'))

% get list of all subject IDs
files = dir([loadpath,'*.mat']);
for f = 1:length(files), allsubj{f}=files(f).name(1:4); end
allsubj = unique(allsubj);

% pull files for only current subject
subj = allsubj{subj};
files = dir([loadpath,subj,'*.mat']);

% Analysis stuff
maxsamps = 10;

freqs = [0.06 6];  % filter cutoffs [lo hi]
newFs = 50;  % new sampling rate if desired (set to empty if no resampling required)

basewin = [-0.05 0.05];  % window for baselining, relative to eliciting event (default = [-0.1 0])
fullwin = [-2 6.2];  % window for plotting full dilation response aligned to trial onset
sampwin = [0 1.4];  % window for plotting dilation response to indiviudal samples (default = [-0.15 1.5])
fbwin = [-0.5 3.5];    % window for plotting dilation response to feeback

artwin = [0.4 5];   % window (aligned to trial onset) for applying artifact rejection
artlen = 1.0;    % length (in secs) of window for definition of 'prolonged' artifact (if any sample from such an artifactual period falls within artwin, trial is rejected)
artp = 0.6;   % proportion of samples for definition of 'noisy' artifact (if length(interp_smp)/length(all_smp)>artp, trial is rejected)

first_deriv = 0;      % take first derivative of pupil signal rather than original time course (will not impose any baseline)
singlesamp_base = 0;  % baseline at trialwise or samplewise levels for regressions
lin_detrend = 0;      % linearly detrend single-trial pupil epochs before regressions
smp_regout_lag = 1;   % in timeseries regressions, determines variables from how many samples prior to sample n will be included as covariates
quad_reg = 0;         % include quandratic term in pupil regressions
outliercut = inf;     % |z|-threshold for throwing out observations from regression models

modeltype = 'H_noise';

% Loop through files
dil_full=[]; dil_samp=[]; bad_ts=[]; dil_fbC=[]; dil_fbE=[]; smpbase_full=[];
LLR_full=[]; LPR_full=[]; psi_full=[]; pCP_full=[]; choices_full=[];
X_samp=[]; Y_samp=[]; pupilPSD=[]; mfit=struct;
errfiles={};
for f = 1:length(files)
    
    fprintf('File %s \n',files(f).name)  % print progress
    try
        % Load behaviour and sample sequences
        s = files(f).name(6);
        b = files(f).name(8:end-4);
        load([behavpath,subj,filesep,'S',s,filesep,'Behaviour',filesep,subj,'_',s,'_',b,'.mat'])
        load([behavpath,subj,filesep,'S',s,filesep,'Sample_seqs',filesep,subj,'_',s,'_',b,'.mat'])
        
        % Converting sample and choice values to appropriate signs for choice regressions
        stimIn = round(stimIn.*-1);
        choices = Behav(:,2)-1;
        pIn = cat(3,normpdf(stimIn,gen.mu(1),gen.sigma(1)),normpdf(stimIn,gen.mu(2),gen.sigma(2)));  % make matrix of sample probabilities for both left (:,:,1) and right (:,:,2) distributions
        LLRin = log(normpdf(stimIn,gen.mu(1),gen.sigma(2))./normpdf(stimIn,gen.mu(2),gen.sigma(1)));
        
        % Calculate model-based variables
        if strcmp(modeltype,'ideal')
            [LPRout,pCP,scaled_prior] = accGlaze_fast(LLRin,gen.H,0,'pCP',pIn);
            
        elseif strcmp(modeltype,'H_noise')
            load([modelpath,modeltype,filesep,'fits',filesep,subj,'_fit.mat'])
            mfit.H = pm_fit(1);
            mfit.noise = pm_fit(2);
            
            [LPRout,pCP,scaled_prior] = accGlaze_fast(LLRin,mfit.H,0,'pCP',pIn);
            
        elseif strcmp(modeltype,'H_noise_IU')
            load([modelpath,modeltype,filesep,'fits',filesep,subj,'_fit.mat'])
            mfit.H = pm_fit(1);
            mfit.noise = pm_fit(2);
            mfit.IU = pm_fit(3);
            
            [LPRout,pCP,scaled_prior] = accGlaze_InconUp_fast(LLRin,mfit.H,mfit.IU,0,'pCP',pIn);
        end
        pCP = log(pCP./(1-pCP));  % logit-transform
        
        % Load eyetracker file
        load([loadpath,files(f).name]);
        pupil = zscore(data.pupil);  % z-scoring within-block
        times = data.times;
        
        % Aligning behavioural and eye-tracker data
        nsampstest=[];
        for t = 1:length(choices)
            nsampstest(t,1) = sum(~isnan(stimIn(t,:)));  % number of samples per trial
        end
        
        if strcmp(files(f).name,'t013_2_4.mat')  % special case where one ET trial has 12 samples somehow
            data.event = data.event([1:80 82:end],:);
            data.eventsmp = data.eventsmp([1:80 82:end],1:10);
            special_file = 1;
        elseif strcmp(files(f).name,'t029_1_7.mat')   % special case where one ET trial has one less sample than it should
            data.event = data.event([1:65 67:end],:);
            data.eventsmp = data.eventsmp([1:65 67:end],1:10);
            special_file = 1;
        elseif strcmp(files(f).name,'t040_1_2.mat')  % special case where one ET trial has 20 samples somehow
            data.event = data.event([1 3:end],:);
            data.eventsmp = data.eventsmp([1 3:end],1:10);
            special_file = 1;
        elseif strcmp(files(f).name,'t040_2_4.mat')   % special case where one ET trial has one less sample than it should
            data.event = data.event([1:50 52:end],:);
            data.eventsmp = data.eventsmp([1:50 52:end],1:10);
            special_file = 1;
        elseif strcmp(files(f).name,'t054_2_5.mat')   % special case where one ET trial has one less sample than it should
            data.event = data.event([1:35 37:end],:);
            data.eventsmp = data.eventsmp([1:35 37:end],1:10);
            special_file = 1;
        elseif strcmp(files(f).name,'t064_2_2.mat')   % special case where one ET trial has one less sample than it should
            data.event = data.event([1:43 45:end],:);
            data.eventsmp = data.eventsmp([1:43 45:end],1:10);
            special_file = 1;
        elseif strcmp(files(f).name,'t101_1_4.mat')  % special case where one ET trial has 15 samples somehow
            data.event = data.event([21 23:end],:);
            data.eventsmp = data.eventsmp([21 23:end],1:10);
            special_file = 1;
        else
            special_file = 0;
        end
        if special_file==1
            nsampstest = nsampstest(data.event(:,10));
            choices = choices(data.event(:,10));
            stimIn = stimIn(data.event(:,10),:);
            pIn = pIn(data.event(:,10),:);
            LLRin = LLRin(data.event(:,10),:);
            LPRout = LPRout(data.event(:,10),:);
            pCP = pCP(data.event(:,10),:);
            scaled_prior = scaled_prior(data.event(:,10),:);
        end
        
        if length(data.event(:,1))>length(choices)  % if # eye-tracker trials exceeds # behavioural trials
            t=1;
            while sum(nsampstest-data.event(t:(t+length(choices)-1),1))~=0
                t = t+1;
            end
            data.event = data.event(t:(t+length(choices)-1),:);
            data.eventsmp = data.eventsmp(t:(t+length(choices)-1),:);
        elseif length(data.event(:,1))<length(choices)  % if # behavioural trials exceeds # eye-tracker trials
            t=1;
            while sum(nsampstest(t:(t+length(data.event(:,1))-1),1)-data.event(:,1))~=0
                t = t+1;
            end
            nsampstest = nsampstest(t:(t+length(data.event(:,1))-1));
            choices = choices(t:(t+length(data.event(:,1))-1));
            stimIn = stimIn(t:(t+length(data.event(:,1))-1),:);
            pIn = pIn(t:(t+length(data.event(:,1))-1),:);
            LLRin = LLRin(t:(t+length(data.event(:,1))-1),:);
            LPRout = LPRout(t:(t+length(data.event(:,1))-1),:);
            pCP = pCP(t:(t+length(data.event(:,1))-1),:);
            scaled_prior = scaled_prior(t:(t+length(data.event(:,1))-1),:);
        end
        
        % Ensuring behavioural and pupil datasets contain same trials
        assert(sum(nsampstest-data.event(:,1))==0,'Sample mismatch for file %s.',files(f).name);
        
        % Downsampling EL data (speeds processing and aligns all datasets to same Fs if some were not recorded @ desired sampling rate)
        pupil = resample(pupil,newFs,data.fsample)';
        data.Xgaze = resample(data.Xgaze,newFs,data.fsample)';    % X-GAZE regressor
        data.Ygaze = resample(data.Ygaze,newFs,data.fsample)';    % Y-GAZE regressor
        times = (0:(length(pupil)-1))./newFs; data.times = times;  % manually creating new times vector
        
        data.event(:,2:4) = round(data.event(:,2:4).*(newFs/data.fsample));
        data.eventsmp = round(data.eventsmp.*(newFs/data.fsample));
        data.badsmp = unique(round(data.badsmp.*(newFs/data.fsample)));  % log of all bad samples that were previously interpolated
        data.badsmp(data.badsmp==0) = [];  % in case a sample index was rounded to zero
        
        data.fsample = newFs;  % replacing stored sampling rate in data structure
        
        % Calculate power spectral density of fixed-duration pupil signal for this block
        if length(pupil) >= data.fsample*9*60
            pupil4psd = pupil(data.fsample*10:data.fsample*9*60);
            [Pxx,F] = periodogram(pupil4psd,[],length(pupil4psd),data.fsample);
            pupilPSD = [pupilPSD Pxx(F>0.01 & F<6)];  % concatenating PSDs from all blocks for this subject
            if ~exist('PSDfreqs','var')
                PSDfreqs = F(F>0.01 & F<6);  % storing frequency vector
            end
        end
        
        % Initialize times vectors
        if f==1
            fulltimes = times(times<=diff(fullwin))+fullwin(1);
            samptimes = times(times<=diff(sampwin))+sampwin(1);
            fbtimes = times(times<=diff(fbwin))+fbwin(1);
            arttimes = times(times<=diff(artwin))+artwin(1);
        end
        
        % Initialize bad samples vector
        if size(pupil,1)>size(pupil,2), pupil=pupil'; end
        badsmps = zeros(size(pupil));
        badsmps(data.badsmp) = 1;   % all artifacts marked with 1
        edges1 = find(diff(badsmps)==1)+1;   % start of artifacts
        edges2 = find(diff(badsmps)==-1)+1;  % end of artifacts
        if badsmps(1)==1, edges1 = [1 edges1]; end
        if badsmps(end)==1, edges2(end+1) = length(badsmps); end
        for a = 1:length(edges1)
            if edges2(a)-edges1(a) > artlen*newFs
                badsmps(edges1(a):edges2(a)) = 2;  % long artifacts marked with 2
            end
        end
        
        % Isolating useable full-sequence trials
        ts=[];
        for t = 1:length(choices)
            if sum(isnan(stimIn(t,:)))==0 && choices(t)<2, ts(end+1) = t; end   %  && Behav(t,6)==0
        end
        
        % Collating useable single trials
        LLR_full = [LLR_full; LLRin(ts,:)];
        LPR_full = [LPR_full; LPRout(ts,:)];
        psi_full = [psi_full; scaled_prior(ts,:)];
        pCP_full = [pCP_full; pCP(ts,:)];
        choices_full = [choices_full; choices(ts)];
        
        % Filter
        [bfilt,afilt] = butter(3, freqs(1)*2/data.fsample, 'high');   % hi-pass
        pupil = filtfilt(bfilt,afilt, pupil);
        
        [bfilt,afilt] = butter(3, freqs(2)*2/data.fsample, 'low');   % lo-pass
        pupil = filtfilt(bfilt,afilt, pupil);
        
        pupil = zscore(pupil);
        rawpupil = pupil;
        if first_deriv  % taking 1st derivative of pupil signal if desired, in z/s
            pupil = diff(pupil).*data.fsample;
        end
        
        % Loop through trials
        for t = ts
            % Pull full trial response
            smp1 = find(times>=times(data.event(t,2))+fullwin(1),1,'first');
            if smp1+length(fulltimes)-1 <= length(pupil)  % making sure epoch lies within range
                if ~first_deriv
                    dil_full(end+1,:) = pupil(smp1:smp1+length(fulltimes)-1)-mean(pupil(times>=times(data.event(t,2))+basewin(1) & times<=times(data.event(t,2))+basewin(2)));
                else dil_full(end+1,:) = pupil(smp1:smp1+length(fulltimes)-1);  % no baseline if first derivative
                end
                % Pull individual sample response
                csmp = size(dil_samp,1)+1;
                for smp = 1:size(data.eventsmp,2)
                    smpbase_full(csmp,smp) = mean(rawpupil(times>=times(data.eventsmp(t,smp))+basewin(1) & times<=times(data.eventsmp(t,smp))+basewin(2)));
                end
                if ~lin_detrend
                    if ~singlesamp_base
                        cbase = mean(pupil(times>=times(data.event(t,2))+basewin(1) & times<=times(data.event(t,2))+basewin(2)));
                    end
                    for smp = 1:size(data.eventsmp,2)
                        smp1 = find(times>=times(data.eventsmp(t,smp))+sampwin(1),1,'first');
                        if ~first_deriv
                            if singlesamp_base
                                dil_samp(csmp,smp,:) = pupil(smp1:smp1+length(samptimes)-1)-mean(pupil(times>=times(data.eventsmp(t,smp))+basewin(1) & times<=times(data.eventsmp(t,smp))+basewin(2)));
                            else dil_samp(csmp,smp,:) = pupil(smp1:smp1+length(samptimes)-1)-cbase;
                            end
                        else dil_samp(csmp,smp,:) = pupil(smp1:smp1+length(samptimes)-1);
                        end
                    end
                else
                    startsmp = find(times>=times(data.event(t,2))+sampwin(1),1,'first'); endsmp = find(times>=times(data.eventsmp(t,end))+sampwin(2),1,'first')+2;
                    ctimes = times(startsmp:endsmp)-min(times(startsmp:endsmp))+sampwin(1);
                    cwave = pupil(startsmp:endsmp); cevents = data.eventsmp(t,:)-startsmp+1;
                    m = regstats(cwave,1:length(cwave),'linear',{'r'});
                    cresids = m.r;
                    if ~singlesamp_base
                        cbase = mean(cresids(ctimes>=basewin(1) & ctimes<=basewin(2)));
                    end
                    for smp = 1:length(cevents)
                        smp1 = find(ctimes>=ctimes(cevents(smp))+sampwin(1),1,'first');
                        if singlesamp_base
                            dil_samp(csmp,smp,:) = cresids(smp1:smp1+length(samptimes)-1)-mean(cresidsL(ctimes>=ctimes(cevents(smp))+basewin(1) & ctimes<=ctimes(cevents(smp))+basewin(2)));
                        else dil_samp(csmp,smp,:) = cresids(smp1:smp1+length(samptimes)-1)-cbase;
                        end
                    end
                end
                for smp = 1:size(data.eventsmp,2)
                    smp1 = find(times>=times(data.eventsmp(t,smp))+sampwin(1),1,'first');
                    X_samp(csmp,smp,:) = data.Xgaze(smp1:smp1+length(samptimes)-1);
                    Y_samp(csmp,smp,:) = data.Ygaze(smp1:smp1+length(samptimes)-1);
                end
                
                % Pull feedback response
                smp1 = find(times>=times(data.event(t,4))+fbwin(1),1,'first');
                if smp1+length(fbtimes)-1 <= length(pupil)  % making sure epoch lies within range
                    if data.event(t,6)==1
                        dil_fbC(end+1,:) = pupil(smp1:smp1+length(fbtimes)-1)-mean(pupil(times>=times(data.event(t,4))+basewin(1) & times<=times(data.event(t,4))+basewin(2)));
                    elseif data.event(t,6)==0
                        dil_fbE(end+1,:) = pupil(smp1:smp1+length(fbtimes)-1)-mean(pupil(times>=times(data.event(t,4))+basewin(1) & times<=times(data.event(t,4))+basewin(2)));
                    end
                end
                
                % Test for artifacts
                smp1 = find(times>=times(data.event(t,2))+artwin(1),1,'first');
                if ~isempty(find(badsmps(smp1:(smp1+length(arttimes)-1))==2, 1)) || length(find(badsmps(smp1:(smp1+length(arttimes)-1))>0))/length(arttimes)>artp
                    bad_ts(end+1) = 1;
                else bad_ts(end+1) = 0;
                end
            else  % just exclude final trial of block if recording was cut off prematurely
                bad_ts(end+1) = 1;
                
                dil_full(end+1,:) = nan;
                dil_samp(end+1,:,:) = nan;
                X_samp(end+1,:,:) = nan;
                Y_samp(end+1,:,:) = nan;
            end
        end
        
    catch
        errfiles{end+1,1} = files(f).name;
    end
end

% Keeping only clean trials
dil_full = dil_full(bad_ts==0,:);
dil_samp = dil_samp(bad_ts==0,:,:);
X_samp = X_samp(bad_ts==0,:,:);
Y_samp = Y_samp(bad_ts==0,:,:);
smpbase_full = smpbase_full(bad_ts==0,:);
pCP_full = pCP_full(bad_ts==0,:);
psi_full = psi_full(bad_ts==0,:);
LLR_full = LLR_full(bad_ts==0,:);
choices_full = choices_full(bad_ts==0);

% Collating trial-averaged responses
dil_full_av = mean(dil_full,1);
dil_samp_av = mean(dil_samp,1);
dil_fbC_av = mean(dil_fbC,1);
dil_fbEav = mean(dil_fbE,1);
ntrials = size(dil_full,1);
p_artifact = 1-(size(dil_full,1)/length(bad_ts));

% Running sample-wise regressions of pupil dilation onto surprise & storing betas & residuals
pupil_r=[]; pupil_r_ext=[];
for i = 1:length(samptimes)
    ts = find_inliers([squeeze(dil_samp(:,:,i)) pCP_full(:,2:end)],outliercut);
    for smp = 1:maxsamps-1
        if quad_reg   % include quadratic term?
            RsP = [nanzscore(pCP_full(:,smp+1)) nanzscore(pCP_full(:,smp+1)).^2];
        else  RsP = nanzscore(pCP_full(:,smp+1)); end
        if smp_regout_lag==0 || smp+1-smp_regout_lag<=1
            smpmodel = regstats(nanzscore(dil_samp(ts,smp+1,i)),[RsP(ts,:) X_samp(ts,smp+1,i) Y_samp(ts,smp+1,i) zscore(smpbase_full(:,smp+1))],'linear',{'beta','r','tstat'});
        else % regressing out surprise on directly preceding samples if desired
            smpmodel = regstats(nanzscore(dil_samp(ts,smp+1,i)),[RsP(ts,:) nanzscore(pCP_full(ts,max([2 smp+1-smp_regout_lag]):smp)) X_samp(ts,smp+1,i) Y_samp(ts,smp+1,i) zscore(smpbase_full(:,smp+1))],'linear',{'beta','r','tstat'});
        end
        dil_pCP_B(smp,i) = smpmodel.beta(2);
        pupil_r(:,smp,i) = nanzscore(smpmodel.r);
    end
    ts = find_inliers([squeeze(dil_samp(:,:,i)) pCP_full(:,2:end) -abs(psi_full(:,2:end))],outliercut);
    for smp = 1:maxsamps-1
        RsP = [nanzscore(pCP_full(:,smp+1)) nanzscore(-abs(psi_full(:,smp+1))) abs(LLR_full(:,smp+1))];
        if smp_regout_lag==0 || smp+1-smp_regout_lag<=1
            smpmodel = regstats(nanzscore(dil_samp(ts,smp+1,i)),[RsP(ts,:) X_samp(ts,smp+1,i) Y_samp(ts,smp+1,i) zscore(smpbase_full(:,smp+1))],'linear',{'beta','r','tstat'});
        else % regressing out surprise on directly preceding samples if desired
            smpmodel = regstats(nanzscore(dil_samp(ts,smp+1,i)),[RsP(ts,:) nanzscore(pCP_full(ts,max([2 smp+1-smp_regout_lag]):smp)) nanzscore(-abs(psi_full(ts,max([2 smp+1-smp_regout_lag]):smp))) nanzscore(LLR_full(ts,max([2 smp+1-smp_regout_lag]):smp)) X_samp(ts,smp+1,i) Y_samp(ts,smp+1,i) zscore(smpbase_full(:,smp+1))],'linear',{'beta','r','tstat'});
        end
        dil_pCP_Bext(smp,i,1:3) = smpmodel.beta(2:4);
        pupil_r_ext(:,smp,i) = nanzscore(smpmodel.r);
    end
end

% compute kernels
[B,~,stats] = glmfit([nanzscore(LLR_full,0,1) nanzscore(LLR_full(:,2:end).*nanzscore(pCP_full(:,2:end)),0,1) nanzscore(LLR_full(:,2:end).*nanzscore(-abs(psi_full(:,2:end))),0,1)],[choices_full ones(length(choices_full),1)],'binomial');
ssB_regU = B(2:maxsamps+1);
ssB_surpU = B(maxsamps+2:end-length(2:maxsamps));
ssB_uncertU = B(end-length(2:maxsamps)+1:end);

% compute PPIs without surprise/uncertainty modulations
for i = 1:length(samptimes)  % PPIs with LLR*surprise interaction terms also included in models
    B = glmfit([nanzscore(LLR_full,0,1) nanzscore(LLR_full(:,2:end).*squeeze(pupil_r(:,:,i)),0,1)],[choices_full ones(length(choices_full),1)],'binomial');
    dil_PPIsimple_B(:,i) = B(end-length(2:maxsamps)+1:end);
    
    B = glmfit([nanzscore(LLR_full,0,1) nanzscore(LLR_full(:,2:end).*squeeze(pupil_r_ext(:,:,i)),0,1)],[choices_full ones(length(choices_full),1)],'binomial');
    dil_PPIsimple_Bext(:,i) = B(end-length(2:maxsamps)+1:end);
end

% compute PPIs
for i = 1:length(samptimes)  % PPIs with LLR*surprise interaction terms also included in models
    B = glmfit([nanzscore(LLR_full,0,1) nanzscore(LLR_full(:,2:end).*nanzscore(pCP_full(:,2:end)),0,1) nanzscore(LLR_full(:,2:end).*nanzscore(-abs(psi_full(:,2:end))),0,1) nanzscore(LLR_full(:,2:end).*squeeze(pupil_r(:,:,i)),0,1)],[choices_full ones(length(choices_full),1)],'binomial');
    dil_PPI_B(:,i) = B(end-length(2:maxsamps)+1:end);
    
    B = glmfit([nanzscore(LLR_full,0,1) nanzscore(LLR_full(:,2:end).*nanzscore(pCP_full(:,2:end)),0,1) nanzscore(LLR_full(:,2:end).*nanzscore(-abs(psi_full(:,2:end))),0,1) nanzscore(LLR_full(:,2:end).*squeeze(pupil_r_ext(:,:,i)),0,1)],[choices_full ones(length(choices_full),1)],'binomial');
    dil_PPI_Bext(:,i) = B(end-length(2:maxsamps)+1:end);
end

% compute PPIs with lasso regression
for i = 1:length(samptimes)  % PPIs with LLR*surprise interaction terms also included in models
    B = lassoglm([nanzscore(LLR_full,0,1) nanzscore(LLR_full(:,2:end).*nanzscore(pCP_full(:,2:end)),0,1) nanzscore(LLR_full(:,2:end).*nanzscore(-abs(psi_full(:,2:end))),0,1) nanzscore(LLR_full(:,2:end).*squeeze(pupil_r(:,:,i)),0,1)],[choices_full],'binomial','Lambda',0.002,'CV',10);
    dil_PPI_B_lasso(:,i) = B(end-length(2:maxsamps)+1:end);
    
    B = lassoglm([nanzscore(LLR_full,0,1) nanzscore(LLR_full(:,2:end).*nanzscore(pCP_full(:,2:end)),0,1) nanzscore(LLR_full(:,2:end).*nanzscore(-abs(psi_full(:,2:end))),0,1) nanzscore(LLR_full(:,2:end).*squeeze(pupil_r_ext(:,:,i)),0,1)],[choices_full],'binomial','Lambda',0.002,'CV',10);
    dil_PPI_Bext_lasso(:,i) = B(end-length(2:maxsamps)+1:end);
end

% save output
if first_deriv
    savename = [savepath,subj,'_d1_',modeltype,'.mat'];
else savename = [savepath,subj,'_raw_',modeltype,'.mat'];
end

save(savename,'dil_full_av','dil_samp_av','dil_fbC_av','dil_fbEav','ntrials','p_artifact',...
    'dil_pCP_B','dil_pCP_Bext','ssB_regU','ssB_surpU','ssB_uncertU','dil_PPIsimple_B','dil_PPIsimple_Bext','dil_PPI_B','dil_PPI_Bext','dil_PPI_B_lasso','dil_PPI_Bext_lasso',...
    'fulltimes','samptimes','fbtimes','modeltype','mfit')

if ~isempty(errfiles), save([savepath,'bad_subj',filesep,subj,'.mat'],'errfiles'), end
