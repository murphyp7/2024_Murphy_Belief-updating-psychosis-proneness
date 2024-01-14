function [newpupil, newXgaze, newYgaze, newblinksmp, nanIdx, p, h, bestchan] = blink_interpolate_1chan(dat, inFile, plotme)
% interpolates blinks and missing data
% Anne Urai, 2016
% Peter Murphy, 2016 - added adaptive interpolation windows

% interpolation parameters
p=[];
p.padding1         = [-0.150 0.150]; % padding before/after EL-defined blinks for initial rough interpolation
p.wfilt            = 11;             % width (in samples) of Hanning window used for lo-pass filtering

p.blinkthres       = 1.5;            % diameter threshold (in mm) for crudely defining blinks in first pass
p.diffthresh       = 0.75;           % z-scored PUPIL derivative threshold that samples before/after a bad period must consistently satisfy   (default = 0.75 for quick blinkers, 0.55 for slow)
p.gazethresh       = 6.0;            % z-scored GAZE derivative threshold that samples before/after a bad period must consistently satisfy   (default = 6.0 - gaze tends to be noisier/higher freq than pupil)
p.t_pre            = 0.05;           % window (s) before a bad period in which all samples must satisfy p.diffthresh   (default = 0.05)
p.t_post           = 0.06;           % window (s) after a bad period in which all samples must satisfy p.diffthresh   (default = 0.06)
p.tmax             = [0.25 0.5];     % max time (s) [start end] points of interpolation window can go from identified bad samples   (default = 0.25 for quick blinkers, 0.5 for slow)
p.diffthresh2      = 5;              % z-scored derivative threshold for identification of bad periods not covered by Eyelink   (default = 3.5)

p.coalesce1        = 0.250;          % merge 2 blinks into 1 if they are below this distance (in s) apart (default = 0.250)
p.coalesce2        = 0.500;          % merge 2 bad periods (IDd from step 2) into 1 if they are below this distance (in s) apart (default = 0.500)

% pick pupil channel (left or right) to use
bestchan = find([var(dat.pupilL) var(dat.pupilR)]==min([var(dat.pupilL) var(dat.pupilR)]));
chanlabels = {'left','right'};

% make copies of some stuff
if bestchan==1
    pupilcopy = dat.pupilL;
else
    pupilcopy = dat.pupilR;
end
Xgazecopy = dat.Xgaze;
Ygazecopy = dat.Ygaze;
nanIdx = [];

% plot if unspecified
if ~exist('plotme', 'var'); plotme = true; end


% ====================================================== %
% STEP 1: INTERPOLATE COARSELY-DEFINED BLINKS USING WIDE FIXED WINDOW
% ====================================================== %

% plot raw time-series
if plotme
    h = figure('Position', [80, 100, 1500, 840]);
    sp1 = subplot(511); hold on
    plot(dat.times,dat.pupilL, 'color', [0.5 0.5 1]);
    plot(dat.times,dat.pupilR, 'color', [1 0.5 0.5]);
    axis tight; box off; ylabel('Raw'); title([inFile,'; best chan=',chanlabels{bestchan}])
    set(gca, 'xtick', []);
else h=[];
end

% get rough blink epochs
blinks = find(pupilcopy<p.blinkthres);
blinksmp = [];

edgesL = find(diff(blinks)~=1);
for b=1:length(edgesL)
    if b==1
        blinksmp(end+1,1:2) = [blinks(1) blinks(edgesL(b))];
    else blinksmp(end+1,1:2) = [blinks(edgesL(b-1)+1) blinks(edgesL(b))];
        if b==length(edgesL)
            blinksmp(end+1,1:2) = [blinks(edgesL(b)+1) blinks(end)];
        end
    end
end

% merge consecutive blinks into 1 if they are X ms together
win          = hanning(p.wfilt);

if ~isempty(blinksmp)
    [~,si] = sort(blinksmp(:,1));
    blinksmp = blinksmp(si,:);
    
%     cblinksmp = blinksmp(1,:);
%     for b = 1:size(blinksmp,1)-1,
%         if (blinksmp(b+1,1) - cblinksmp(end,2) < p.coalesce1*dat.fsample)
%             cblinksmp(end,2) = blinksmp(b+1,2);
%         else
%             cblinksmp(end+1,:) = blinksmp(b+1,:);
%         end
%     end
%     blinksmp = cblinksmp; clear cblinksmp
        
    % pad the blinks
    padblinksmp(:,1) = round(blinksmp(:,1) + p.padding1(1) * dat.fsample);
    padblinksmp(:,2) = round(blinksmp(:,2) + p.padding1(2) * dat.fsample);
    
    % avoid idx outside range
    if any(padblinksmp(:) < 1), padblinksmp(padblinksmp < 1) = 1; end
    if any(padblinksmp(:) > length(dat.pupilL)), padblinksmp(padblinksmp > length(dat.pupilL)) = length(dat.pupilL); end
    
    % interpolate
    [pupilcopy,~] = interp_nans(pupilcopy,padblinksmp);
    [Xgazecopy2,~] = interp_nans(Xgazecopy,padblinksmp);
    [Ygazecopy2,~] = interp_nans(Ygazecopy,padblinksmp);
    
    % check to make sure all nans have been dealt with
    assert(~any(isnan(pupilcopy)));
    
    % low-pass filter so later-derived threshold is on same scale as signal it's being applied to
    pupilcopy    = filter2(win.',pupilcopy,'same');
    
    
    % ====================================================== %
    % STEP 2: USE DERIVATIVE OF STEP 1 OUTPUT TO DEFINE INTERPOLATION WINDOWS
    % ====================================================== %
    
    % low-pass filter original time-series (otherwise particularly noisy measurements will yield very bad results)
    if bestchan==1
        pupilsmooth  = filter2(win.',dat.pupilL,'same');
    else
        pupilsmooth  = filter2(win.',dat.pupilR,'same');
    end
    
    % calculate raw-unit derivative threshold from roughly interpolated time-series
    raw_dthresh = p.diffthresh*std(diff(pupilcopy(p.wfilt:end-p.wfilt+1)));
    raw_Xthresh = p.gazethresh*std(diff(Xgazecopy2(p.wfilt:end-p.wfilt+1)));
    raw_Ythresh = p.gazethresh*std(diff(Ygazecopy2(p.wfilt:end-p.wfilt+1)));
    
    % use threshold to find window start/end points
    padblinksmp = [];
    for b = 1:size(blinksmp,1)
        s1 = blinksmp(b,1);   % starting sample
        s2 = blinksmp(b,2);   % ending sample
        if s1-round(p.t_pre*dat.fsample)<p.wfilt  % if this starting sample is very close to start of timeseries
            s1 = 1;  % just take first sample as starting point
        else  % otherwise, searching for appropriate starting sample by applying threhsold to pupil derivative
            %while max(abs(diff(pupilsmooth((s1-round(p.t_pre*dat.fsample)):s1-1))))>raw_dthresh && s1-round(p.t_pre*dat.fsample)>p.wfilt && (s1-blinksmp(b,1))>-p.tmax*dat.fsample
            while (max(abs(diff(pupilsmooth((s1-round(p.t_pre*dat.fsample)):s1-1))))>raw_dthresh ||...
                    max(abs(diff(dat.Xgaze((s1-round(p.t_pre*dat.fsample)):s1-1))))>raw_Xthresh ||...
                    max(abs(diff(dat.Ygaze((s1-round(p.t_pre*dat.fsample)):s1-1))))>raw_Ythresh)...
                    && s1-round(p.t_pre*dat.fsample)>p.wfilt && (s1-blinksmp(b,1))>-p.tmax(1)*dat.fsample
                s1 = s1-1;
            end
            if s1-(p.t_pre*dat.fsample) == p.wfilt
                s1 = 1;
            end
        end
        if s2+(p.t_post*dat.fsample)>length(pupilsmooth)-p.wfilt+1  % if this ending sample is very close to end of timeseries
            s2 = length(pupilsmooth);  % just take last sample as ending point
        else  % otherwise, searching for appropriate starting sample by applying threhsold to pupil derivative
            %while max(abs(diff(pupilsmooth((s2+2:s2+round(p.t_post*dat.fsample))))))>raw_dthresh && s2+round(p.t_post*dat.fsample)<length(pupilsmooth)-p.wfilt+1 && (s2-blinksmp(b,2))<p.tmax*dat.fsample
            while (max(abs(diff(pupilsmooth((s2+2:s2+round(p.t_post*dat.fsample))))))>raw_dthresh ||...
                    max(abs(diff(dat.Xgaze((s2+2:s2+round(p.t_post*dat.fsample))))))>raw_Xthresh ||...
                    max(abs(diff(dat.Ygaze((s2+2:s2+round(p.t_post*dat.fsample))))))>raw_Ythresh)...
                    && s2+round(p.t_post*dat.fsample)<length(pupilsmooth)-p.wfilt+1 && (s2-blinksmp(b,2))<p.tmax(2)*dat.fsample
                s2 = s2+1;
            end
            if s2+round(p.t_post*dat.fsample) == length(pupilsmooth)-p.wfilt+1
                s2 = length(pupilsmooth);
            end
        end
        padblinksmp(b,1:2) = [s1 s2];
    end
    
    % interpolate
    [dat.pupilL,nanIdx] = interp_nans(dat.pupilL,padblinksmp);
    [dat.pupilR,~] = interp_nans(dat.pupilR,padblinksmp);
    [dat.Xgaze,~] = interp_nans(dat.Xgaze,padblinksmp);
    [dat.Ygaze,~] = interp_nans(dat.Ygaze,padblinksmp);
    
    % check to make sure all nans have been dealt with
    assert(~any(isnan(dat.pupilL)));
    assert(~any(isnan(dat.pupilR)));
    
    % plot initial interpolation pass
    if plotme
        plot(dat.times, dat.pupilL, 'b');
        plot(dat.times, dat.pupilR, 'r');
        axis tight; box off; ylabel('Interp');
        set(gca, 'xtick', []);
    end
end

% ====================================================== %
% STEP 3: USE DERIVATIVE OF STEP 2 OUTPUT TO IDENTIFY REMAINING BAD SAMPLES
% ====================================================== %

% low-pass filtering
if bestchan==1
    pupilsmooth  = filter2(win.',dat.pupilL,'same');
else
    pupilsmooth  = filter2(win.',dat.pupilR,'same');
end

% identify periods with derivative that exceeds harsh threshold
pupildiff = zscore(diff(pupilsmooth));   % calculate derivative of filtered pupil signal
Xgazediff = zscore(diff(dat.Xgaze));   % calculate derivative of X-gaze position
Ygazediff = zscore(diff(dat.Ygaze));   % calculate derivative of Y-gaze position
badsmp = find(abs(pupildiff) > p.diffthresh2*std(pupildiff));  % get positions of all outlying data points from chosen pupil channel
badsmp = badsmp(badsmp>=p.wfilt & badsmp<=length(pupilsmooth)-p.wfilt+1);  % chopping off first and last n samples, which are contaminated by filtering

if ~isempty(badsmp)
    
    badpts = [badsmp(1) badsmp(find(abs(diff(badsmp))>1)+1)]';  % get samples where periods of outlying data points begin
    badpts(:,2) = [badsmp(abs(diff(badsmp))>1) badsmp(end)]';  % get samples where periods of outlying data points end
    
    % merge 2 windows into 1 if they are X ms together (since there will usually be a gap between down- and up-turns in derivative)
    cbadpts = badpts(1,:);
    for b = 1:size(badpts,1)-1,
        if badpts(b+1,1) - cbadpts(end,2) < p.coalesce2 * dat.fsample,
            cbadpts(end,2) = badpts(b+1,2);
        else
            cbadpts(end+1,:) = badpts(b+1,:);
        end
    end
    badpts = cbadpts; clear cbadpts
    
    %     for b = 1:size(badpts, 1)-1,
    %         if badpts(b+1, 1) - badpts(b, 2) < p.coalesce2 * dat.fsample,
    %             badpts(b, 2) = badpts(b+1, 2);
    %             badpts(b+1, :) = nan;
    %         end
    %     end
    %     badpts(isnan(nanmean(badpts, 2)), :) = []; % remove those duplicates
    
    % calculate derivative thresholds for this round
    raw_dthresh = p.diffthresh*std(abs(diff(pupilsmooth(p.wfilt:end-p.wfilt+1))));
    raw_Xthresh = p.gazethresh*std(abs(diff(dat.Xgaze(p.wfilt:end-p.wfilt+1))));
    raw_Ythresh = p.gazethresh*std(abs(diff(dat.Ygaze(p.wfilt:end-p.wfilt+1))));
    
    % use threshold to find window start/end points
    padblinksmp = [];
    for b = 1:size(badpts,1)
        s1 = badpts(b,1);   % starting sample
        s2 = badpts(b,2);   % ending sample
        if s1-round(p.t_pre*dat.fsample)<p.wfilt  % if this starting sample is very close to start of timeseries
            s1 = 1;  % just take first sample as starting point
        else  % otherwise, searching for appropriate starting sample by applying threhsold to pupil derivative
            while (max(abs(diff(pupilsmooth((s1-round(p.t_pre*dat.fsample)):s1-1))))>raw_dthresh ||...
                    max(abs(diff(dat.Xgaze((s1-round(p.t_pre*dat.fsample)):s1-1))))>raw_Xthresh ||...
                    max(abs(diff(dat.Ygaze((s1-round(p.t_pre*dat.fsample)):s1-1))))>raw_Ythresh)...
                    && s1-round(p.t_pre*dat.fsample)>p.wfilt && (s1-badpts(b,1))>-p.tmax(1)*dat.fsample
                s1 = s1-1;
            end
            if s1-(p.t_pre*dat.fsample) == p.wfilt
                s1 = 1;
            end
        end
        if s2+(p.t_post*dat.fsample)>length(pupilsmooth)-p.wfilt+1  % if this ending sample is very close to end of timeseries
            s2 = length(pupilsmooth);  % just take last sample as ending point
        else  % otherwise, searching for appropriate starting sample by applying threhsold to pupil derivative
            while (max(abs(diff(pupilsmooth((s2+2:s2+round(p.t_post*dat.fsample))))))>raw_dthresh ||...
                    max(abs(diff(dat.Xgaze((s2+2:s2+round(p.t_post*dat.fsample))))))>raw_Xthresh ||...
                    max(abs(diff(dat.Ygaze((s2+2:s2+round(p.t_post*dat.fsample))))))>raw_Ythresh)...
                    && s2+round(p.t_post*dat.fsample)<length(pupilsmooth)-p.wfilt+1 && (s2-badpts(b,2))<p.tmax(2)*dat.fsample
                s2 = s2+1;
            end
            if s2+round(p.t_post*dat.fsample) == length(pupilsmooth)-p.wfilt+1
                s2 = length(pupilsmooth);
            end
        end
        padblinksmp(b,1:2) = [s1 s2];
    end
    
    % interpolate
    if bestchan==1
        old_dat = dat.pupilL;
    else
        old_dat = dat.pupilR;
    end
    [dat.pupilL,new_nans] = interp_nans(dat.pupilL,padblinksmp);
    [dat.pupilR,~] = interp_nans(dat.pupilR,padblinksmp);
    [dat.Xgaze,~] = interp_nans(dat.Xgaze,padblinksmp);
    [dat.Ygaze,~] = interp_nans(dat.Ygaze,padblinksmp);
    nanIdx(end+1:end+length(new_nans)) = new_nans;
    
    % check to make sure all nans have been dealt with
    assert(~any(isnan(dat.pupilL)));
    assert(~any(isnan(dat.pupilR)));
    
    % Plotting derivative and final interpolated timeseries
    if plotme, sp2 = subplot(512); hold on;
        plot(dat.times(11:end-10), Xgazediff(10:end-10), 'color', [0 1 0]);
        plot(dat.times(11:end-10), Ygazediff(10:end-10), 'color', [1 1 0]);
        plot(dat.times(11:end-10), pupildiff(10:end-10),'color', [0.15 0.15 0.15]);  % trimming few edge samples because these can be very extreme
        plot(dat.times(new_nans), zeros(1,length(new_nans)), '.');
        box off; ylabel('Derivative');
        set(gca, 'xtick', []); ylim([-max(abs(pupildiff(11:end-10)))*1.02 max(abs(pupildiff(11:end-10)))*1.02]);
        
        sp3 = subplot(513); hold on;
        plot(dat.times, old_dat, 'color', [0.5 0.5 0.5]);
        if bestchan==1
            plot(dat.times, dat.pupilL, 'b');
            ylim([min(dat.pupilL)*0.9 max(dat.pupilL)*1.08]);
        else
            plot(dat.times, dat.pupilR, 'r');
            ylim([min(dat.pupilR)*0.9 max(dat.pupilR)*1.08]);
        end
        box off; ylabel('Clean');
        
        sp4 = subplot(514); hold on;
        plot(dat.times, Xgazecopy, 'color', [0.5 0.5 0.5]);
        plot(dat.times, dat.Xgaze, 'color', [0 1 0]);
        if min(dat.Xgaze)~=max(dat.Xgaze)
            ylim([min(dat.Xgaze)*0.8 max(dat.Xgaze)*1.2]);
        else ylim([-1 1].*(min(dat.Xgaze)+1))
        end
        box off; ylabel('X-gaze');
        
        sp5 = subplot(515); hold on;
        plot(dat.times, Ygazecopy, 'color', [0.5 0.5 0.5]);
        plot(dat.times, dat.Ygaze, 'color', [1 1 0]);
        if min(dat.Ygaze)~=max(dat.Ygaze)
            ylim([min(dat.Ygaze)*0.8 max(dat.Ygaze)*1.2]);
        else ylim([-1 1].*(min(dat.Ygaze)+1))
        end
        box off; ylabel('Y-gaze');
        
        try
            linkaxes([sp1 sp2 sp3 sp4 sp5], 'x');
            set([sp1 sp2 sp3 sp4 sp5], 'tickdir', 'out');
        end
        xlim([dat.times(1) dat.times(end)]);
    end
    
    newblinksmp = padblinksmp;
else newblinksmp = [];
end

% specify output
if bestchan==1
    newpupil = dat.pupilL;
else
    newpupil = dat.pupilR;
end
newXgaze = dat.Xgaze;
newYgaze = dat.Ygaze;
bestchan = chanlabels{bestchan};



