% Script for two-alternative sequential evidence samples drawn from one
% of 2 distributions with same std but different means, and mean can switch
% at fixed hazard rate.

% muE = generative distribution means
% sigmaE = generative distribution standard deviations
% rangeE = generative distribution range for truncation
% H = hazard rate of distribution switches
% ntrials = number of simulated trials
% nsamps = number evidence samples per trial

function LLRin = Gen_LLR(muE,sigmaE,rangeE,H,pshort,nsamps,ntrials)

% Generate sequences of distribution switch positions
pswitch = [ones(ntrials,1) rand([ntrials,nsamps-1])];  % first draw uniformly distributed random probabilities b/w 0 and 1 that will determine switch positions (first sample is never a switch)
pswitch(pswitch>H) = 0; pswitch(pswitch~=0) = 1;  % binarize matrix to mark only switch positions

% Generate sequences of which distributions will be drawn from at which times (1=left, 2=right)
distseqs = zeros(ntrials,nsamps); dists = [1 2];
for t = 1:ntrials
    if t<=ntrials/2, cdist = 1; else cdist = 2; end  % making sure each distribution is starting distribution an equal number of times
    s = 0;
    while s<nsamps
        s = s+1;
        if ~pswitch(t,s), distseqs(t,s) = cdist;  % if switch in distribution has not occured
        else cdist = dists(dists~=cdist);  % if switch in distribution has occured
            distseqs(t,s) = cdist;
        end
    end
end

% Generate actual evidence samples from distribution sequences
stimIn = distseqs;
stimIn(distseqs==1) = round(muE(1)+(randn(size(stimIn(distseqs==1))).*sigmaE(1)));
stimIn(distseqs==2) = round(muE(2)+(randn(size(stimIn(distseqs==2))).*sigmaE(2)));
stimIn(stimIn<rangeE(1)) = rangeE(1); stimIn(stimIn>rangeE(2)) = rangeE(2);  % in case drawn values exceed range limits

% Shuffle trial order
shuforder = randperm(ntrials);
pswitch = pswitch(shuforder,:);
stimIn = stimIn(shuforder,:);
distseqs = distseqs(shuforder,:);

% Compute LLRs
LLRin = log(normpdf(stimIn,muE(2),sigmaE(2))./normpdf(stimIn,muE(1),sigmaE(1)));

% Store generative distributions @ sequence end (i.e. correct choices; 1=left, 2=right)
fdists = distseqs(:,end);

% Truncate sequence length of random selection of trials
tlengths = randsample(2:(nsamps-1),round(ntrials*pshort),true);  % randomly draw lengths of truncated sequences
ttrials = sort(randsample(ntrials,round(ntrials*pshort),false))';  % randomly draw trials to be truncated
for t = 1:length(ttrials)
    pswitch(ttrials(t),tlengths(t)+1:end) = nan;
    stimIn(ttrials(t),tlengths(t)+1:end) = nan;
    distseqs(ttrials(t),tlengths(t)+1:end) = nan;
    fdists(ttrials(t)) = distseqs(ttrials(t),tlengths(t));
    LLRin(ttrials(t),tlengths(t)+1:end) = nan;
end