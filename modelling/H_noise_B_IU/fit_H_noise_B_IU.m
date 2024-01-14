% Uses Glaze model to generate finite number of observed choices for a
% given combination of H (hazard rate), B (gain) and noise parameters, and
% attempts to recover parameters using particle swarm optimization

% ntrials = number of trials to simulate
% i = iteration number
% pm(1:3) = [H, B, noise]

function [pm_fit,err] = fit_H_noise_B_IU(subj)

% Set global variables to pass to PSO algorithm
global LLRin choices nsamps

% Pick subject
allsubj = {'t010';'t011';'t012';'t013';'t014';'t015';'t016';'t017';'t018';'t019';...
           't020';'t021';'t022';'t023';'t024';'t025';'t026';'t027';'t028';'t029';...
           't030';'t031';'t032';'t033';'t034';'t035';'t036';'t037';'t038';'t039';...
           't040';'t041';'t043';'t044';'t042';'t045';'t046';'t047';'t048';'t049';...
           't050';'t051';'t052';'t053';'t054';'t055';'t056';'t057';'t058';'t059';...
           't060';'t061';'t062';'t063';'t064';'t065';'t066';'t067';'t068';'t069';...
           't070';'t071';'t072';'t073';'t074';'t075';'t076';'t077';'t078';'t079';...
           't080';'t081';'t082';'t083';'t084';'t085';'t086';'t087';'t088';'t089';...
           't090';'t091';'t092';'t093';'t094';'t095';'t096';'t097';'t098';'t099';...
           't100';'t101';'t102';'t103';'t104';'t105'};
subj = allsubj{subj}

% Path stuff
addpath /mnt/homes/home024/pmurphy/Surprise_scz/modelling/H_noise_B_IU/
addpath(genpath('/mnt/homes/home024/pmurphy/Surprise_accumulation/Simulations/particle_swarm_Glaze'))
loadpath = '/mnt/homes/home028/gmonov/SCZ/Data/decision_making/';
savepath = '/mnt/homes/home024/pmurphy/Surprise_scz/modelling/H_noise_B_IU/fits/';

% Define parameter estimation settings
range_p.H = [0.0001 0.5];  % parameter bounds
range_p.noise = [0.0001 10];
range_p.B = [0.01 3.5];
range_p.incon = [0.01 10];

mv = [0.005;...    % maximum particle velocities
    0.2;...
    0.04;...
    0.04]';

seeds.H = [0.1 0.04];  % good seed distributions for a selection of particles - [mean sd]
seeds.noise = [1.5 0.5];
seeds.B = [1 0.3];
seeds.incon = [1 0.5];

% Seed random number generator
seed = round(sum(100*clock)); %never the same seed
rand('state', seed);

% Load data for current subject
fprintf('Pulling behavioural data for subject %s...\n',subj)

fsess = dir([loadpath,subj,filesep]);
LLRin=[]; choices=[];
for s = 1:length(fsess)-2
    fblock = dir([loadpath,subj,filesep,fsess(s+2).name,filesep,'Behaviour',filesep,'*.mat']);
    fsmp = dir([loadpath,subj,filesep,fsess(s+2).name,filesep,'Sample_seqs',filesep,'*.mat']);
    for b = 1:length(fblock)
        load([fblock(b).folder,filesep,fblock(b).name])
        load([fsmp(b).folder,filesep,fsmp(b).name])
        
        % Converting sample and choice values to appropriate signs for choice regressions
        stimIn = round(stimIn.*-1);
        Cchoices = Behav(:,2)-1;
        
        % Convert stimulus values to LLRs
        if size(stimIn,1)>length(Cchoices), stimIn = stimIn(1:length(Cchoices),:); end  % trimming stimIn trials in case they exceed .mat trials (happens if block was terminated prematurely)
        CLLRin = log(normpdf(stimIn,gen.mu(2)*-1,gen.sigma(2))./normpdf(stimIn,gen.mu(1)*-1,gen.sigma(1)));
        
        % Concatenate
        LLRin = [LLRin; CLLRin(Cchoices==0 | Cchoices==1,:)];
        choices = [choices; Cchoices(Cchoices==0 | Cchoices==1)];
    end
end

for t = 1:size(LLRin,1)  % Calculating number of samples per trial, to be passed on to objective function calculation for fast indexing
    nsamps(t,1) = length(find(~isnan(LLRin(t,:))));
end

% Defining PSO options (see pso_Trelea_vectorized.m for details)
P(1)=0;    P(2)=1500;     P(3)=200;    P(4:13)=[1.6 1.9 0.9 0.4 400 1e-25 250 NaN 0 1];
% display  n_iterations  n_particles       acceleration, inertia, tolerance, etc

% Seeding first n particles with parameters drawn from realistic distributions
n_seeded = 75;
PSOseedValue=[];

PSOseedValue(1:n_seeded,end+1) = seeds.H(1)+(randn(n_seeded,1).*seeds.H(2));
if ~isempty(find(PSOseedValue(:,end)<range_p.H(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.H(1)),end) = range_p.H(1); end
if ~isempty(find(PSOseedValue(:,end)>range_p.H(2))), PSOseedValue(find(PSOseedValue(:,end)>range_p.H(2)),end) = range_p.H(2); end

PSOseedValue(1:n_seeded,end+1) = seeds.noise(1)+(randn(n_seeded,1).*seeds.noise(2));
if ~isempty(find(PSOseedValue(:,end)<range_p.noise(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.noise(1)),end) = range_p.noise(1); end
if ~isempty(find(PSOseedValue(:,end)>range_p.noise(2))), PSOseedValue(find(PSOseedValue(:,end)>range_p.noise(2)),end) = range_p.noise(2); end

PSOseedValue(1:n_seeded,end+1) = seeds.B(1)+(randn(n_seeded,1).*seeds.B(2));
if ~isempty(find(PSOseedValue(:,end)<range_p.B(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.B(1)),end) = range_p.B(1); end
if ~isempty(find(PSOseedValue(:,end)>range_p.B(2))), PSOseedValue(find(PSOseedValue(:,end)>range_p.B(2)),end) = range_p.B(2); end

PSOseedValue(1:n_seeded,end+1) = seeds.incon(1)+(randn(n_seeded,1).*seeds.incon(2));
if ~isempty(find(PSOseedValue(:,end)<range_p.incon(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.incon(1)),end) = range_p.incon(1); end
if ~isempty(find(PSOseedValue(:,end)>range_p.incon(2))), PSOseedValue(find(PSOseedValue(:,end)>range_p.incon(2)),end) = range_p.incon(2); end

% Concatenating parameter ranges
par_range = [range_p.H; range_p.noise; range_p.B; range_p.incon];

% Randomly seeding remaining particles within prespecified bounds
PSOseedValue(size(PSOseedValue,1)+1:P(3),1:size(PSOseedValue,2)) = normmat(rand([P(3)-n_seeded,size(PSOseedValue,2)]),par_range',1);

% Running PSO routine
[output,te,tr] = pso_Trelea_vectorized_Glaze('CE_H_noise_B_IU',length(mv),mv,par_range,0,P,'goplotpso',PSOseedValue);

% Store parameter estimates
pm_fit = output(1:end-1)
err = output(end)

% Save generative stats and fitted parameters for this iteration
save([savepath,subj,'_fit.mat'],'pm_fit','err','gen','P','PSOseedValue','seeds','range_p','te','tr')

clear global
