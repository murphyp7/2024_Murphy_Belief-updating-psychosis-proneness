% Uses Glaze model to generate finite number of observed choices for a
% given combination of H (hazard rate), noise and IU parameters, and
% attempts to recover parameters using particle swarm optimization

% ntrials = number of trials to simulate
% i = iteration number
% pm(1:3) = [H, noise, IU]

function [pm_fit,err] = Glaze_parameter_recoveryIU(iter)

nIters = 150;  % number of iterations for parameter recovery

% make vectors of parameter combinations and pick one for current job
ntrials = [1376, 10000];  % full range: [1250, 2150, 10000]
Hs = [0.03 0.06 0.1];  % full range: [0.08 0.1 0.2]
noises = [0.001 2];  % full range: [0.001 1 3]
IUs = [1 1.4];

pm_in=[];
for t = 1:length(ntrials)
    for h = 1:length(Hs)
        for n = 1:length(noises)
            for i = 1:length(IUs)
                pm_in(end+1,:) = [ntrials(t) Hs(h) noises(n) IUs(i)];
            end
        end
    end
end

nts = pm_in(iter,1);
H = pm_in(iter,2);
noise = pm_in(iter,3);
IU = pm_in(iter,4);

% Set global variables to pass to PSO algorithm
global LLRin choices nsamps

% Path stuff
addpath /mnt/homes/home024/pmurphy/Surprise_scz/param_recovery
addpath(genpath('/mnt/homes/home024/pmurphy/Surprise_accumulation/Simulations/particle_swarm_Glaze'))
savepath = '/mnt/homes/home024/pmurphy/Surprise_scz/param_recovery/sims/';

% Define task generative distributions
genH = 0.1;
muE = [-17 17];     % means
sigmaE = [27 27];   % standard deviations
rangeE = [-90 90];  % range for truncation
maxsamps = 10;        % number evidence samples per trial
pshort = 0.25;      % proportion of trunctated sequences

% Define parameter estimation settings
range_p.H = [0.0001 0.5];  % parameter bounds
range_p.noise = [0.0001 10];
range_p.incon = [0.01 10];

mv = [0.005;...    % maximum particle velocities
    0.2;...
    0.04]';

seeds.H = [0.1 0.04];  % good seed distributions for a selection of particles - [mean sd]
seeds.noise = [1.5 0.5];
seeds.incon = [1 0.5];

% Seed random number generator
seed = round(sum(100*clock)); %never the same seed
rand('state', seed);

% Loop through each iteration
for i = 1:nIters
    
    % Generate stimulus sequences and convert to LLRs
    LLRin = Gen_LLR(muE,sigmaE,rangeE,genH,pshort,maxsamps,nts);
    
    for t = 1:size(LLRin,1)  % Calculating number of samples per trial, to be passed on to objective function calculation for fast indexing
        nsamps(t,1) = length(find(~isnan(LLRin(t,:))));
    end
    
    % Generate choice probabilities give ground-truth parameters
    CPs = Glaze_sim_fast_IU(LLRin,nsamps,H,1,noise,IU);
    
    % Generate noisy binary choices given choice probabilities
    choices = zeros(nts,1);
    rnds = rand(nts,1);
    choices(CPs>=0.5 & rnds<=CPs) = 1;
    choices(CPs<0.5 & rnds<CPs) = 1;
    
    % Defining PSO options (see pso_Trelea_vectorized.m for details)
    P(1)=0;  P(2)=1500;     P(3)=200;    P(4:13)=[1.6 1.9 0.9 0.4 400 1e-25 250 NaN 0 1];
    % display  n_iterations  n_particles       acceleration, inertia, tolerance, etc
    
    % Seeding first n particles with parameters drawn from realistic distributions
    n_seeded = 75;  % default: 75
    PSOseedValue=[];
    
    PSOseedValue(1:n_seeded,end+1) = seeds.H(1)+(randn(n_seeded,1).*seeds.H(2));
    if ~isempty(find(PSOseedValue(:,end)<range_p.H(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.H(1)),end) = range_p.H(1); end
    if ~isempty(find(PSOseedValue(:,end)>range_p.H(2))), PSOseedValue(find(PSOseedValue(:,end)>range_p.H(2)),end) = range_p.H(2); end
    
    PSOseedValue(1:n_seeded,end+1) = seeds.noise(1)+(randn(n_seeded,1).*seeds.noise(2));
    if ~isempty(find(PSOseedValue(:,end)<range_p.noise(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.noise(1)),end) = range_p.noise(1); end
    if ~isempty(find(PSOseedValue(:,end)>range_p.noise(2))), PSOseedValue(find(PSOseedValue(:,end)>range_p.noise(2)),end) = range_p.noise(2); end
    
    PSOseedValue(1:n_seeded,end+1) = seeds.incon(1)+(randn(n_seeded,1).*seeds.incon(2));
    if ~isempty(find(PSOseedValue(:,end)<range_p.incon(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.incon(1)),end) = range_p.incon(1); end
    if ~isempty(find(PSOseedValue(:,end)>range_p.incon(2))), PSOseedValue(find(PSOseedValue(:,end)>range_p.incon(2)),end) = range_p.incon(2); end

    % Concatenating parameter ranges
    par_range = [range_p.H; range_p.noise; range_p.incon];
    
    % Randomly seeding remaining particles within prespecified bounds
    PSOseedValue(size(PSOseedValue,1)+1:P(3),1:size(PSOseedValue,2)) = normmat(rand([P(3)-n_seeded,size(PSOseedValue,2)]),par_range',1);
    
    % Running PSO routine
    [output,~,~] = pso_Trelea_vectorized_Glaze('CE_H_noise_IU',length(mv),mv,par_range,0,P,'goplotpso',PSOseedValue);
    
    % Store parameter estimates
    pm_fit(i,1:length(output)-1) = output(1:end-1);
    err(i,1) = output(end);
end

% Save generative stats and fitted parameters for this iteration
save([savepath,'Param_rec_H',num2str(H),'_noise',num2str(noise),'_IU',num2str(IU),'_t',num2str(nts),'.mat'],'pm_fit','err','muE','sigmaE','rangeE','nsamps','P')


