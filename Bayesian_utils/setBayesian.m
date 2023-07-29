function setBayesian(dataSource,whichQuals,doParallel,startDir,seedChoice)

% setHLM sets up the Bayesian model

% This function takes the following inputs:

% dataSource    - which data source to estimate data on; non-timed (1) or timed (2)
% whichQuals    - sets the quality of run; min 1 and max 5
% doParallel    - whether to run chains sequentially (0) or in parallel (1)
% seedChoice    - set whether to do manual seed choice (1), or random seed (2)

%% Specifies qualities to be selected from
nBurnin      = [1e2,1e3,1e4,2e4,4e4];  %from 100 to 40k
nSamples     = [5e1,5e2,5e3,1e4,2e4];  %from 50 to 20k
nChains      = [4,4,4,4,4];            %Keep this to 4
nThin        = 10;                     %thinnning factor, 1 = no thinning, 2=every 2nd etc.

%% Specifies subjects, trials and directory_name
switch dataSource
    case {1}, subjList = 1:40;
    case {2}, subjList = 1:41;
end

nTrials = 40; %each participant did 40 trials in each condition (additive and multiplicative)

computeBayesian(dataSource,nBurnin(whichQuals),nSamples(whichQuals),nThin,nChains(whichQuals),subjList,doParallel,startDir,nTrials,seedChoice)
