function setBayesian(dataSource,whichJAGS,whichQuals,doParallel,startDir,seedChoice)

% setHLM sets up multiple HLM models to run sequentially according to inputs

% This function takes the following inputs:

% dataSource    - which data source to estimate data on; non-timed (1) or timed (2)
% whichJAGS     - which copy of matjags to run on. this allows parallel jobs to run as long as they use different matjags
% whichQuals    - sets the order of qualities to run
% doParallel    - whether to run chains in parallel
% seedChoice    - set whether to do manual seed choice (1), or random seed (2)

%% Specifies qualities to be selected from
numRuns      = 1;     %how many separate instances of an MCMC to run
nBurnin      = [1e2,1e3,1e4,2e4,4e4];  %from 100 to 40k
nSamples     = [5e1,5e2,5e3,1e4,2e4];  %from 50 to 20k
nChains      = [4,4,4,4,4];            %Keep this to 4
nThin        = 10;                     %thinnning factor, 1 = no thinning, 2=every 2nd etc.

%% Specifies subjects, trials and directory_name
switch dataSource
    case {1}, subjList = 1:41;
    case {2}, subjList = 1:40;,
end %dataSource
nTrials = 80;

%% Runs HLMs sequentially
computeBayesian(dataSource,nBurnin(whichQuals),nSamples(whichQuals),nThin,nChains(whichQuals),subjList,whichJAGS,doParallel,startDir,nTrials,seedChoice)
end
