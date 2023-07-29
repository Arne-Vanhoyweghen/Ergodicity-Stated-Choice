% runHLM provides inputs to setHLM for a given run on the cluster

% It provides the following inputs when calling setHLM:

% dataSource    - which data source to estimate data on; non-timed (1) or timed (2)
% whichJAGS     - which copy of matjags to run on. this allows parallel jobs to run as long as they use different matjags
% whichQuals    - sets the order of qualities to run
% doParallel    - whether to run chains in parallel
% seedChoice    - set whether to do manual seed choice (1), or random seed (2)

% The idea is that this is written into by the user, then called by a
% cluster job via the terminal:

%%Specify startpath
restoredefaultpath
[startDir,~] = fileparts(mfilename('fullpath'));  %specify your starting directory here (where this script runs from)
addpath(fullfile(startDir,'/Bayesian_utils'));

%% Specify variables
dataSource = 1
whichQuals = 1;
doParallel = 0;
seedChoice = 1;

%% Call setHLM
setBayesian(dataSource,whichQuals,doParallel,startDir,seedChoice)
