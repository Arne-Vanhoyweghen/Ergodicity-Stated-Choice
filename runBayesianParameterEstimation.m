% This script provides inputs to setBayesian (in the folder Bayesian_utils) for a given run

% It provides the following inputs when calling setHLM:

% dataSource    - which data source to estimate data on; non-timed (1) or timed (2)
% whichQuals    - sets the quality of run; min 1 and max 5
% doParallel    - whether to run chains sequentially (0) or in parallel (1)
% seedChoice    - set whether to do manual seed choice (1), or random seed (2)

% The idea is that this is written into by the user, then called by a
% cluster job via the terminal:

%%Specify startpath
restoredefaultpath
[startDir,~] = fileparts(mfilename('fullpath'));  %specify your starting directory here (where this script runs from)
addpath(fullfile(startDir,'/Bayesian_utils'));

%% Specify variables
dataSource = 1:2;
whichQuals = 1;
doParallel = 0;
seedChoice = 1;

%% Call setHLM
for source = dataSource
    setBayesian(source,whichQuals,doParallel,startDir,seedChoice)
end
