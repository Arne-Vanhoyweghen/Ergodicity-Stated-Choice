function computeBayesian(dataSource,nBurnin,nSamples,nThin,nChains,subjList,doParallel,startDir,nTrials,seedChoice)

% computeBayesian read in and manipulate all data,
% as well as specify variables needed to run the Bayesian model through JAGS.

% dataSource    - which data source to estimate data on; non-timed (1) or timed (2)
% nBurnin       - specifies how many burn in samples to run
% nSamples      - specifies the number of samples to run
% nThin         - specifies the thinnning number
% nChains       - specifies number of chains
% subjList      - list of subject numbers to include
% doParallel    - sets whether to run chains sequentially or in parallel
% startDir      - root directory for the repo
% nTrials       - number of trials in experiment

%% Set paths
cd(startDir);%move to starting directory
matjagsdir=fullfile(startDir,'/Bayesian_utils/matjags');
addpath(matjagsdir)
jagsDir=fullfile(startDir,'/Bayesian_utils/JAGS');
addpath(jagsDir)
dataDir=fullfile(startDir,'/data');

%% Choose & load data
switch dataSource
    case{1}, load(fullfile(dataDir, 'timed_data.mat')); source = 'timed_data';
    case{2}, load(fullfile(dataDir, 'non_timed_data.mat')); source  = 'non_timed_data';
end %dataSource

%% Choose JAGS file
modelName = 'JAGS_isoelastic_parameter_estimation';

%% Set key variables
nConditions=2;%number of dynamics
doDIC=0;%compute Deviance information criteria? This is the hierarchical equivalent of an AIC, the lower the better
nSubjects=length(subjList);%number of subjects

%% Set bounds of hyperpriors
%hard code the upper and lower bounds of hyperpriors, typically uniformly
%distributed. These values will be imported to JAGS.

% beta - prior on log since cannot be less than 0;
muLogBetaL=-2.3;muLogBetaU=3.4; %bounds on mean of the distribution of log beta
sigmaLogBetaL=0.01;sigmaLogBetaU=1.6;%bounds on the std of the distribution of log beta

% eta
muEtaL=-5; muEtaU=5; %bounds on the mean of the distribution of eta
sigmaEtaL=0.01; sigmaEtaU=1.6; %bounds on the std of the distribution of eta

%% Print information for user
disp('**************');
disp(['Mode: ', modelName])
disp(['dataSource: ', source])
disp(['started: ',datestr(clock)])
disp('**************');

%% Initialise matrices
%initialise matrices with nan values of size subjects x conditions x trials
dim = nan(nSubjects,nConditions,nTrials); %specify the dimension
choice = dim; %initialise choice data matrix
dwLU=dim; dwLL=dim; dwRU=dim; dwRL=dim;%initialise wealth increments
w=dim;%initialise wealth


%% Compile choice & gamble data
trialInds = 1:nTrials;
for c = 1:nConditions
    switch c
        case {1} %eta = 0
            choice(:,c,trialInds)=choice_add(:,trialInds);
            dwLU(:,c,trialInds)=x1_1_add(:,trialInds);
            dwLL(:,c,trialInds)=x1_2_add(:,trialInds);
            dwRU(:,c,trialInds)=x2_1_add(:,trialInds);
            dwRL(:,c,trialInds)=x2_2_add(:,trialInds);
            w(:,c,trialInds)=wealth_add(:,trialInds);

        case {2}% eta=1
            choice(:,c,trialInds)=choice_mul(:,trialInds);
            dwLU(:,c,trialInds)=x1_1_mul(:,trialInds);
            dwLL(:,c,trialInds)=x1_2_mul(:,trialInds);
            dwRU(:,c,trialInds)=x2_1_mul(:,trialInds);
            dwRL(:,c,trialInds)=x2_2_mul(:,trialInds);
            w(:,c,trialInds)=wealth_mul(:,trialInds);
    end %switch
end %c

%% Nan check
disp([num2str(length(find(isnan(choice)))),'_nans in choice data']);%nans in choice data is okay
disp([num2str(length(find(isnan(dwLU)))),'_nans in gambles Left Upper matrix'])% nans in gamble matrices is not
disp([num2str(length(find(isnan(dwLL)))),'_nans in gambles Left Lower matrix'])
disp([num2str(length(find(isnan(dwRU)))),'_nans in gambles Right Upper matrix'])
disp([num2str(length(find(isnan(dwRL)))),'_nans in gambles Right Lower matrix'])
disp([num2str(length(find(isnan(w)))),'_nans in wealth matrix'])

%% Configure data structure for graphical model & parameters to monitor
%everything you want jags to use
dataStruct = struct(...
            'nSubjects', nSubjects,'nConditions',nConditions,'nTrials',nTrials,...
            'w',w,'dwLU',dwLU,'dwLL',dwLL,'dwRU',dwRU,'dwRL',dwRL,'y',choice,...
            'muLogBetaL',muLogBetaL,'muLogBetaU',muLogBetaU,'sigmaLogBetaL',sigmaLogBetaL,'sigmaLogBetaU',sigmaLogBetaU,...
            'muEtaL',muEtaL,'muEtaU',muEtaU,'sigmaEtaL',sigmaEtaL,'sigmaEtaU',sigmaEtaU);
%all parameters you want JAGS to output
for i = 1:nChains
    monitorParameters = {'mu_eta','tau_eta','sigma_eta',...
                            'mu_log_beta','tau_log_beta','sigma_log_beta',...
                            'beta_i', 'beta_g','eta_i', 'eta_g'};
    S=struct; init0(i)=S; %sets initial values as empty so randomly seeded
end %i

%% Run JAGS sampling via matJAGS
tic;fprintf( 'Running JAGS ...\n' ); % start clock to time % display

[samples, stats] = matjags( ...
    dataStruct, ...                           % Observed data
    fullfile(jagsDir, [modelName '.txt']), ...% File that contains model definition
    init0, ...                                % Initial values for latent variables
    whichJAGS,...                             % Specifies which copy of JAGS to run on
    'doparallel' , doParallel, ...            % Parallelization flag
    'nchains', nChains,...                    % Number of MCMC chains
    'nburnin', nBurnin,...                    % Number of burnin steps
    'nsamples', nSamples, ...                 % Number of samples to extract
    'thin', nThin, ...                        % Thinning parameter
    'dic', doDIC, ...                         % Do the DIC?
    'monitorparams', monitorParameters, ...   % List of latent variables to monitor
    'savejagsoutput' , 1 , ...                % Save command line output produced by JAGS?
    'verbosity' , 1 , ...                     % 0=do not produce any output; 1=minimal text output; 2=maximum text output
    'cleanup' , 0 ,...                        % clean up of temporary files?
    'rndseed',seedChoice);                    % Randomise seed; 1=no; 2=yes

toc % end clock

%% Save stats and samples
disp('saving samples...')
save(fullfile(dataDir, append('Bayesia_parameter_estimation','_',source)),'stats','samples','-v7.3')

%% Print readouts
disp('stats:'),disp(stats)%print out structure of stats output
disp('samples:'),disp(samples);%print out structure of samples output
try
    rhats=fields(stats.Rhat);
    for lp = 1: length(rhats)
        disp(['stats.Rhat.',rhats{lp}]);
        eval(strcat('stats.Rhat.',rhats{lp}))
    end
catch
    disp('no field for stats.Rhat')
end
