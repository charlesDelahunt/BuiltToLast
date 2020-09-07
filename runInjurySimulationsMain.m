%
% Main script to test effects of injury and training on a moth brain.
% Reproduces experiments described in 'Built to Last', by C.B. Delahunt, % P.D. Maia, and J.N. Kutz
% Outputs: EN timecourse (recommended); one matfile per run/moth instance with experiment results; 
%          statistical plots (if wished).
%
% To run:
%   1. Edit USER ENTRIES
%   2. Run
%
% Order of events:
%   Within the loop over number of simulations, for each simulation:
%       1. Create a stream of odor stimuli according to experiment needs.
%       2. Create injury vectors.
%       3. Generate a new moth instance according to a pre-defined template in 2 steps: 
%           a) call 'specifyModelParamsMnist_fn' 
%           b) create connection matrices via 'initializeConnectionMatrices_fn'
%       4. Load the experiment parameters.
%       5. Run the simulation with 'sdeInjuryWrapperFeb_fn'
%       6. Plot results, save results to matfile.
%
% Dependencies: Matlab, Statistics and machine learning toolbox, Signal processing toolbox

% Copyright (c) 2018 - 2020, Charles B. Delahunt.  delahunt@uw.edu
% MIT License

%-------------------------------------------------

close all
clear  

%% USER ENTRIES: 
 
% Flags to show various images:   
showThumbnailsUsed = 1; % 1 -> Show projections of experiment odors onto glomeruli. 0 -> ignore.
showENPlots = [0, 1];  %  arg1 == 1 -> show statistical plots of EN response 
% changes. arg2 == 1 -> plot the EN timecourse. argin* == 0 -> ignore.
%---------

% Experiment parameters:

% Experiment type, choose one of:
% 'oneOdor'; 'trainedOdorVsControlOdor'; 'trainedOdorPlusNoiseVsControlNoise';
experimentType = 'oneOdor'; 

inputNoiseLevelList = 0.4; % [0.2, 0.4, 0.6, 0.8, 1.0];  % relative magnitude of noise, used in 
                         % experiment with {odor+noise, noise only}
                         
% Choose the number of training puffs. For 'oneOdor' experiment, numPuffsInTrainRound1 == 0.
numPuffsInTrainRound1 = 2;  
numPuffsInTrainRound2 = 8; 

% Choose injury type:  'Fas' or 'Ablate'
injuryType = 'Fas'; 

% Define injury region:  'RNs' or 'PNs'
injuryRegion = 'RNs'; % 
  
% Injury level:
injuryFrList = 0.4; % sweep: [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6] 

% Other:
normalizeOdorClassMags = 1;  % 0 or 1. For 2 odors:  if 1, the two classes will have roughly equal 
                   % summed input mags into the AL. For odor + noise vs noise this is a dummy.

odorSpread = 0.3;  % Fraction of gloms that an odor hits (mean).
odorStrengthStd = 0.3;  % std dev of odor strength to each glom (as fraction of mean strength)
inputNoiseSpread = 0.5;  % fraction of gloms that input noise in {odor + noise, noise} will affect 
                         % (the gloms change with each puff)
%---------

% Moth characteristics-under-test:
numPIperGList =  1;  % for RN injury: 1.  PN injury sweep: [0 1 2 3 5]; 

alNoiseList =  1; %   [0.01,  0.33, 0.67, 1.0, 1.33];       
%---------
  
numRuns = 1;  % how many runs you wish to do with this moth template.
%---------

% To save results if wished:
saveAllNeuralTimecourses = 0; % 0 -> save only EN (ie readout) timecourses.  Caution: 1 -> very high
                              % memory demands, hinders longer runs.
 
saveResultsImageFolder = []; % String. If non-empty, images will be saved here, if showENPlots is 
                             % also non-zero).
experimentPrefix = 'prefix';  % will get key params appended to it.
saveResultsDataFolder = 'experimentResultsFolder';

% END OF USER ENTRIES
%----------------------------------------------------------------------------------------------

%% Misc book-keeping 

nG = 60;

classLabels = 1:2;  % Two experiments: (i) 1 = odor + noise, 2 = noise only;
                    % (ii) 1 = odor1, 2 = odor2.
nC = length(classLabels);  % = number of distinct odors.

numVal = 11;  % number of puffs used in standard batches (eg baseline sets), per odor 

% make a vector of class one, the only one trained on :
trClasses = ones( classLabels(1), numPuffsInTrainRound1 + numPuffsInTrainRound2);  

% Experiment set-up details, applies to both {odor1, odor2} and {odor + noise, noise only}
experimentFn = @twoOdorsInjuryExperimentParamsFeb_fn; 
% % For one odor, no pre-training. Injury then training:
% experimentFn = @oneOdorInjuryPlusTrainExperimentNoMiscOcto_31may;
 
octoMag2017 = 0.6; % relevant for one-odor experiment (ie if using 2017 params).
stimMag2017 = 7; 

% Because two-odor experiments require pre-training (to create discrimination), 
% new (c 2019) parameters for generating the moth instance differ slightly from old (c 2017) 
% parameters used in one odor experiments. The two key  differences are: 2019 has different hebTau 
% values, and includes weight decay.  
use2017Params =  strcmpi(experimentType, 'oneOdor');

injuryTypeDistribution = [];
if strcmpi(injuryType,'Fas')
    injuryTypeDistribution = [0.35, 0.35, 0.15];  % other 0.25 remain healthy, are not listed.
end
if strcmpi(injuryType,'Ablate')
    injuryTypeDistribution = [1, 0, 0];   % Ablation
end
if isempty(injuryTypeDistribution)
    disp('Error: set injuryType = FAS or Ablate')
end
    
odorStreamCreateFn = [];
if strcmpi(experimentType, 'oneOdor')
    numPuffsInTrainRound1 = 0;
    inputNoiseLevelList = 0;
    odorStreamCreateFn = @odorAndAddedNoiseCreatePuffs_fn;
    modelParamsFn = @specifyModelParamsSymmetricPnsPins_2017values_Fn;
    experimentTag = 'oneOdor'; 
end
if strcmpi(experimentType, 'trainedOdorVsControlOdor')
    inputNoiseLevelList = 0;
    odorStreamCreateFn = @twoOdorsCreatePuffs_fn;
    modelParamsFn = @specifyModelParamsSymmetricPnsPins_Fn;
    experimentTag = 'twoOdors';
end
if strcmpi(experimentType, 'trainedOdorPlusNoiseVsControlNoise')
	odorStreamCreateFn = @odorAndAddedNoiseCreatePuffs_fn;
    modelParamsFn = @specifyModelParamsSymmetricPnsPins_Fn;
    experimentTag = 'odorVsNoise';
end
if isempty(odorStreamCreateFn)
    disp(['Error: Set experimentType = oneOdor, trainedOdorVsControlOdor, or ',...
          ' trainedOdorPlusNoiseVsControlNoise'])
end
%----------------------------------
 
% Loop through the number of simulations specified:
counter = 0; 
for alNoise = alNoiseList
for numPIperG =  numPIperGList
for injuryFr = injuryFrList
for inputNoiseLevel = inputNoiseLevelList
for run =  1:numRuns 
    
    counter = counter + 1;

    if true  % condition to do run
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
          
        resultsFilename =  [experimentPrefix, '_', experimentTag, '_', injuryType, 'To',...
            injuryRegion, '_qn' ,num2str(numPIperG),...
            '_alNoi', strrep(num2str(alNoise),'.','p'), '_inputNoi',...
            strrep(num2str( inputNoiseLevel) ,'.','p' ), '_injFr',...
            strrep(num2str( injuryFr) ,'.','p' ), '_run', num2str(run),...
            '_', num2str(counter) ];
        disp(['Running ' resultsFilename])
        
        %% Create a stream of puffs for use in experiment.
        % For two odor experiments, order of puffs is: standard (baseline), trainRound1, standard  
        % (post-train1), [injury occurs], standard (post-injury), trainRound2, standard (post-train2).  
        % Only odor1 is trained. There is only one EN.
        % For one odor experiments, trainRound1 has 0 puffs, and odor2 puffs have 0 magnitude.

        numPuffs = 4*numVal + numPuffsInTrainRound1 + numPuffsInTrainRound2;

        odorStreams = odorStreamCreateFn(nG, numPuffs, odorSpread,...
            odorStrengthStd, inputNoiseLevel, inputNoiseSpread,...
            normalizeOdorClassMags);
        % odorStreams is nG (ie nR) x numPuffs x nC array.
        % OdorStreams is the projection of the odor puffs onto the RNs.

        % % show the final versions of thumbnails to be used, if wished:
        if showThumbnailsUsed
            titleString = 'Input odor(s) projection(s) onto RNs'; 
            normalize = 1;
            showFeatureArrayThumbnailsInjuryFeb_fn(odorStreams,...
                showThumbnailsUsed, normalize, titleString );  
                % argin2 = number of images per class to show.
        end

        % Define the experiment parameters, including book-keeping for time-stepped evolutions, eg
        % when octopamine occurs, windowing of Firing rates, etc:
        experimentParams = experimentFn( trClasses, classLabels, numVal, [numPuffsInTrainRound1,...
                                         numPuffsInTrainRound2] );
        
        % 2017 had different stimMag and octoMag:
        if use2017Params
            experimentParams.octoMag = octoMag2017;
            experimentParams.classMags = stimMag2017*ones(size(experimentParams.classMags));
        end

        %% Create injury vectors: 

        % define injuryParams for relevant neuron types: [ablateFr, reflectFr, filterFr]        
        % june 19:  24-hour pie chart is fed into fancyFilterParams.  

        if strcmpi(injuryRegion,'RNs')
            % RN pie injury: 
            injuryFn = @createInjuryFlags_1;     
            createInjuryParams.rn_pie = [ 0, 0, 1 ]; % for ease, do not change this. To model 
            % ablation, use injuryFr*[1,0,0] in fancyFilterParams.rn_pie below
            createInjuryParams.fancyFilterHandle = @rnCombinedDamageFn;
            fancyFilterParams.rn_pie = injuryFr*injuryTypeDistribution; % FAS injury types 
            fancyFilterParams.scalingFactor = 15;  % (used for LP filter component only). 
            % we need to know the expected maxFR from the antennae:
            temp = odorStreams(:);
            temp = temp(temp > 0);
            antennaeMax = (mean(experimentParams.classMags(:)))*  (mean(temp) + 2*std(temp));
            fancyFilterParams.antennaeMax = antennaeMax; % (used for LP filter component only) 
             % Max inputs from antennae are stimMag*(mean + 2*std of odorStreams).
             % The LP injury filter is a quadratic that behaves properly in the range [0,16]. 
             % So inputsFromAntennae / antennaeMax * scalingFactor has range [0, 1]*scalingFactor  
             % = [0, 15] as inputs to the quadratic during odor puffs.
            createInjuryParams.fancyFilterParams = fancyFilterParams;
        end
        if strcmpi(injuryRegion, 'PNs')
            % 2. PN pie injury:
            injuryFn = @createInjuryFlags_1;
            createInjuryParams.pn_pie = injuryFr*injuryTypeDistribution;
            createInjuryParams.pin_pie = injuryFr*injuryTypeDistribution;
            createInjuryParams.fancyFilterHandle = @pedroLowPassFilter_v2;
            fancyFilterParams.pnMax = 0.7; % (used for LP filter component only). empirical max = 
                                           % 0.7 to 0.8. % '1' was used in 2017 experiments
            fancyFilterParams.pinMax = 0.7; % pin = qn
            fancyFilterParams.scalingFactor = 10;  % (used for LP filter component only). 
            createInjuryParams.fancyFilterParams = fancyFilterParams;
        end

        %% Create a moth:

        % a) load template params with specify_params_fn:
        modelParams = modelParamsFn('dummy', 0); % argin2 = 0 -> don't save the generated 
                    % modelParams into a matfile. argin1 = name of file to save to.
        % modelParams is a struct
    
        % Modify post-injury neural noise in RNs if RNs are injured: Noise is higher post-injury
        % since there are fewer RNs contributing to the total, so lower effects from averaging:
        if strcmpi(injuryRegion, 'RNs')
            effectiveInjuryFr = min( 0.99, sum( [1, 0.5, 0.5].* fancyFilterParams.rn_pie ) );  
                      % ie get the effective total injury via 1*ablateFr + 0.5*reflectFr +
                      % 0.5*lowpassFr. This should come to about 0.6*injuryFr.
                      % the min( ) is to prevent zero case.
            modelParams.postInjuryNoiseRnRatio = 1/sqrt(1-effectiveInjuryFr);    
        else    % case: PN injury, so no modification
            modelParams.postInjuryNoiseRnRatio = 1; 
        end

        modelParams.numPIperG = numPIperG;
        modelParams.nPI = numPIperG*modelParams.nG;

        % c) Now populate the moth's connection matrices using the modelParams:
        modelParams = initializeConnectionMatricesSymmetricPnsPins(modelParams);  

        modelParams.trueClassLabels = classLabels;     % misc parameter tagging along
        modelParams.saveAllNeuralTimecourses = saveAllNeuralTimecourses;

        % get the injury params:
        injuryVectors = injuryFn( modelParams, createInjuryParams );

        %---------------------------------------------------------------------------------

        %% 4. run this experiment as sde time-step evolution:

        simResults = sdeInjuryWrapperFeb_fn( modelParams, experimentParams, odorStreams,... 
                                             injuryVectors ); 

        %-----------------------------------

        %% Experiment Results: EN behavior, classifier calculations: 

        if ~isempty(saveResultsImageFolder)
            if ~exist(saveResultsImageFolder)
                mkdir(saveResultsImageFolder)
            end
        end
        % Process the sim results to group EN responses by class and time:
        r = viewENresponsesInjuryFeb_fn(simResults, modelParams, experimentParams, showENPlots,...
                                        classLabels, resultsFilename, saveResultsImageFolder);

    
        % append the accuracy results, and other run data, to the first entry of r:
        r(1).modelParams = modelParams;  % will include all connection weights of this moth
        r(1).experimentParams = experimentParams; 
        r(1).K2Efinal = simResults.K2Efinal;

        if ~isempty(saveResultsDataFolder) 
            if ~exist(saveResultsDataFolder, 'dir' )
                mkdir(saveResultsDataFolder)
            end 
            % format of results filenames:  [experimentPrefix] _qn0_alNoi1_injFr0p6_run1_45.mat 

            save( fullfile(saveResultsDataFolder, resultsFilename ) , 'r')
        end
    end % if counter > N
end % for run  
end % for inputNoiseLevel
end % for injuryFr
end % for numPIperG
end % for ALnoise


% MIT license:
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
% associated documentation files (the "Software"), to deal in the Software without restriction,  
% including without limitation the rights to use, copy, modify, merge, publish, distribute,
% sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is 
% furnished to do so, subject to the following conditions: 
% The above copyright notice and this permission notice shall be included in all copies or 
% substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
% INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
% PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
% COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN 
% AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION 
% WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


