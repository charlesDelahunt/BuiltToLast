
function [expParams]  = twoOdorsInjuryExperimentParamsFeb_fn( trClasses, classLabels, numVal,...
                                                              numPuffsForTrain1AndForTrain2)
%
% This function defines parameters of a time-evolution experiment: overall timing, stim timing and  
% strength, octo timing and strength, lowpass window parameter,  etc.
%
% This experiment has two odors, one trained and one control. Two flavors: 
%                   (i) odor + ambient noise (training odor) vs pure ambient noise (control odor). 
%                   (ii) odor 1 (training odor) vs odor 2 (control odor).
% These odors are created by 'twoOdorsCreatePuffs_fn' or 'odorAndAddedNoiseCreatePuffs_fn', and are 
% contained in 'featureArray' in the main run script.
% The goal is to see if (a) injury makes it harder to discriminate an odor from ambient
%         noise; and (b) if subsequent training makes it easier.
%
% Inputs: 
%       1. trClasses = 1 x numTrainingPuffs vector: class labels of trained class (de facto == 1)
%       2. classLabels = [1,2]
%       3. numVal = number of puffs in each validation batch
%       4. numPuffsForTrain1AndForTrain2 = 1 x 2 vector: how many puffs to train on 
%          ie [numTrain1, numTrain2]
% Output:
%   1. expParams: struct with experiment info.

% Copyright (c) 2018-2020 Charles B. Delahunt.  delahunt@uw.edu
% MIT License

%----------------------------------------------------------------

% Order of time periods:
%       1. no event period: allow system to settle to a steady state spontaneous FR baseline
%       2. baseline period: deliver a group of odor+noise, then a group of  noise only
%       3. no event buffer 
%       5. training period:  deliver odor+noise + octopamine + allow hebbian updates 
%       6. no event buffer 
%       7. post-training period: deliver a group of odor+noise, then a group of  noise only
%       8. injury: apply injury
%       9. post-injury period: deliver a group of odor+noise, then a group of  noise only
%      10: 2nd training: deliver odor+noise + octopamine + allow hebbian updates 
%      11: post-2nd-training period: deliver a group of odor+noise, then a group of  noise only 

stimMag =  3; %   2019: = 3  But this does not work with 2017 hebbian params. 2017: = 7.  5; 
% (stimMag = stim magnitudes as passed into AL. See original version in smartAsABug codebase)
stimLength = 0.22;
classLabels = [1, 2];   % 1 = odor, 2 = ambient noise

% abbreviations:
nC = length(classLabels);  
nV = numVal;
nT = numPuffsForTrain1AndForTrain2; % 2 x 1 vector

%% Define the time span and events:
% Generate a vector of timepoints to indicate when puffs occur.

step = 3;  % the time between digits (3 seconds)
trStep = step + 2;  % allow more time between training digits  

expParams.simStart = -30;  % use negative start-time for convenience (artifact)

% Baseline period:
% do a loop, to allow gaps between class groups:
baselineValTimes = [];  
startTime = 30;
gap = 10;
for i = 1:nC
    baselineValTimes = [ baselineValTimes, startTime : step : startTime + (nV - 1)*step  ];  
    startTime = max(baselineValTimes) + gap;
end 
endOfBaseline = max(baselineValTimes) + 25;  % include extra buffer before training

% Training period:
if nT(1) > 0
    train1Times = endOfBaseline : trStep : endOfBaseline + ( nT(1) - 1)*trStep;   
    endOfTrain1 = max(train1Times) + 1; 
    startOfTrain1 = min(train1Times) - 1;
else
    train1Times = [];
    endOfTrain1 = endOfBaseline + 2;
    startOfTrain1 = endOfBaseline + 1;
end

% Val period:
% do a loop, to allow gaps between class groups: 
postTrain1ValTimes = [];
startTime = endOfTrain1 + 25;   % include buffer before Validation
for i = 1:nC
    postTrain1ValTimes = [ postTrain1ValTimes, startTime : step : startTime + (nV - 1)*step  ];  
    startTime = max(postTrain1ValTimes) + gap;
end
endOfVal = max(postTrain1ValTimes) + 4;  

% Injury applied
injuryTime = endOfVal;

% post-injury validation period:
postInjuryValTimes = [];
startTime = injuryTime + 25;  % to allow polling of spontaneous activity
for i = 1:nC
    postInjuryValTimes = [ postInjuryValTimes, startTime : step : startTime + (nV - 1)*step  ]; 
    startTime = max(postInjuryValTimes) + gap;
end
endOfPostInjuryVal = max(postInjuryValTimes) + 4;  

% 2nd train period:
% Training period:
startTime = endOfPostInjuryVal + 25;
if nT(2) > 0
    train2Times = startTime : trStep : startTime + ( nT(2)  - 1)*trStep; % vector of timepoints, 
    % one digit every 'trStep' seconds
    endOfTrain2 = max(train2Times) + 1;
    startOfTrain2 = min(train2Times) - 1;
else
    train2Times = [];
    endOfTrain2 = startTime + 2;
    startOfTrain2 = startTime + 1;
end
 
% post-train2 validation period:
postTrain2ValTimes = [];
startTime = endOfTrain2 + 25; % include buffer  
for i = 1:nC
    postTrain2ValTimes = [ postTrain2ValTimes,  startTime : step : startTime + (nV - 1)*step  ];  
    startTime = max(postTrain2ValTimes) + gap;
end
endOfPostTrain2Val = max(postTrain2ValTimes) + 4;  

%% assemble vectors of stimulus data for export:

% Assign the classes of each stim. Assign the baseline and val in blocks, and the training stims in 
% the order passed in:
whichClass = zeros( size ( [baselineValTimes, train1Times, postTrain1ValTimes,...
                            postInjuryValTimes, train2Times, postTrain2ValTimes ] ) );
for c = 1:nC    
    whichClass( (c-1)*nV+ 1: c*nV )  = classLabels(c);  % baseline (pre-train1)
    whichClass( nV*nC + nT(1) + (c-1)*nV + 1 : nV*nC + nT(1) + c*nV ) = ...
                                                                      classLabels(c);  % post-train1
    whichClass( nV*nC + nT(1) + nV*nC + (c-1)*nV + 1 :  nV*nC + nT(1) + nV*nC + c*nV ) = ...
                                                                      classLabels(c);  % post-injury
    whichClass( nV*nC + nT(1) + nV*nC + nV*nC + nT(2) + (c-1)*nV + 1 : ...
                nV*nC + nT(1) + nV*nC + nV*nC + nT(2) + c*nV  ) = classLabels(c);  % post-train2
end
% training puffs are of classLabel 1 only:
whichClass( nV*nC + 1 : nV*nC + nT(1) ) = classLabels(1); % train1 follows one set of val puffs
whichClass( nV*nC*3 + nT(1) + 1 : nV*nC*3  + nT(1) + nT(2) ) = classLabels(1); % train2 happens 
% after 3 sets of val puffs and the train1 puffs 

expParams.whichClass = whichClass;

stimStarts =  [baselineValTimes, train1Times, postTrain1ValTimes, postInjuryValTimes,...
               train2Times, postTrain2ValTimes ] ; 
expParams.stimStarts = stimStarts; % starting times
expParams.durations = stimLength*ones( size( stimStarts ) );  % durations
expParams.classMags = stimMag*ones( size( stimStarts ) );  % magnitudes 

% octopamine input timing:
expParams.octoMag = 1;          % 2017: = 0.6
expParams.octoStart = [ train1Times, train2Times ];   
expParams.durationOcto = 1;

% heb timing: Hebbian updates are enabled 25% of the way into the stimulus, and
% last until 75% of the way through (ie active during the peak response period)
expParams.hebDelay = 1.25*stimLength;
expParams.hebStarts = [train1Times, train2Times] + expParams.hebDelay;
expParams.lengthOfHeb = 2.5*stimLength;
expParams.hebDurations = (expParams.lengthOfHeb)*ones(size( [train1Times, train2Times] ));
% expParams.startTrain = min(expParams.hebStarts);
% expParams.endTrain = max(expParams.hebStarts) + max(expParams.TrainDurations);

%% Other time parameters required for time evolution book-keeping:

% The numbers 1,2,3 do refer to time periods where spont responses are allowed to settle before 
% recalibration.
expParams.startPreNoiseSpontMean1 = -25;     
expParams.stopPreNoiseSpontMean1 = -15;
% Currently no change is made in start/stopSpontMean2. So spontaneous behavior may be stable in 
% this range.
expParams.startSpontMean2 = -10;
expParams.stopSpontMean2 = -5;
% currently, spontaneous behavior is steady-state by startSpontMean3.
expParams.startSpontMean3 = 0;
expParams.stopSpontMean3 = 10;

% % these maybe not necessary
% expParams.preTrain1PollTime = min(train1Times) - 2;
% expParams.postTrain1PollTime = max(train1Times) + 2;
% expParams.preInjuryPollTime = injuryTime - 2;
% expParams.postInjuryPollTime = injuryTime + 2;
% expParams.preTrain2PollTime = min(train2Times) - 2;
% expParams.postTrain2PollTime = max(train2Times) + 2;

% timePoints for plotting EN responses:
% spontaneous response periods, pre and post training, to view effect of training on spontaneous FRs:
expParams.preTrain1SpontStart = expParams.stopSpontMean3 + 2;
expParams.preTrain1SpontStop = startOfTrain1 - 3;
expParams.postTrain1SpontStart = endOfTrain1+ 5;
expParams.postTrain1SpontStop = min(postTrain1ValTimes) - 3;
expParams.postInjurySpontStart = injuryTime + 5;
expParams.postInjurySpontStop  = min(postInjuryValTimes) - 3;
expParams.postTrain2SpontStart = endOfTrain2+ 5;
expParams.postTrain2SpontStop = min(postTrain2ValTimes) - 3;

expParams.startTrain1 = startOfTrain1;
expParams.endTrain1 = endOfTrain1;
expParams.startTrain2 = startOfTrain2;
expParams.endTrain2 = endOfTrain2;
 
% hamming filter window parameter (= width of transition zone in seconds). The lp filter is applied 
% to odors and to octo
expParams.lpParam =  0.12;

expParams.simStop = max(stimStarts) + 10;

expParams.injuryTime = injuryTime;


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
