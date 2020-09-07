
function [ results ] = viewENresponsesInjuryFeb_fn( simResults, modelParams, experimentParams,...
    showPlots, classLabels, resultsFilename, saveImageFolder)

% This version is designed to process an experiment that has two training
% sessions and one injury. So relevant periods for data collection are:
% preTrain1
% postTrain1 && preInjury
% postInjury && preTrain2
% postTrain2

% View readout neurons (EN):
% Color-code them dots by class and by concurrent octopamine.
% Collect stats: median, mean, and std of FR for each digit, pre- and post-training.
% Throughout, digits may be referred to as odors, or as odor puffs.
% 'Pre' = naive. 'Post' = post-training
%
% Inputs:
%   1. simResults: output of sdeWrapper_fn.m
%   2. modelParams: struct of this moth's parameters
%   3. experimentParams: struct with timing and digit class info from the experiment.
%   4. showPlots: 1 x 2 vector. First entry: show changes in accuracy, multiple plots. 
%                 2nd entry: show EN timecourses.
%   5. classLabels: 1 to 10
%   6. resultsFilename:  to generate image filenames if saving. Optional argin
%   7. saveImageFolder:  where to save. If this = [], images will not be saved. Optional argin.
%
% Outputs (as fields of resultsStruct):
%   1. preTrain1MeanResp = numENs x numOdors matrix = mean of EN responses pre-training
%   2. preTrain1StdResp = numENs x numOdors matrix = std of EN responses pre-training
%   1. postTrain1MeanResp = numENs x numOdors matrix = mean of EN responses post-training1
%   2. postTrain1StdResp = numENs x numOdors matrix = std of EN responses  post-training1
%   1. postInjuryMeanResp = numENs x numOdors matrix = mean of EN responses  post injury
%   2. postInjuryStdResp = numENs x numOdors matrix = std of EN responses  post injury
%   1. postTrain2MeanResp = numENs x numOdors matrix = mean of EN responses  post-training2
%   2. postTrain2StdResp = numENs x numOdors matrix = std of EN responses  post-training2
%   3. ditto for post etc
%   4. percentChangeInMeanResp = 1 x numOdors vector
%   5. trained = list of indices corresponding to the odor(s) that were trained
%   6. preTrain1SpontMean = mean(preTrain1Spont);
%   7. preTrain1SpontStd = std(preTrain1Spont);
%   8. postTrain1SpontMean = mean(postTrain1Spont);
%   9. postTrain1SpontStd = std(postTrain1Spont);

% Copyright (c) 2018-2020 Charles B. Delahunt.  delahunt@uw.edu
% MIT License

%-----------------------------------------------------

if nargin < 7, saveImageFolder = []; end
if nargin < 6, resultsFilename = []; end

if ~isempty(saveImageFolder)
    if ~exist( saveImageFolder, 'dir')
        mkdir( saveImageFolder )
    end
end

nE = modelParams.nE;

% various timepoints for sequestering groups of EN responses:
stopSpontMean3 = experimentParams.stopSpontMean3; 
simStop = experimentParams.simStop; 
injuryTime = experimentParams.injuryTime;

startTrain1 = experimentParams.startTrain1;
endTrain1 = experimentParams.endTrain1;
startTrain2 = experimentParams.startTrain2;
endTrain2 = experimentParams.endTrain2; 

preTrain1SpontStart = experimentParams.preTrain1SpontStart;
preTrain1SpontStop = experimentParams.preTrain1SpontStop;
postTrain1SpontStart = experimentParams.postTrain1SpontStart;
postTrain1SpontStop = experimentParams.postTrain1SpontStop;
postInjurySpontStart = experimentParams.postInjurySpontStart;
postInjurySpontStop  = experimentParams.postInjurySpontStop;
postTrain2SpontStart = experimentParams.postTrain2SpontStart;
postTrain2SpontStop = experimentParams.postTrain2SpontStop; 

% 1. data from experimentParams:
simStart = experimentParams.simStart;
classMags = experimentParams.classMags;
stimStarts = experimentParams.stimStarts; % to get timeSteps from very start of sim
stimStarts = stimStarts.*(classMags > 0);  % ie only use non-zero puffs
whichClass = experimentParams.whichClass;
whichClass = whichClass.*(classMags > 0);
hebStarts = experimentParams.hebStarts;

classList = sort(unique(whichClass));
numClasses = length(classList);

E = simResults.E;   % # timesteps x #ENs
T = simResults.T;   % # timesteps x 1
octoHits = simResults.octoHits;
if max(octoHits) > 0
    octoTimes = T(octoHits > 0);
else
    octoTimes = [];
end

%% 

colors = { 'r', 'b'  }; %  for 2 classes
% concurrent octopamine will be marked with yellow x's.


% calc spont stats:
preTrain1Spont = E( T> preTrain1SpontStart & T < preTrain1SpontStop );
postTrain1Spont = E( T > postTrain1SpontStart & T < postTrain1SpontStop );
postInjurySpont = E( T > postInjurySpontStart & T < postInjurySpontStop );
postTrain2Spont =  E( T > postTrain2SpontStart & T < postTrain2SpontStop );
preTrain1SpontMean = mean(preTrain1Spont);
preTrain1SpontStd = std(preTrain1Spont);
postTrain1SpontMean = mean(postTrain1Spont);
postTrain1SpontStd = std(postTrain1Spont);
postInjurySpontMean = mean(postInjurySpont);
postInjurySpontStd = std(postInjurySpont);
postTrain2SpontMean = mean(postTrain2Spont);
postTrain2SpontStd = std(postTrain2Spont);

% Set regions to examine:


% Make one stats plot per EN. Loop through ENs (but there is only one EN):

for enInd = 1:nE
    thisEnResponse = E(:, enInd );
    %% Calculate pre- and post-train odor response stats:
    % Assumes that there is at least 1 sec on either side of an odor without octo
    
    for i = 1:length(stimStarts)
        t = stimStarts(i);
        % Note: to find no-octo stimStarts, there is a certain amount of machinery in order to mesh 
        % with the timing data from the experiment. For some reason octoTimes are not stored 
        % exactly as listed in format short mode. So we need to use abs difference > small thresh, 
        % rather than ~ismember(t, octoTimes):
        small = 1e-8;
        
        % assign no-octo, PRE-train1 response val (or -1):
        preTrain1OdorResp(i) = -1;  % as flag
        if min(abs(octoTimes - t)) > small && t < startTrain1
            preTrain1OdorResp(i) = max( thisEnResponse ( T > t-1 & T < t+1 ) ); % *Resp = max EN 
                                                            % value within +/- 1 sec of odor start
        end
        % assign no-octo, POST-train1 response val (or -1):
        postTrain1OdorResp(i) = -1;
        if ~isempty(octoTimes)
            if min(abs(octoTimes - t)) > small && t > endTrain1 && t < injuryTime
                postTrain1OdorResp(i) = max( thisEnResponse ( T > t-1 & T < t+1 ) );
            end
        end
        % assign no-octo, POST-injury response val (or -1):
        postInjuryOdorResp(i) = -1;
        if ~isempty(octoTimes)
            if min(abs(octoTimes - t)) > small && t > injuryTime && t < startTrain2
                postInjuryOdorResp(i) = max( thisEnResponse ( T > t-1 & T < t+1 ) );
            end
        end
        % assign no-octo, POST-train2 response val (or -1):
        postTrain2OdorResp(i) = -1;
        if ~isempty(octoTimes)
            if min(abs(octoTimes - t)) > small && t > endTrain2  
                postTrain2OdorResp(i) = max( thisEnResponse ( T > t-1 & T < t+1 ) );
            end
        end
    end
    
    % calc no-octo stats for each odor, pre and post train: There are 2 classes (odors)
    for k = 1:numClasses
        preTrain1Full = preTrain1OdorResp( preTrain1OdorResp >= 0 & whichClass == classList(k) );   
        % ">= 0" culls the -1 values created above.
        postTrain1Full = postTrain1OdorResp( postTrain1OdorResp >= 0 & whichClass == classList(k) );
        postInjuryFull = postInjuryOdorResp( postInjuryOdorResp >= 0 & whichClass == classList(k) );
        postTrain2Full = postTrain2OdorResp( postTrain2OdorResp >= 0 & whichClass == classList(k) ); 
        
        % calculate the averaged sniffs of each sample: SA means 'sniffsAveraged'
        preTrain1SA = zeros(1, length( preTrain1Full ) );  % this will contain the average 
                                                        % responses over all sniffs for each sample.
        for i = 1:length(preTrain1SA)
            preTrain1SA(i) = mean( preTrain1Full( (i-1) + 1 : i ) );
        end
        postTrain1SA = zeros( 1, length(postTrain1Full) );
        for i = 1:length(postTrain1SA)
            postTrain1SA(i) = mean( postTrain1Full( (i-1) + 1 : i ) );
        end        
        postInjurySA = zeros( 1, length(postInjuryFull) );
        for i = 1:length(postInjurySA)
            postInjurySA(i) = mean( postInjuryFull( (i-1) + 1 : i ) );
        end
        postTrain2SA = zeros( 1, length(postTrain2Full) );
        for i = 1:length(postTrain2SA)
            postTrain2SA(i) = mean( postTrain2Full( (i-1) + 1 : i ) );
        end
        
        if isempty(preTrain1SA)
            preTrain1MeanResp(k) = -1;
            preTrain1MedianResp(k) = -1;
            preTrain1StdResp(k) = -1;
            preTrain1NumPuffs(k) = 0;
        else
            preTrain1MeanResp(k) = mean(preTrain1SA);
            preTrain1MedianResp(k) = median(preTrain1SA);
            preTrain1StdResp(k) = std(preTrain1SA);
            preTrain1NumPuffs(k) = length(preTrain1SA);
        end
        if isempty(postTrain1SA)
            postTrain1MeanResp(k) = -1;
            postTrain1MedianResp(k) = -1;
            postTrain1StdResp(k) = -1;
            postTrain1NumPuffs(k) = 0;
        else
            postTrain1MeanResp(k) = mean(postTrain1SA);
            postTrain1MedianResp(k) = median(postTrain1SA);
            postTrain1StdResp(k) = std(postTrain1SA);
            postTrain1NumPuffs(k) = length(postTrain1SA);
        end
        if isempty(postInjurySA)
            postInjuryMeanResp(k) = -1;
            postInjuryMedianResp(k) = -1;
            postInjuryStdResp(k) = -1;
            postInjNumPuffs(k) = 0;
        else
            postInjuryMeanResp(k) = mean(postInjurySA);
            postInjuryMedianResp(k) = median(postInjurySA);
            postInjuryStdResp(k) = std(postInjurySA);
            postInjNumPuffs(k) = length(postInjurySA);
        end
        if isempty(postTrain2SA)
            postTrain2MeanResp(k) = -1;
            postTrain2MedianResp(k) = -1;
            postTrain2StdResp(k) = -1;
            postTrain2NumPuffs(k) = 0;
        else
            postTrain2MeanResp(k) = mean(postTrain2SA);
            postTrain2MedianResp(k) = median(postTrain2SA);
            postTrain2StdResp(k) = std(postTrain2SA);
            postTrain2NumPuffs(k) = length(postTrain2SA);
        end
    end % for k
    
    % to plot +/- 1 std of % change in meanResp, we want the std of our
    % estimate of the mean = stdResp/sqrt(numPuffs). Make this calc:
    preTrain1StdMeanEst = preTrain1StdResp./sqrt(preTrain1NumPuffs);
    postTrain1StdMeanEst = postTrain1StdResp./sqrt(postTrain1NumPuffs);
    postInjuryStdMeanEst = postInjuryStdResp./sqrt(postInjNumPuffs);
    postT2StdMeanEst = postTrain2StdResp./sqrt(postTrain2NumPuffs);
     
    % for plotting (x axis), easy since we have two odor only.
    xInds = 1:2;
    postOffset = xInds + 0.25;
    
    % change in response due to Train1: for each of the two odors
    percentChangeDueToTrain1InMeanResp = 100*(postTrain1MeanResp - preTrain1MeanResp )./ ...
        preTrain1MeanResp;
    percentChangeDueToTrain1InNoiseSubtractedMeanResp =...
        100*(postTrain1MeanResp - preTrain1MeanResp - preTrain1SpontMean)./preTrain1MeanResp;    
    percentChangeDueToTrain1InMedianResp = 100*(postTrain1MedianResp - preTrain1MedianResp )./ ...
        preTrain1MedianResp;
    percentChangeDueToTrain1InNoiseSubtractedMedianResp =...
        100*(postTrain1MedianResp - preTrain1MedianResp - preTrain1MedianResp)./preTrain1MedianResp;
    % due to injury (taking post-train1 as baseline):
    percentChangeDueToInjuryInMeanResp = 100*(postInjuryMeanResp - postTrain1MeanResp )./ ...
        postTrain1MeanResp;
    percentChangeDueToInjuryInNoiseSubtractedMeanResp =...
        100*(postInjuryMeanResp - postTrain1MeanResp - postTrain1SpontMean)./postTrain1MeanResp;    
    percentChangeDueToTrain1InMedianResp = 100*(postTrain1MedianResp - postTrain1MedianResp )./ ...
        postTrain1MedianResp;
    percentChangeDueToTrain1InNoiseSubtractedMedianResp =...
        100*(postTrain1MedianResp - postTrain1MedianResp - postTrain1MedianResp)./postTrain1MedianResp;
    % due to Train2 (taking post-injury as baseline):
    percentChangeDueToTrain2InMeanResp = 100*(postTrain2MeanResp - postInjuryMeanResp )./ ...
        postInjuryMeanResp;
    percentChangeDueToTrain2InNoiseSubtractedMeanResp =...
        100*(postTrain2MeanResp - postInjuryMeanResp - postInjurySpontMean)./postInjuryMeanResp;    
    percentChangeDueToTrain2InMedianResp = 100*(postTrain2MedianResp - postInjuryMedianResp )./ ...
        postInjuryMedianResp;
    percentChangeDueToTrain2InNoiseSubtractedMedianResp =...
        100*(postTrain2MedianResp - postInjuryMedianResp - postInjuryMedianResp)./ ...
        postInjuryMedianResp;
    % due to Injury and Train2 (taking post-Train1 as baseline)
    percentChangeDueToInjuryAndTrain2InMeanResp =...
        100*(postTrain2MeanResp - postTrain1MeanResp )./postTrain1MeanResp;
    percentChangeDueToInjuryAndTrain2InNoiseSubtractedMeanResp =...
        100*(postTrain2MeanResp - postTrain1MeanResp - postTrain1SpontMean)./postTrain1MeanResp;    
    percentChangeDueToInjuryAndTrain2InMedianResp =...
        100*(postTrain2MedianResp - postTrain1MedianResp )./postTrain1MedianResp;
    percentChangeDueToInjuryAndTrain2InNoiseSubtractedMedianResp =...
        100*(postTrain2MedianResp - postTrain1MedianResp - postTrain1MedianResp)./ ...
        postTrain1MedianResp;
    % create cell of xticklabels:
    trueXLabels = cell(size(classLabels));
    for j = 1:length(classLabels)
        trueXLabels{j} = num2str(mod(classLabels(j),10) );   % the 'mod' turns 10 into 0
    end
    %% plot stats if wished:
%     if showPlots(1)
%         % Removed.        
%     end % if showPlots
 % Save plot code:
    if ~isempty(saveImageFolder) && showPlots(1)
        saveas( thisFig, fullfile(saveImageFolder,[resultsFilename '_en' num2str(enInd),...
            '.png']), 'png')
    end
    
    %---------------------------------------------------------------------------------
    
    % store results in a struct:
    r(enInd).preTrain1OdorResp = preTrain1OdorResp;  % preserves all the sniffs for each stimulus
    r(enInd).postTrain1OdorResp = postTrain1OdorResp;
    r(enInd).postInjuryOdorResp = postInjuryOdorResp;
    r(enInd).postTrain2OdorResp = postTrain2OdorResp;
    r(enInd).preTrain1SniffsAved = preTrain1SA;   % the averaged sniffs for each stimulus
    r(enInd).postTrain1SniffsAved = postTrain1SA;
    r(enInd).postTrain2SniffsAved = postTrain2SA;
    r(enInd).postInjurySniffsAved = postInjurySA;
    r(enInd).odorClass = whichClass;
    
    r(enInd).percentChangeDueToTrain1InMeanResp = percentChangeDueToTrain1InMeanResp;  % key stat
    r(enInd).percentChangeDueToInjuryInMeanResp = percentChangeDueToInjuryInMeanResp;  % key stat
    r(enInd).percentChangeDueToTrain2InMeanResp = percentChangeDueToTrain2InMeanResp;  % key stat
    r(enInd).percentChangeDueToInjuryAndTrain2InMeanResp = ...
        percentChangeDueToInjuryAndTrain2InMeanResp;  % key stat
    
    r(enInd).percentChangeDueToTrain1InNoiseSubtractedMeanResp = ...
        percentChangeDueToTrain1InNoiseSubtractedMeanResp;
    r(enInd).relativeChangeDueToTrain1InNoiseSubtractedMeanResp = ...
        percentChangeDueToTrain1InNoiseSubtractedMeanResp / ...
        percentChangeDueToTrain1InNoiseSubtractedMeanResp(enInd);
    r(enInd).percentChangeDueToTrain1InMedianResp = percentChangeDueToTrain1InMedianResp;
    r(enInd).percentChangeDueToTrain1InNoiseSubtractedMedianResp = ...
        percentChangeDueToTrain1InNoiseSubtractedMedianResp;
    r(enInd).relativeChangeDueToTrain1InNoiseSubtractedMedianResp = ...
        ( (postTrain1MedianResp - preTrain1MedianResp - preTrain1MedianResp )./preTrain1MedianResp )...
        / ( (postTrain1MedianResp(enInd) - preTrain1MedianResp(enInd) - preTrain1MedianResp )./ ...
        preTrain1MedianResp(enInd) );
    
    r(enInd).trained = enInd;
    
    % The following variables are the most important:
    
    % EN odor responses, pre and post training.
    % Vectors of length numClasses (number of distinct odors. First entry is trained odor)
    r(enInd).preTrain1MeanResp = preTrain1MeanResp;
    r(enInd).preTrain1StdResp = preTrain1StdResp;
    r(enInd).postTrain1MeanResp = postTrain1MeanResp;
    r(enInd).postTrain1StdResp = postTrain1StdResp;
    r(enInd).postInjuryMeanResp = postInjuryMeanResp;
    r(enInd).postInjuryStdResp = postInjuryStdResp;
    r(enInd).postT2MeanResp = postTrain2MeanResp;
    r(enInd).postT2StdResp = postTrain2StdResp;
    % spont responses, pre and post training:
    r(enInd).preTrain1SpontMean = mean(preTrain1Spont);
    r(enInd).preTrain1SpontStd = std(preTrain1Spont);
    r(enInd).postTrain1SpontMean = mean(postTrain1Spont);
    r(enInd).postTrain1SpontStd = std(postTrain1Spont);
    r(enInd).postInjurySpontMean = mean(postInjurySpont);
    r(enInd).postInjurySpontStd = std(postInjurySpont);
    r(enInd).postTrain2SpontMean = mean(postTrain2Spont);
    r(enInd).postTrain2SpontStd = std(postTrain2Spont);
    
    results = r;
    
end % for enInd = 1:numClasses

%% Plot distance between odor 1 and odor 2 at each stage

if showPlots(1)
    colorList = { 'b', 'r', 'g' };
    figure
        subplot(2,1,1)
        hold on
        for i = 1:length(preTrain1MeanResp)
            thisColor = colorList{i};
            errorbar( [ preTrain1MeanResp(i), postTrain1MeanResp(i), postInjuryMeanResp(i),...
                postTrain2MeanResp(i)], [preTrain1StdResp(i), postTrain1StdResp(i),...
                postInjuryStdResp(i), postTrain2StdResp(i) ], thisColor )
        end
        % now plot spontaneous noise:
        thisColor = colorList{3};
        errorbar( [ preTrain1SpontMean, postTrain1SpontMean, postInjurySpontMean,...
            postTrain2SpontMean], [ preTrain1SpontStd, postTrain1SpontStd, postInjurySpontStd,...
            postTrain2SpontStd], thisColor )
         
        title(['Distance between trained odor and control at each stage, ',...
               strrep(resultsFilename,'_','-')])
        
        % percent changes:
        denom = 0.01*preTrain1MeanResp;
        subplot(2,1,2)
        hold on
        for i = 1:length(preTrain1MeanResp)
            thisColor = colorList{i};
            errorbar( [ preTrain1MeanResp(i) / denom(i), postTrain1MeanResp(i) / denom(i),...
                        postInjuryMeanResp(i) / denom(i), postTrain2MeanResp(i) / denom(i) ],...
                      [ preTrain1StdResp(i) / denom(i), postTrain1StdResp(i) / denom(i),...
                      postInjuryStdResp(i) / denom(i), postTrain2StdResp(i) / denom(i) ], thisColor)
        end
        % now plot spontaneous noise:
        thisColor = colorList{3};
        spontDenom = 0.01*preTrain1SpontMean;
        errorbar([ preTrain1SpontMean / spontDenom, postTrain1SpontMean / spontDenom,...
            postInjurySpontMean / spontDenom, postTrain2SpontMean / spontDenom ],...
            [ preTrain1SpontStd / spontDenom, postTrain1SpontStd / spontDenom ,...
                postInjurySpontStd / spontDenom , postTrain2SpontStd / spontDenom  ], thisColor )
        legend('trained odor', 'control odor','spont noise')
        xlabel('naive, post-train1, post-injury, post-train2' )
        ylabel('mean EN response' )
        ylim([ -50, 300])
        title('Percent changes for main odor, control, and spontaneous noise')
         
end        

if showPlots(1)
    plotDeltaMeansFn(preTrain1MedianResp, postTrain1MedianResp, preTrain1MeanResp, ...
        postTrain1MeanResp, preTrain1StdResp, postTrain1StdResp, preTrain1SpontMean,...
        postTrain1SpontMean, 'due to Train1' )
    plotDeltaMeansFn(postTrain1MedianResp, postInjuryMedianResp, postTrain1MeanResp,...
        postInjuryMeanResp, postTrain1StdResp, postInjuryStdResp, postTrain1SpontMean,...
        postInjurySpontMean, 'due to Injury' )
    plotDeltaMeansFn(postInjuryMedianResp, postTrain2MedianResp, postInjuryMeanResp,...
        postTrain2MeanResp, postInjuryStdResp, postTrain2StdResp, postInjurySpontMean,...
        postTrain2SpontMean, 'due to Train2' )
    plotDeltaMeansFn(postTrain1MedianResp, postTrain2MedianResp, postTrain1MeanResp,...
        postTrain2MeanResp, postTrain1StdResp, postTrain2StdResp, postTrain1SpontMean,...
        postTrain2SpontMean, 'due to Injury and Train2' )
 
 % Save plot code:
    if ~isempty(saveImageFolder)
        saveas( thisFig, fullfile(saveImageFolder,[resultsFilename, '_en', num2str(enInd),...
            '.png']), 'png')
    end
end
%% Plot EN timecourses normalized by mean digit response:
 
% labels = whichClass;

if showPlots(2)
    scrsz = get(0,'ScreenSize');
    % go through each EN: 
    for enInd = 1:nE           % recal EN1 targets digit class 1, EN2 targets digit class 2, etc
        % NOTE: since there is only one EN, we take a notational shortcut and use r.field instead 
        % of r(enInd).field
        enFig2 = figure('Position',[scrsz(1), scrsz(2), scrsz(3)*0.8, scrsz(4)*0.4 ]); % make a new figure  
  
        hold on
        
        % plot timecourse of EN:
        zeroTimeInd = find(T >= 0, 1);
        plot(T(zeroTimeInd:end), E(zeroTimeInd:end) )
        
        plot(injuryTime, 0, 'r.','markersize', 20)
        plot([injuryTime, injuryTime], [0, 0.1],'r', 'linewidth', 2)
        
        % mark stims: 
        inds = find(whichClass == 1);
        plot( stimStarts(inds), zeros(size(inds)), 'b.','markersize', 20 )
        
        inds = find(whichClass == 2);
        plot( stimStarts(inds), zeros(size(inds)), 'g.','markersize', 20 )
         
        plot( hebStarts, 0.01*ones(size(hebStarts)),'y.','markersize', 20 ) 
        
        % Labels:
        title([ strrep(resultsFilename,'_','-'), ' EN timecourse'])
        xlabel(['Time. Blue dots = main odor puffs, green dots = control puffs, ', ...
            'red line = time of injury, yellow dots = training puffs'],'fontsize', 14,...
            'fontweight', 'bold')
        ylabel('EN response (dimensionless)','fontsize', 14, 'fontweight', 'bold')
 
 % Save EN timecourse:
        if ~isempty(saveImageFolder)
            saveas( enFig2, fullfile(saveImageFolder,[resultsFilename '_enTimecourses',...
                num2str(enInd) '.png']), 'png')
        end

    end  % for enInd = 1:nE
    
end % if showPlots(2)

function [] = plotDeltaMeansFn( preMedian, postMedian, preMean, postMean, preStd, postStd,...
                                preSpont, postSpont, dataLabel )

        xInds = [1, 2];
        postOffset = xInds + 0.25;
        enInd = 1;
        
        del = texlabel('Delta');
        scrsz = get(0,'ScreenSize');
        thisFig = figure('Position',[scrsz(1), scrsz(2), scrsz(3)*0.4, scrsz(4)*0.4 ]);
        % medians
        subplot(2,3,1)
        hold on,
        grid on,
        
        % raw median, pre and post:
        plot(xInds, preMedian,'*b')
        plot(postOffset, postMedian,'bo','markerfacecolor','b')
        %     plot(pre, preMeanResp + preStdResp,'+g')
        %     plot(post, postMeanResp + postStdResp,'+g')
        %     plot(pre, preMeanResp - preStdResp,'+g')
        %     plot(post, postMeanResp - postStdResp,'+g')
        
        % make the home EN of this plot red:
        plot(enInd, preMedian(enInd), 'ro')
        plot(enInd + 0.25, postMedian(enInd), 'ro','markerfacecolor','r')
        title(['EN ' num2str(enInd), ' median +/- std' ' ' dataLabel ])
        xlim([0,max(xInds) + 1])
        ylim( [0  1.1*max([ preMedian, postMedian ]) ] )
        
        % connect pre to post with lines for clarity:
        for j = xInds
            lineColor = 'b';
            if j == 1, lineColor = 'r'; end
            plot( [ xInds(j),postOffset(j) ], [ preMedian(j), postMedian(j) ], lineColor )
        end
        
        % percent change in medians:
        subplot(2,3,2)
        hold on,
        grid on,
        plot(xInds, 100*(postMedian - preMedian )./preMedian,...
            'bo','markerfacecolor','b')
        % mark the trained odors in red:
        plot(enInd, 100*(postMedian(enInd) - preMedian(enInd) )./preMedian(enInd),...
            'ro','markerfacecolor','r')
        title( [ '% ' del ' median' ])
        xlim([0,max(xInds) + 1])
        % ylim([-50,400])
        
        % relative changes in median, ie control/trained:
        subplot(2,3,3)
        pn = sign( postMedian(enInd) - preMedian(enInd) );
        hold on,
        grid on,
        plot(1:2, pn*( (postMedian - preMedian )./preMedian) / ...
            ( (postMedian(enInd) - preMedian(enInd) )./preMedian(enInd) ),...
            'bo','markerfacecolor','b')
        % mark the trained odors in red:
        plot(enInd, pn*1, 'ro','markerfacecolor','r')
        title(['relative ' del ' median relative to ', del, 'trainedOdorMedian'])
        xlim([0,max(xInds) + 1])
        % ylim([0,2])
        
        %------------------------------------------------------------------------
        % means
        % raw means, pre and post:
        subplot(2,3,4)
        hold on,
        grid on,
        errorbar(xInds, preMean, preStd, 'bo')
        errorbar(postOffset, postMean, postStd,'bo','markerfacecolor','b')
        errorbar(enInd, preMean(enInd), preStd(enInd), 'ro')
        errorbar(enInd + 0.25, postMean(enInd),  postStd(enInd), 'ro','markerfacecolor','r')
        title(['EN ' num2str(enInd), ' mean +/- std'])
        xlim([0,max(xInds) + 1])
        ylim( [0  1.1*max([ preMean, postMean ]) + max([preStd, postStd]) ] )
        
        % plot spont:
        errorbar(xInds(1), preMean(1), preStd(1), 'mo')
        errorbar(postOffset(1), postMean(1), postStd(1),'mo','markerfacecolor','m')
        
        % percent change in means
        subplot(2,3,5)
        hold on,
        grid on,
        plot(xInds, 100*postMean./preMean,'bo','markerfacecolor','b')
        % mark the trained odors in red:
        plot(enInd, 100*postMean(enInd)/preMean(enInd),'ro','markerfacecolor','r')
        title([ '% ' del ' mean' ] )
        xlim([0,max(xInds) + 1])
        
        % relative percent changes
        subplot(2,3,6)
        pn = sign( postMean(enInd) - preMean(enInd) );
        hold on,
        grid on,
        plot(xInds, pn*postMean  / preMean(enInd),'bo','markerfacecolor','b')
        % mark the trained odors in red:
        plot(enInd, pn, 'ro','markerfacecolor','r')
        title(['relative ' del ' mean relative to ', del, 'trainedOdorMean'])
        xlim([0,max(xInds) + 1])
% end of plotDeltaMeansFn



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




















