function [fA] = odorAndAddedNoiseCreatePuffs_fn(nG, numPuffs, odorSpread, odorStrengthStd,...
                                                inputNoiseLevel, inputNoiseSpread, normalizeMags)

% generates random odor projections, sufficient puffs for the experiment,
% for two classes: (i) odor + added noise, (ii) noise only
% Inputs:
%     1. nG = number of glomeruli
%     2. numPuffs = how many puffs of the trained odor we need
%     3. odorSpread = fraction of gloms the odor will hit on average
%     4. odorStrengthStd = std dev of odor strength, AND std of noise strength
%     5. inputNoiseLevel = what fraction of odor magnitude should noise magnitude be. 
%     6. inputNoiseSpread = fraction of gloms the noise will hit (the particular gloms hit change 
%        with each puff)
%     7. normalizeMags =  dummy, since there is only one actual odor. 
%     
%     Output:
%     1. fA = nG x numPuffs x numOdors array
% 
%     Odor inputs to each RN have magnitude 1 +/- std
%     Noise inputs have magnitude 'noiseLevel' +/- std*noiseLevel

% Copyright (c) 2018-2020 Charles B. Delahunt.  delahunt@uw.edu
% MIT License

%--------------------------------------------------------------------
 
numOdors = 2;

odor1Targets = rand(nG, 1) < odorSpread; % 1s and 0s.

% modify to give varying strengths of odor to Gloms:
odor1ToG = ( odor1Targets + odorStrengthStd*randn(size(odor1Targets)) ).*odor1Targets; % the last 
                                                                         % term ensures 0s stay 0s
odor1ToG = max(0, odor1ToG); % To prevent any negative weights odor1ToG is a nG x 1 col matrix, with
                             % some zero entries and some non-zero entries which have mean = 1

% make two instances of added noise as a matrix nG x numPuffs. one instance gets
% added to odor1, one instance becomes odor2:
noisePuffTargets = rand(nG, numPuffs) < inputNoiseSpread; 
noisePuffs = inputNoiseLevel*ones(nG, numPuffs) + ...
             odorStrengthStd*inputNoiseLevel*randn(nG, numPuffs);
noisePuffs = noisePuffs.*noisePuffTargets;

odor1 = odor1ToG + noisePuffs;

noisePuffTargets = rand(nG, numPuffs) < inputNoiseSpread; 
noisePuffs = inputNoiseLevel*ones(nG, numPuffs) + ...
             odorStrengthStd*inputNoiseLevel*randn(nG, numPuffs);
noisePuffs = noisePuffs.*noisePuffTargets;

odor2 = noisePuffs;

fA = zeros(nG, numPuffs, numOdors);
fA(:, :, 1) = odor1;
fA(:, :, 2) = odor2;



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