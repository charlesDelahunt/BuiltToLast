function [fA] = twoOdorsCreatePuffs_fn(nG, numPuffs, odorSpread, odorStrengthStd,...
                                       inputNoiseLevel, inputNoiseSpread, normalizeMags)

% generates random odor projections, sufficient puffs for the experiment,
% for two classes: (i) odor1 (which gets trained, (ii) odor2
% Inputs:
%     1. nG = number of glomeruli
%     2. numPuffs = how many puffs of the trained odor we need
%     3. odorSpread = fraction of gloms the odor will hit on average
%     4. odorStrengthStd = std dev of odor strength, AND std of noise strength
%     5. inputNoiseLevel = what fraction of odor magnitude should noise magnitude be. 
%     6. inputNoiseSpread = fraction of gloms the noise will hit (the particular gloms hit change 
%        with each puff)
%     7. normalizeMags = 0 or 1. If 1, normalize summed input (into AL) mag of odor2  to match mag 
%        of odor 1. 
%     
%     Output:
%     1. fA = nG x numPuffs x numOdors array

% All odors projection onto gloms has mean magnitude = 1, unless odor2 gets normalized.

% Copyright (c) 2018-2020 Charles B. Delahunt.  delahunt@uw.edu
% MIT License

%-----------------------------------------------------

numOdors = 2;

odor1Targets = rand(nG, 1) < odorSpread; % 1s and 0s.
% modify to give varying strengths of odor to Gloms:
odor1ToG = ( odor1Targets + odorStrengthStd*randn(size(odor1Targets)) ).*odor1Targets; % the last 
% term ensures 0s stay 0s
odor1 = max(0, odor1ToG); % to prevent any negative weights

odor2Targets = rand(nG, 1) < odorSpread; % 1s and 0s.     
% modify to give varying strengths of odor to Gloms:
odor2ToG = ( odor2Targets + odorStrengthStd*randn(size(odor2Targets)) ).*odor2Targets; % the last 
% term ensures 0s stay 0s
odor2 = max(0, odor2ToG); % to prevent any negative weights

if normalizeMags
    odor2 = odor2*( sum(odor1ToG) / sum(odor2ToG) );
end

fA = zeros(nG, numPuffs, numOdors);
fA(:, :, 1) = repmat(odor1, [1, numPuffs, 1] );
fA(:, :, 2) = repmat(odor2, [1, numPuffs, 1] );


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

