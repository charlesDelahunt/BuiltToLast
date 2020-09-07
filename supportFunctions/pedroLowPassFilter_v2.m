% low-pass filter based on Pedro's work:
function [fr_lp] = pedroLowPassFilter_v2(frRaw, params )
% inputs: 
%   1. frIn = pre-filter firing rates
%   2. params has three fields:
%             a. scalingFactor = to get inputs into [0, 10] range.
%             b. maxFR = maximum FR for the type of neuron being filtered. For
%                      PNs, empirical value  (2017 experiments) is 0.7 to 0.8
%             c. randnVector: vector of randn(size(LP-affected neurons)), fixed at time of injury.
%                        
% outputs: 
%   1. fr_lp = low pass filtered outputs.

% method: scale frIn to [0, 10]. then pass it through a quadratic, add a
% variance term, and cap upper and lower levels:

% Copyright (c) 2018-2020 Charles B. Delahunt.  delahunt@uw.edu
% MIT License

%--------------------------------------------------------------------
 
scalingFactor = params.scalingFactor;
maxFR = params.maxFR;
randnVector = params.randnVector;

% Compute low-passed FR:
sig   = 0.4611;  % from Pedro
frSc = frRaw * scalingFactor / maxFR;  % scaled to [0, 10 to 15]
fr_lpSc = - 0.029*frSc.^2 + 0.93*frSc - 0.21 + sig*randnVector;
fr_lp = fr_lpSc / scalingFactor * maxFR;
fr_lp = max(fr_lp, 0);
fr_lp = min(fr_lp, frRaw);


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