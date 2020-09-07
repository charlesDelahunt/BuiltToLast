% filtering function based on Pedro's work:
function [fr_lp] = pedroLowPassFilterForRNs(frIn, scalingFactor, frMax )
% inputs: frIn = pre-filter firing rates (vector)
%	      scalingFactor = amount to multiply frIn by to bring it into a good regime for Pedro's 
%						  formula: If frIn*scalingFactor < 5 -> LP has approx no effect. If > 16,
%						  LP filter effect decreases (ie peak of quadratic is at approx 16).
%         frMax = maximum FR for the type of neuron being filtered. For PNs, this
%                 could be 1/dt. Note that the initial value of frMax gets
%                 updated once stimuli start hitting, so its effect in this
%                 function is to make input FRs ~ [0, 1] (then the
%                 scalingFactor makes frSc ~ [0, scalingFactor ] ).
% outputs: fr_lp = low pass filtered outputs (vector).

% Key point: since frIn = RNs, we assume there are 500 individual axons, so
% the sigma variation will average out in the summed total. So set sig = 0.

% method: scale frIn to [0, 10]. then pass it through a quadratic, add a
% variance term, divide by the scaling term (eg 10), and cap upper and
% lower levels.

% CAUTION: behavior becomes aberrant if frSc > 15

% Copyright (c) 2018-2020 Charles B. Delahunt.  delahunt@uw.edu
% MIT License

%--------------------------------------------------------------------
 
% Compute low-passed FR:
sig   = 0; % we assume that over many antennae input neurons the random variations cancel out
frSc = frIn*scalingFactor/frMax;  % scaled firing rates
fr_lpSc = - 0.029*frSc.^2 + 0.93*frSc - 0.21 + sig*randn(size(frSc)); % low-passed scaled FRs
fr_lp = fr_lpSc/scalingFactor*frMax; % low-passed FRs (after reversal of the scaling)
fr_lp = max(fr_lp, 0);
fr_lp = min(fr_lp, frIn);


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
 
