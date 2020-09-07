function [ injuryVectors ] = createInjuryFlags_1( modelParams, injuryParams)
%
% Describes an injury regime.
% inputs:
%   modelParams = 'modelParams' after connection matrices have been
%   constructed.
%   injuryParams = placeholder argin
% The function will internally specify injuryPieChart = 1 x 3 vector, sum <= 1, for each neuron type
% affected. 
%   first entry = fraction of neurons to be ablated; 2nd = fraction to be
%   reflected; 3rd = fraction to be filtered. 
%   NOTE: these fractions don't sum to 1, because healthy neurons are not listed (default = healthy)
%
% Outputs:
%   injuryVectors: struct with fields rn_fasFlags, pn_fasFlags, etc, as needed.

% *_fasFlags: vector of same length as the neuron type's FR vector, with ints:
%       0 = no injury:      xOut = xIn
%       1 = ablation:       xOut = 0
%       2 = reflection:     xOut = 0.5*xIn
%       3 = fancy filter:   xOut = @fancyHandle(xIn, params) 

% Copyright (c) 2018 - 2020, Charles B. Delahunt.  delahunt@uw.edu
% MIT License

%-------------------------------------------------

% initialize injuryFlag vectors:  0 entries = healthy (default)
rn_fasFlags = zeros(1, modelParams.nR);
pn_fasFlags = zeros(1, modelParams.nP);
pin_fasFlags = zeros(1, modelParams.nPI);
ln_fasFlags = zeros(1, modelParams.nG);
kc_fasFlags = zeros(1, modelParams.nK);
fancyFilterHandle = [];
fancyFilterParams = [];

%-------------------------------------------------------------------

% define fancy filter handle and params here:
if isfield(injuryParams,'fancyFilterHandle')
    fancyFilterHandle = injuryParams.fancyFilterHandle;
end
if isfield(injuryParams,'fancyFilterParams')
    fancyFilterParams = injuryParams.fancyFilterParams;
end

%----------------------------------------------------------------

% update each injuryFlag vector in turn
% These flags give the type of injury each neuron suffers. RNs ARE A SPECIAL CASE!
if isfield(injuryParams,'rn_pie')
    % NOTE: RNs are a special case because there are 30k, distributed evenly across the antennae. So
    % RNs should be injured using ONLY t(3) plus a filter Fn = *constant; NOT t(1) and t(2)
    t = injuryParams.rn_pie;          % 1 x 3 vector
    t = cumsum(t);
    v = rand(size(rn_fasFlags));
    rn_fasFlags(v < t(1)) = 1;            % ablate 
    rn_fasFlags(v > t(1) & v < t(2)) = 2;   % reflect
    rn_fasFlags(v > t(2) & v < t(3)) = 3;   % fancy filter
end
    
if isfield(injuryParams,'pn_pie')
    t = injuryParams.pn_pie;          % 1 x 3 vector
    t = cumsum(t);
    v = rand(size(pn_fasFlags));
    pn_fasFlags(v < t(1)) = 1;            % ablate 
    pn_fasFlags(v > t(1) & v < t(2)) = 2;   % reflect
    pn_fasFlags(v > t(2) & v < t(3)) = 3;   % fancy filter
    % Used for low-pass-injury neurons, to give variability in LP injuries:
    fancyFilterParams.pnSigVectorForLP = ( randn(size(find(pn_fasFlags == 3))) )';  % transpose into
                                                                                    % a col vector
end

if isfield(injuryParams,'pin_pie')
    t = injuryParams.pin_pie;          % 1 x 3 vector
    t = cumsum(t);
    v = rand(size(pin_fasFlags));
    pin_fasFlags(v < t(1)) = 1;            % ablate 
    pin_fasFlags(v > t(1) & v < t(2)) = 2;   % reflect
    pin_fasFlags(v > t(2) & v < t(3)) = 3;   % fancy filter
    fancyFilterParams.pinSigVectorForLP = ( randn(size(find(pin_fasFlags == 3))) )';  % transpose
                                                                                % into a col vector
end

if isfield(injuryParams,'ln_pie')
    t = injuryParams.ln_pie;          % 1 x 3 vector
    t = cumsum(t);
    v = rand(size(ln_fasFlags));
    ln_fasFlags(v < t(1)) = 1;            % ablate 
    ln_fasFlags(v > t(1) & v < t(2)) = 2;   % reflect
    ln_fasFlags(v > t(2) & v < t(3)) = 3;   % fancy filter
    fancyFilterParams.lnSigVectorForLP = ( randn(size(find(ln_fasFlags == 3))) )';  % transpose 
                                                                                % into a col vector
end

if isfield(injuryParams,'kc_pie')
    t = injuryParams.kc_pie;          % 1 x 3 vector
    t = cumsum(t);
    v = rand(size(kc_fasFlags));
    kc_fasFlags(v < t(1)) = 1;            % ablate 
    kc_fasFlags(v > t(1) & v < t(2)) = 2;   % reflect
    kc_fasFlags(v > t(2) & v < t(3)) = 3;   % fancy filter
    fancyFilterParams.kcSigVectorForLP = ( randn(size(find(kc_fasFlags == 3))) )';  % transpose into
                                                                                    % a col vector
end

%--------------------------------------------------------------------

% bundle into argout:

injuryVectors.rn_fasFlags = rn_fasFlags;
injuryVectors.pn_fasFlags = pn_fasFlags;
injuryVectors.pin_fasFlags = pin_fasFlags;
injuryVectors.ln_fasFlags = ln_fasFlags;
injuryVectors.kc_fasFlags = kc_fasFlags;
injuryVectors.fancyFilterHandle = fancyFilterHandle;
injuryVectors.fancyFilterParams = fancyFilterParams;
 

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