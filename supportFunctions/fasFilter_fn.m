function [ xOut ] = fasFilter_fn( xIn, fasFlags, fancyHandle, lowPassParams )
% fasFilter(neuronFRs, fasFlags) applies injury types to neuron instantaneous FRs, after inputs
% have been sigmoided and before FRs have been updated.
% 1. xIn = vector output of the sigmoid fn that squashes summed inputs
% 2. fasFlags = vector of flags that say what type of FAS injury these neurons have suffered:
%       0 = no injury:      xOut = xIn
%       1 = ablation:       xOut = 0
%       2 = reflection:     xOut = 0.5*xIn
%       3 = fancy filter:   xOut = @fancyHandle(xIn, params) 
% 3. fancyHandle = function handle for a filtering fn (optional arg)
% 4. lowPassParams = struct with various fields for use with fancyHandle. (optional arg)

% Copyright (c) 2018 - 2020, Charles B. Delahunt.  delahunt@uw.edu
% MIT License

%-------------------------------------------------

xOut = xIn;
xOut(fasFlags == 1) = 0;
xOut(fasFlags == 2) = 0.5*xIn(fasFlags == 2);

% if some xIns require fancy filtering, check if there is a handle:
assert( sum(fasFlags == 3) == 0 || (sum(fasFlags == 3) > 0 &&  ~isempty(fancyHandle) ),...
    'fasFilter_fn.m: fasFlags == 3 but there is no fancy filtering fn handle');

% Apply the low pass filter to the low-passed neurons only, if there are any:
fancyFilterFn = fancyHandle; 
in = xIn(fasFlags == 3);
if ~isempty(in)
    out = fancyFilterFn( in, lowPassParams );
    xOut(fasFlags == 3) = out;
end

end


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
