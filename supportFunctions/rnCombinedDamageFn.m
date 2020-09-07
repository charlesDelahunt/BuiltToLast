function [ rnOut ] = rnCombinedDamageFn( rnIn, lowPassParams )
% rnCombinedDamageFn: applies the pie chart to antenna inputs to the RNs to get an aggregate
% output value (considering all damage types, applied across all 500 inputs to a given RN).
% Inputs:    rn = vector of RN FRs
%            fancyFilterParams has three fields:
%               antennaeMax = maximum FR for antennae inputs (~3)
%               scalingFactor = multiplies antennae inputs to get the values into [0, 10 to 15] for 
%                            Pedro's LP filter to give appropriate behavior. 
%						     Pedro uses scalingFactor = 10 given raw inputs in [0, 1]
%						     Outside of [0 15], Pedro's LP filter equation misbehaves
%							 
%            rn_pie = pieChart for RN: 1 x 3 vector, [ ablationFr, reflectionFr, lowPassFr ]
% NOTE: Ablation-only is a special case, rn_pie = [ n, 0, 0 ] where n = ablation fraction

% Copyright (c) 2018-2020 Charles B. Delahunt.  delahunt@uw.edu
% MIT License

%--------------------------------------------------------------------
 
rn_pie = lowPassParams.rn_pie;

ablateFraction = rn_pie(1);
reflectFraction = rn_pie(2);
lowPassFraction = rn_pie(3);
transmitFraction = 1 - (ablateFraction + reflectFraction + lowPassFraction);  % the undamaged 
                                                                              % fraction

% ie if there is a low-pass component, calculate how it will affect the LP portion of the signal:
if lowPassFraction > 0 
    resultOfLowPassOnFullRnIn = pedroLowPassFilterForRNs(rnIn, lowPassParams.scalingFactor,...
                                                         lowPassParams.maxFR); 
else
    resultOfLowPassOnFullRnIn = 0;
end

rnOut = rnIn.*(1*transmitFraction + 0*ablateFraction + 0.5*reflectFraction) + ...
                                                       resultOfLowPassOnFullRnIn*lowPassFraction;


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

