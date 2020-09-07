% Function that generates the various connection matrices, given a modelParams struct.

function [modelParams] = initializeConnectionMatricesSymmetricPnsPins(modelParams)

% input: 'modelParams', a struct
% output: 'params', a struct that includes connection matrices and all info
%         necessary to FR evolution and plotting.
%
% step 1: unpack input struct 'modelParams'
% step 2: build the matrices
% step 3: pack the matrices into a struct 'params' for output
% These steps are kept separate for clarity of step 2.

% Copyright (c) 2018-2020 Charles B. Delahunt.  delahunt@uw.edu
% MIT License

%--------------------------------------------------------------------

%% step 1: unpack modelParams:
% no editing necessary in this section.
nG = modelParams.nG;
nP = modelParams.nP;
nR = modelParams.nR;
nPI = modelParams.nPI;
nK = modelParams.nK;
nS = modelParams.nS;
nE = modelParams.nE;
numPperG = modelParams.numPperG;
numPIperG = modelParams.numPIperG;

tauR = modelParams.tauR;
tauP = modelParams.tauP;
tauPI = modelParams.tauPI;
tauL = modelParams.tauL;
tauK = modelParams.tauK;

cR = modelParams.cR;
cP = modelParams.cP;
cPI = modelParams.cPI;
cL = modelParams.cL;
cK = modelParams.cK;

if isfield(modelParams,'spontRdistFlag') 
    spontRdistFlag = modelParams.spontRdistFlag;
else
    spontRdistFlag = 1; % gaussian 
end
spontRmu = modelParams.spontRmu;
spontRsig = modelParams.spontRsig; 
% for gamma dist only:
if isfield(modelParams,'spontRbase') 
    spontRbase = modelParams.spontRbase;
end

RperSFrMu = modelParams.RperSFrMu; 
S2Rmu = modelParams.S2Rmu;
S2Rsig = modelParams.S2Rsig;
R2Gmu = modelParams.R2Gmu;
R2Gsig = modelParams.R2Gsig;
R2Pmult = modelParams.R2Pmult;  
R2Psig = modelParams.R2Psig;
R2PImult = modelParams.R2PImult;
R2PIsig = modelParams.R2PIsig;
R2Lmult = modelParams.R2Lmult; 
R2Lsig = modelParams.R2Lsig;

L2Gfr = modelParams.L2Gfr;
L2Gmu = modelParams.L2Gmu;
L2Gsig = modelParams.L2Gsig;
killIntraStimLwts = modelParams.killIntraStimLwts; 

L2Rmult = modelParams.L2Rmult;  
L2Rsig = modelParams.L2Rsig;

L2Pmult = modelParams.L2Pmult; 
L2Psig = modelParams.L2Psig;
L2PImult = modelParams.L2PImult;
L2PIsig = modelParams.L2PIsig;
L2Lmult = modelParams.L2Lmult;
L2Lsig = modelParams.L2Lsig;

GsensMu = modelParams.GsensMu;
GsensSig = modelParams.GsensSig;

KperEfrMu = modelParams.KperEfrMu;
K2Emu = modelParams.K2Emu;
K2Esig = modelParams.K2Esig;

octo2Gmu = modelParams.octo2Gmu;
octo2Gsig = modelParams.octo2Gsig;
octo2Pmult = modelParams.octo2Pmult;
octo2Psig = modelParams.octo2Psig;
octo2PImult = modelParams.octo2PImult;
octo2PIsig = modelParams.octo2PIsig;
octo2Lmult = modelParams.octo2Lmult;
octo2Lsig = modelParams.octo2Lsig;
octo2Rmult = modelParams.octo2Rmult;
octo2Rsig = modelParams.octo2Rsig;
octo2Kmu = modelParams.octo2Kmu;
octo2Ksig = modelParams.octo2Ksig;
octo2Emu = modelParams.octo2Emu;
octo2Esig = modelParams.octo2Esig;
    
noiseR = modelParams.noiseR;
RnoiseSig = modelParams.RnoiseSig;
postInjuryNoiseRnRatio = modelParams.postInjuryNoiseRnRatio;   % since RN neural noise is increased 
                                                               % by RN damage.
noiseP = modelParams.noiseP; 
PnoiseSig = modelParams.PnoiseSig;
noisePI = modelParams.noisePI;
PInoiseSig = modelParams.PInoiseSig;
noiseL = modelParams.noiseL;
LnoiseSig = modelParams.LnoiseSig;
noiseK = modelParams.noiseK;
KnoiseSig = modelParams.KnoiseSig;
noiseE = modelParams.noiseE;
EnoiseSig = modelParams.EnoiseSig;
noiseScalingMu = modelParams.noiseScalingMu;
noiseScalingSig = modelParams.noiseScalingSig; 

KperPfrMu = modelParams.KperPfrMu; 
KperPIfrMu = modelParams.KperPIfrMu; 
P2Kmu = modelParams.P2Kmu;
P2Ksig = modelParams.P2Ksig;
PI2Kmu = modelParams.PI2Kmu; 
PI2Ksig = modelParams.PI2Ksig;
kGlobalDampFactor = modelParams.kGlobalDampFactor;
kGlobalDampSig = modelParams.kGlobalDampSig;
hebMaxPK = modelParams.hebMaxPK;
hebMaxPIK = modelParams.hebMaxPIK;
hebMaxKE = modelParams.hebMaxKE;

makeStimsOrthogonalFlag = modelParams.makeStimsOrthogonalFlag;
modelParamsFilename = modelParams.modelParamsFilename;

%-------------------------------------------------------------------------

%% Step 2: generate connection matrices
  
S2Rbinary = rand(nR, nS) < RperSFrMu; % 1s and 0s.
if makeStimsOrthogonalFlag 
    % remove any overlap in the active odors, by keeping only one non-zero entry in each row:
    b = S2Rbinary;
    for i = 1:nR 
        row = b(i,:);
        if sum(row) > 1 
            c = find(row == 1);
            t = ceil( rand(1,1)*length(c) );   % pick one index to be non-zero
            b(i,:) = 0;
            b( i, c(t) ) = 1;
        end
    end
    S2Rbinary = b;
end

% modify to give varying strengths of S to Rs:
S2R = ( S2Rmu*S2Rbinary + S2Rsig*randn(size(S2Rbinary)) ).*S2Rbinary; % the last term ensures 
                                                                      % 0s stay 0s
S2R = max(0, S2R); % to prevent any negative weights

% spontaneous FRs for Rs:
% % gaussian version:
    % Rspont = spontRmu*ones(nG, 1) + spontRsig*randn(nG, 1);
    % Rspont = max(0, Rspont);
% gamma version:
    a = spontRmu/spontRsig;
    b = spontRmu/a;  % = spontRsig
    g = makedist( 'gamma', 'a', a, 'b', b );    
    Rspont = spontRbase + random( g, [nR, 1] );
    
 % NOTE RE EFFECT OF RN INJURY ON RN SPONT:
 % We assume that spontaneous activity of RNs 'Rspont' is not caused by ORNs
 % spontaneous firing , but noiseR is controlled by ORN firing. So the destruction of ORN -> RN 
 % axons will increase the amount of noise as controlled by noiseR and noiseRsig, but will NOT
 % affect Rspont. If we assumed the opposite, then injury to ORN -> RN axons would reduce spontRmu 
 % and increase spontRsig, as well as increasing noiseR.    

% R2G connection vector. nG x 1 col vector:
R2G  = max( 0, R2Gmu*ones(nG, 1) + R2Gsig*randn(nG, 1) ); % col vector, each entry is strength of an
                                           % R in its G the last term prevents negative R2G effects.
% now make R2P, etc, all are cols nG x 1:

% R2P:
% There are multiple PNs in each glom. Make R affect all these the same by
% setting groups of rows of R2P equal:
temp = zeros(nP,nG);
mask = zeros(nP,nG);    % to record which vals are non-zero
for i = 1:nG
    temp((i-1)*numPperG + 1:i*numPperG, i ) = R2G(i);
    mask((i-1)*numPperG + 1:i*numPperG, i ) = 1;
end

R2P = ( R2Pmult + R2Psig*randn(nP, nG) ).*temp;
R2P = R2P.*mask;
% R2P is nP x nG, with groups of correlated rows. Each row has exactly one non-zero entry.
% The PNs within a glom still vary by R2Psig.

% R2L

R2L = ( R2Lmult + R2Lsig*randn(nG, 1) ).*R2G;

% R2PI:
% There are multiple PINs in each glom. Make R affect all these the same by setting groups of rows 
% of R2PI equal:
temp = zeros(nPI,nG);
mask = zeros(nPI,nG);    % to record which vals are non-zero
for i = 1:nG
    temp((i-1)*numPIperG + 1:i*numPIperG, i ) = R2G(i);
    mask((i-1)*numPIperG + 1:i*numPIperG, i ) = 1;
end
R2PI = ( R2PImult + R2PIsig*randn(nPI, nG) ).*temp;
R2PI = R2PI.*mask;

% Construct L2G = nG x nG matrix of lateral neurons. This is a precursor to L2P etc
L2G = max( 0, L2Gmu + L2Gsig*randn(nG) ); % kill any vals < 0
% set diagonal = 0:
L2G = L2G - diag(diag(L2G));

% if wished, set any inhib connections between stim receivers to 0:
if killIntraStimLwts 
    for i = 1:nS  
        inds = find(S2R(:, i) > 0); %i'th col is i'th stim
        L2G(inds,inds) = 0;
    end
end

% Check if enough of these values = 0:
numZero = sum(L2G(:) == 0) - nG;  % ignore the diagonal zeroes
numToKill = floor( (1-L2Gfr)*(nG^2 - nG) - numZero ); 
if numToKill > 0  % case: we need to set more vals to 0 to satisfy frLN constraint:
    L2G = L2G(:);
    randList = rand(size(L2G) ) < numToKill/(nG^2 - nG - numZero);
    L2G (L2G > 0 & randList == 1) = 0;
end
L2G = reshape(L2G,[nG,nG]);
% Structure of L2G:
% L2G(i,j) = the synaptic LN weight going to G(i) from G(j),
% ie the row gives the 'destination glom', the col gives the 'source glom'

% gloms vary widely in their sensitivity to gaba (Hong, Wilson 2014). 
% multiply the L2* vectors by Gsens + GsensSig:
gabaSens = GsensMu + GsensSig*randn(nG,1);
L2GgabaSens = L2G.*repmat(gabaSens,[1,nG]);   % ie each row is multiplied by a different value, 
                                              % since each row represents a destination glom
% make versions of L2*gabaSens for PNs and PINs, keeping their shared-glom structure: 
% P:
temp = zeros( nP, nG);
for i = 1:nG
    temp((i-1)*numPperG+1:i*numPperG,: ) = repmat(L2GgabaSens(i,:), [numPperG, 1 ] );
end
L2PgabaSens = temp;  % this gaba vector is nP x nL, but with groups of rows equal.
% PI:
temp = zeros( nPI, nG);
for i = 1:nG
    temp((i-1)*numPIperG+1:i*numPIperG,: ) = repmat(L2GgabaSens(i,:), [numPIperG, 1 ] );
end
L2PIgabaSens = temp;  % this gaba vector is nP x nL, but with groups of rows equal.

% now generate all the L2etc matrices

L2R = max( 0, ( L2Rmult + L2Rsig*randn(nG) ).*L2GgabaSens );  % the last term keeps 0 entries = 0 
L2P = max( 0, ( L2Pmult + L2Psig*randn(nP, nG) ).*L2PgabaSens ); 
L2L = max( 0, ( L2Lmult + L2Lsig*randn(nG) ).*L2GgabaSens );
L2PI = max( 0, ( L2Lmult + L2PIsig*randn(nPI, nG) ).*L2PIgabaSens ); % Masked by G2PI later

%% AL -> KC connections:
% Ps (excit):
% no effort is made to link the weights of PNs from the same glom.
P2KconnMatrix = rand(nK, nP) < KperPfrMu; % each col is a P, and a fraction of the entries will = 1.
        % different cols (PNs) will have different numbers of 1's (~binomial dist).
P2K = max (0, P2Kmu + P2Ksig*randn(nK, nP) ); % all >= 0
P2K = P2K.*P2KconnMatrix; 
% cap P2K values at hebMaxP2K, so that hebbian training never decreases wts:
P2K = min(P2K, hebMaxPK);
% PKwt maps from the Ps to the Ks. Given firing rates P, PKwt gives the
% effect on the various Ks
% It is nK x nP with entries >= 0.

% PIs (inhib): 
PI2KconnMatrix = rand(nK, nPI) < KperPIfrMu; % same method as for P2K
PI2K = max (0, PI2Kmu + PI2Ksig*randn(nK, nPI) ); % all >= 0
PI2K = PI2K.*PI2KconnMatrix; 
% cap P2K values at hebMaxP2K, so that hebbian training never decreases wts:
PI2K = min(PI2K, hebMaxPK);

% K2E (excit):
K2EconnMatrix = rand(nE, nK) < KperEfrMu; % each col is a K, and a fraction of the entries will = 1.
        % different cols (KCs) will have different numbers of 1's (~binomial dist).
K2E = max (0, K2Emu + K2Esig*randn(nE, nK) ); % all >= 0
K2E = K2E.*K2EconnMatrix;  
K2E = min(K2E, hebMaxKE);
% K2E maps from the KCs to the ENs. Given firing rates KC, K2E gives the effect on the various ENs.
% It is nE x nK with entries >= 0.

%% octopamine to Gs and to Ks:
octo2G = max( 0, octo2Gmu + octo2Gsig*randn(nG, 1) );  % intermediate step 
octo2K = max( 0, octo2Kmu + octo2Ksig*randn(nK, 1) );
% each of these is a col vector with entries >= 0

% octo2P requires fussing because there are multiple PNs per glom:
octo2P = zeros(nP,1); 
mask = zeros(nP,1);    % to record which vals are non-zero
for i = 1:nG
    octo2P((i-1)*numPperG + 1:i*numPperG ) = octo2G(i);
end
octo2P = max(0, octo2Pmult*octo2P + octo2Psig*randn(nP, 1) ); % effect of octo on P, includes gaussian variation from P to P

% octo2PI also requires fussing because there are multiple PINs per glom:
octo2PI = zeros(nPI,1); 
mask = zeros(nPI,1);    % to record which vals are non-zero
for i = 1:nG
    octo2PI((i-1)*numPIperG + 1:i*numPIperG ) = octo2G(i);
end
octo2PI = max(0, octo2PImult*octo2PI + octo2PIsig*randn(nPI, 1) ); % effect of octo on P, includes gaussian variation from P to P

octo2L = max(0, octo2Lmult*octo2G + octo2Lsig*randn(nG, 1) );
octo2R = max(0, octo2Rmult*octo2G + octo2Rsig*randn(nG, 1) );

octo2E = max(0, octo2Emu + octo2Esig*randn(nE, 1) );                                          

%% Noise:

% % each neuron has slightly different noise levels for sde use. Define noise vectors for each type:

% % Gaussian versions:
% noiseRvec = epsRsig + RnoiseSig*randn(nR, 1);
% noiseRvec = max(0, noiseRvec);   % remove negative noise entries
% noisePvec = epsPsig + PnoiseSig*randn(nP, 1);
% noisePvec = max(0, noisePvec);
% noiseLvec = epsLsig + LnoiseSig*randn(nG, 1);
% noiseLvec = max(0, noiseLvec);
% noisePIvec = noisePI + PInoiseSig*randn(nPI, 1);
% noisePIvec = max(0, noisePIvec);

% % gamma versions:
if noiseR > 0
    a = noiseR/RnoiseSig;
    b = noiseR/a;
    g = makedist( 'gamma', 'a', a, 'b', b );
    % in order to recalculate noiseR after injury (which affects param b),
    % we need to reuse the same random draws from the gamma dist:
    seedVal = randi(1e6, 1);
    rng(seedVal);
    noiseRvec = random(g,[nR,1]);
    noiseRvec(noiseRvec > 15) = 0;   % experiment to see if just outlier noise vals boost KC noise
    % now post-injury:
    % 'a' does not change because the multiplier due to injury cancels in top and bottom.
    b = RnoiseSig*postInjuryNoiseRnRatio;
    g = makedist( 'gamma', 'a', a, 'b', b );
    rng(seedVal);
    postInjuryNoiseRvec = random(g,[nR,1]);
    postInjuryNoiseRvec(postInjuryNoiseRvec > 15) = 0;   % experiment to see if just outlier noise 
                                                         % vals boost KC noise
    
    rng('shuffle'); % reset random seed to totally random
else
    noiseRvec = zeros(nR,1);
    postInjuryNoiseRvec = zeros(nR,1);
end

if noiseP > 0
    a = noiseP/PnoiseSig;
    b = noiseP/a;
    g = makedist( 'gamma', 'a', a, 'b', b );
    noisePvec = random(g,[nP,1]);
    noisePvec(noisePvec > 15) = 0;   % experiment to see if outlier noise vals boost KC noise
else
    noisePvec = zeros(nP,1);
end 
if noisePI > 0
    a = noisePI/PInoiseSig;
    b = noisePI/a;
    g = makedist( 'gamma', 'a', a, 'b', b );
    noisePIvec = random(g,[nPI,1]);
    noisePIvec(noisePIvec > 15) = 0;   % experiment to see if outlier noise vals boost KC noise
else
    noisePIvec = zeros(nPI,1);
end
if noiseL > 0
    a = noiseL/LnoiseSig;
    b = noiseL/a;
    g = makedist( 'gamma', 'a', a, 'b', b );
    noiseLvec = random(g,[nG,1]);
else
    noiseLvec = zeros(nG,1);
end
% KCs, ENs: = 0
noiseKvec = noiseK + KnoiseSig*randn(nK, 1);    % = 0 anyway
noiseKvec = max(0, noiseKvec);
noiseEvec = noiseE + EnoiseSig*randn(nE, 1);    % = 0 anyway
noiseEvec = max(0, noiseEvec );

% in general, noise is scaled to neuron's mean spont FR. To vary from this,
% make vectors that modify the exact correlation:
noiseScalingR = max( 0, noiseScalingMu + noiseScalingSig*randn( nR, 1 ) );
noiseScalingP = 1; % max( 0, 1); %noiseScalingMu + noiseScalingSig*randn( nP, 1 ) ); 
noiseScalingPI = 1; %  max( 0, 1); %noiseScalingMu + noiseScalingSig*randn( nPI, 1 ) ); 
noiseScalingL = 1; %  max( 0, 1); % noiseScalingMu + noiseScalingSig*randn( nG, 1 ) ); 
noiseScalingK = 1; %  max( 0,  1); %noiseScalingMu + noiseScalingSig*randn( nK, 1 ) );
noiseScalingE = 1;
kGlobalDampVec = kGlobalDampFactor + kGlobalDampSig*randn(nK,1);  % each KC is affected a bit 
                                                                  % differently by LH inhibition.
%-------------------------------------------------------------------------------

%% append these matrices to 'modelParams' struct:
% no editing necessary in this section

modelParams.modelParamsFilename = modelParamsFilename;

modelParams.S2R = S2R;
modelParams.R2P = R2P;
modelParams.R2PI = R2PI;
modelParams.R2L = R2L;
modelParams.octo2R = octo2R;
modelParams.octo2P = octo2P;
modelParams.octo2PI = octo2PI;
modelParams.octo2L = octo2L;
modelParams.octo2K = octo2K;
modelParams.octo2E = octo2E;
modelParams.L2P = L2P;
modelParams.L2L = L2L;
modelParams.L2PI = L2PI;
modelParams.L2R = L2R;
modelParams.P2K = P2K;
modelParams.PI2K = PI2K;
modelParams.K2E = K2E;
modelParams.Rspont = Rspont;  % col vector

modelParams.noiseRvec = noiseRvec;
modelParams.postInjuryNoiseRvec = postInjuryNoiseRvec;
modelParams.noisePvec = noisePvec;
modelParams.noisePIvec = noisePIvec;
modelParams.noiseLvec = noiseLvec;
modelParams.noiseKvec = noiseKvec;
modelParams.noiseEvec = noiseEvec;
modelParams.kGlobalDampVec = kGlobalDampVec;
modelParams.noiseScalingR = noiseScalingR;
modelParams.noiseScalingP = noiseScalingP;
modelParams.noiseScalingPI = noiseScalingPI;
modelParams.noiseScalingL = noiseScalingL;
modelParams.noiseScalingK = noiseScalingK;
modelParams.noiseScalingE = noiseScalingE;
 

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