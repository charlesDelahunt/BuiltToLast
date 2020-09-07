% Set model params and save into a matfile with the name given in argin. Also, output the struct  
% itself. These parameters are used to create connection matrices, and govern how dynamics evolve. 
% They do not create the actual matrices (this is done in 'initializeConnectionMatrices.m')
% FUNCTION VERSION

function [out] = specifyModelParamsSymmetricPnsPins_Fn(paramsFilename, saveFlag)

% the following abbreviations are used:
% n = number of 
% G = glomerulus (so eg nG = number of glomeruli)
% R = response neuron (from antennae): this concept is not used. We
%           use the stim -> glomeruli connections directly
% P = excitatory projection neuron. note sometimes P's stand in for gloms
%       in indexing, since they are 1 to 1
% PI = inhibitory projection neuron
% L = lateral neuron (inhibitory)
% K = kenyon cell (in MB)
% S = stimulus (S)
% fr = fraction
% mu = mean
% sig = standard deviation
% _2_ = synapse connection to, eg P2Kmu = mean synapse strength from PN to KC
% octo = octopamine delivery neuron

% for dynamics equations used in ode solver, see 'dynamics_MBAL_v3.m':

% General structure of synaptic strength matrices:
% rows give the 'from' a synapse
% cols give the 'to'. 
% so M(i,j) is the strength from obj(i) to obj(j)

% Copyright (c) 2018-2020 Charles B. Delahunt.  delahunt@uw.edu
% MIT License

%-----------------------------------------------------

% choose a modelParamsFilename, to save this configuration:
modelParamsFilename = paramsFilename;

% below, 'G' stands for glomerulus. Glomeruli are not explicitly part of the equations, but
% the matrices for LN interconnections, PNs, and RNs are indexed according to the G

nG = 60; 
numPperG = 5;
nP = nG*numPperG; % Pn = # excitatory Pn.

numPIperG = 0;   % param we can vary, so PI:P = 0, 1:5, 2:5, 3:5, 4:5, 5:5 
nPI = nG*numPperG; 

nR = nG; % RNs, one per glom
% for now assume no pheromone gloms. Can add later, along with special Ps, PIs, and Ls
nK = 2000;  % number of kenyon cells (MB). moth has ~4k, so reduced to account for non-olfactory
nS = 2; % number of distinct stimuli (ie odors)

% note that outputs P and PI only affect ML
nE = 1;    % extrinsic neurons in eg beta-lobe of MB, ie read-out neurons

%-------------------------------------------------

% Set param values:

% time constants for the diff eqns:
% CAUTION: an important param: the same param set will give different results if only this one is
% changed. Assume 3*(1/time constant) = 1+ seconds (3 t.c -> 95% decay)
tau = 7;   % >=5 is good for odor (back to equilibrium in 1+ seconds) 
tauR = tau;
tauP = tau;
tauPI = tau;
tauL = tau;
tauK = tau;
tauE = tau;

% sigmoid range parameter for the diff eqns:
C = 10.5;  
cR = C;  
cP = C;
cPI = C;
cL = C;
cK = C;  

% sigmoid slope param for the diff eqns:
% slope of sigmoid at zero = C*slopeParam/4
desiredSlope = 1; % = 1.5 previous to 22 sept;
slopeParam = desiredSlope*4/C;   % a convenience variable. 'desiredSlope' is the one to adjust. 

% sigmoid param to make range 0:C or -C/2:C/2
symmetricAboutZero = 1;  % '0' means: [-inf:inf] -> [0, C], 0 -> C/2. '1' means: [-inf:inf] ->
                         % [-c/2,c/2], 0 -> 0. 
sigFactor = 0.3;   % this value multiplies the mean value of a connection matrix to give the STD.
% currently it applies to all connection matrices, ie it is a global parameter.

% define spontaneous steady-state FRs of RNs:
spontRdistFlag = 2;   % 1 = gaussian, 2 = gamma + base
% % if using gaussian dist:
    % spontRmu = 2.75;  
    % spontRsig = 0.05*spontRmu;  
    % spontRbase = 0; % not used
% if using gamma dist:
spontRbase = 0.3; % to give Rspont = base + gammaDist
spontRmu = 0.1;
spontRsig = 5*spontRmu;  % the shape param

% Rspont gets generated, using a distribution as defined in initializeConnectionMatrices.m
RspontAsInput = 1;  % 1 means it is used as an input to Rinputs (pre-sigmoid) in sde_fn. 
                    % 0 means it is used in the dR eqn, as the steady state that R relaxes to.

% CAUTION RE USE OF 'SIG' PARAMS FOR P, L, PI, R: typically the effect of an input to G, as
% passed through to P, L, PI, and R = mult*effectOnG + sig*NormalDist. mult and sig are different 
% for P, L, PI, R. We need to be careful that sig is scaled properly relative to mult*effectOnG, not
% just relative to mult.

% -----------------------------------------------------
% Stim characteristics:

% distribution of stim->glom number non-zero connections, and strength:
RperSFrMu = 0.25;  %  v good: 0.5; % 0.6; % 0.30;  % for orthogonal stims, try 0.6. for overlapping 
% stims, try 0.3 ave fraction of Rs each S feeds into. The actual # will vary as a binomial distrib
% mu and sig define mean and std of connection matrix entries:
% KEY POINT: R2L, R2P, and R2PI should be highly correlated, since all
% depend on the R2G connectivity strength. For this reason, use an
% intermediate connectivity matrix R2G to give the base levels of
% connectivity. Then derive R2P, R2L, R2PI from L2G.
S2Rmu = 0.25; % controls how strongly a given S will affect G's
S2Rsig = sigFactor*S2Rmu;

makeStimsOrthogonalFlag = 0; % 1 will result in fewer S2R entries per odor, if there are many odors.

% ------------------------------------------------------------------

% R characteristics:
R2Gmu = 4.5;  % controls how strongly R's affect their G's
R2Gsig =  sigFactor*R2Gmu;  
% used to create R2P, using as base R2G:
R2Pmult =  4;    % this multiplies R2G values, to give R2P values
R2Psig = 0.1*R2Pmult;   % variation in R2P values, beyond the conversion from R2G values
R2PImult = R2Pmult;
R2PIsig = R2Psig;
R2Lmult = 1.3*0.75; 
R2Lsig = 0.1*R2Lmult;   

% KEY POINT: L2R, L2L, L2P, and L2PI should be highly correlated, since all
% depend on the L2G connectivity strength. For this reason, use an
% intermediate connectivity matrix L2G to give the base levels of
% connectivity. Then derive L2P, L2R, L2L from L2G.
% distribution of LNs, ie inhib glom-glom synaptic strengths:
L2Gfr =  0.8;  % the fraction of possible LN connections (glom to glom) that are non-zero.   
               % Hong: fr = most or all
L2Gmu = 5;   % 5 = original value                                                                                 
L2Gsig = 0.1*L2Gmu; %                                                                 
killIntraStimLwts = 0;  % flag, 1 says to remove inhib connections between G's affected by same stim
                        % CAUTION: this is deprecated. Do not set = 1.

% sensitivity of G to gaba, ie sens to L, varies substantially with G:
GsensMu = 0.45; 
GsensSig = 0.2*GsensMu;   % Hong-Wilson says PNs (also Gs) sensitivity to gaba varies as N(0.4,0.2).
% But this N( ) is the end-result of many interactions, not a direct application of GsensMu and 
% GsensSig. Note: GsensSig expresses variation in L effect, so if we assume that P, L, and R are all
% similarly gaba-resistent within a given G, set L2Psig etc below = 0.

%---------------------------------------------------

% define effect of LNs on various components:

% used to create L2R, using as base L2G:
L2Rmult =  6/ (nG*L2Gfr);  
L2Rsig =  0.1*L2Rmult;                                                                        

% used to create L2P, using as base L2G:
L2Pmult =  1/ (nG*L2Gfr);    
L2Psig = 0.1*L2Pmult;  % variation in L2P values, beyond the conversion from L2G values            
% used to create L2PI, using as base L2G:
L2PImult =  L2Pmult;
L2PIsig =  0.1*L2PImult;  
% used to create L2L, using as base L2G:
L2Lmult = 2/ (nG*L2Gfr);   
L2Lsig =  0.1*L2Lmult;  

%---------------------------------------------------------------------

% distributions that describe how Ps connect to Ks, ie synapses:
% note that during construction of the synapses, these numbers (#K:#P, #P:#PI) 
% need to become integers by some mechanism (double-check this).
% a) excit PNs. We need # of KCs connected to each PN:
KperPfrMu = 0.2/numPperG; % the mean fraction of KCs a given 'PN' connects to. Note that multiple 
                  % PNs come out of a single glom. 
                  % assuming 5 true PNs per glom, and 2000 KCs, 0.2 in numerator means: 
                  % each glom -> 0.2*2000 = 400 KC, so each true PN -> 80
                  % KCs, so there are 300*80 PN:KC links = 24000 links, so
                  % each KC is linked to 24k/2k = 12 true PNs (see Turner 2008).
                  % The actual # will vary as a binomial distribution
% b) inhib PNs. We need the # of Gs feeding into each PI and # of Ks the PI goes to:
KperPIfrMu = KperPfrMu;  % NOTE! this val is determined by numP, not numPI. ie PIs
                         % will behave the same as Ps, even if there are fewer of them.

% distribution that describes the strength of PN->KC synapses:
% these params apply to both normal and inhib PNs. These are plastic, via
% PN stimulation during octo stim.
P2Kmu = 3;  % less than 4 tends to kill KC outputs (if orth stims and normal sigmoid)
P2Ksig = 0.05*P2Kmu;  % narrower band of values

PI2Kmu = P2Kmu;  % positive since these are synaptic strengths 
PI2Ksig = P2Ksig;

%---------------------------------------------------------------

% KC -> EN connections:

KperEfrMu = 1;  % what fraction of KCs attach to a given EN
K2Emu = 0.1;
K2Esig = 0; 

%------------------------------------------------------------

% the KCs are damped by a global factor depending on total PNs or KCs FR.
% The driving signal (PNs or KCs) is specified in the dynamics functions.
% kGlobalDampFactor = 0 means no damping.
kGlobalDampFactor = 1;  
kGlobalDampSig = 0; % allows variation of damping effect per KC. 
sparsityTarget= 0.05;
octoSparsityTarget = 3*sparsityTarget;

%---------------------------------------------------------------

% distribution of octopamine -> glom strengths (small variation):
% NOTES:
%   1. these values assume octoMag = 1
%   2. if using the multiplicative effect, we multiply the inputs to neuron N by (1 +
%      octo2N) for positive inputs, and (1-octo2N) for negative inputs. So we
%      probably want octo2N to be close to 0. It is always non-neg by
%      construction in 'initializeConnectionMatrices.m'

% In the context of dynamics eqns, octo reduces a neuron's responsivity to negative inputs, but 
% maybe it does this less strongly than it increases the neuron's responsivity to pos inputs:
octoNegDiscount = 0.5; % < 1 means less strong effect on neg inputs. 

% since octo strengths are correlated with a glom, use same method as for L:
octo2Gmu = 1;  
octoSigFactor = 0.3;
octo2Gsig = 0.3*octo2Gmu;  

% Per Jeff, octo affects R, and may also affect P, PI, and L.  
% used to create octo2R, using as base octo2G.
octo2Rmult = 1.75; 
octo2Rsig = 0.1*octo2Rmult; % 0;
% used to create octo2P, using as base octo2G:
octo2Pmult = 0; % 0.02;                                 %  0.8; % 4;  
% comment: octo2P = 0 makes most sense if we only want to stimulate
% response to odor inputs. stimulating P's responsiveness amplifies all signals, inc noise.
octo2Psig = octoSigFactor*octo2Pmult;
% used to create octo2PI, using as base octo2G:
octo2PImult = octo2Pmult;            
octo2PIsig = octo2Psig; 
% used to create octo2L, using as base octo2G:
octo2Lmult = 1.85;  % must be > 0 if Rspont not affected by octo in sde_dynamics (ie octo only 
       % affects R's reaction to odor), since jeff's data shows that some Rspont decrease with octo.                              
octo2Lsig = 0.1*octo2Lmult; % 0;

% end of AL-specific octo param specs.

% distribution of octopamine -> KC strengths (small variation):
octo2Kmu = 0; % 6;  
octo2Ksig = octoSigFactor*octo2Kmu;

octo2Emu = 0;  % for completeness, not used.
octo2Esig = 0; 

%------------------------------------------------------------------------

% Hebbian learning parameters:
hebTauPK = 100; % learning rate for P2K weights. Higher means slower
hebTauPIK = hebTauPK;
hebTauKE = 1000; % learning rate for K2E weights
hebMaxPK = P2Kmu + 5*P2Ksig;
hebMaxPIK = PI2Kmu + 5*PI2Ksig;
hebMaxKE = 1.5*K2Emu + 5*K2Esig;

% powers applied to weight the FRs during heb:
hebPowerP = 1;  % 2;
hebPowerPI = hebPowerP;
hebPowerKin = 1;
hebPowerKout = 1;
hebPowerE = 1; 

dieBackTauKE = 2.5*hebTauKE; % 1/decay rate. Higher means slower decay.  
% dieBackTauKE and hebTauKE want to be in balance.
 
% hebTauPK =  5e3*goal;  % learning rate for P2K weights in MNIST moth template. Higher means slower.  
% Decay: There is no decay for P2K weights in MNIST moth template, nor in original injury template
dieBackTauPK = 250*hebTauPK; %  If > 0, divide this fraction of gains evenly among all nonzero 
                             % weights, and subtract.
dieBackTauPIK = dieBackTauPK;  

%------------------------------------------------------------------------

% add some noise parameters.
% distributions for random variation in P,R,L and K vectors at each step of the simulation:
% these control random variations in the various neurons that are applied at
% each time step as 'epsG' and 'epsK'. This might serve to break up potential oscillations.
noise = 1; % may, mothProject, mothML 
noiseSigFactor = 0.3;  % this FR variation can be spec'ed out by noiseScalingSig below, if wished. 

% following get gamma dist noise:
noiseR = noise;
RnoiseSig = noiseSigFactor*noiseR;  
noiseP = noise; 
PnoiseSig = noiseSigFactor*noiseP;
noiseL = noise;
LnoiseSig = noiseSigFactor*noiseL;
noisePI = noise;
PInoiseSig = noiseSigFactor*noisePI;

noiseK = 0; %  adding noise to KCs muddies repeatability (duh); 
KnoiseSig = noiseSigFactor*noiseK;
noiseE = 0;
EnoiseSig = 0;

% Not currently used:
% another way to inject the noise parameters (effect is the same as using 'noise' and 
% 'noiseSigFactor'):
% Initially, noise is scaled to match mean spont FR.
% But Jeff's data shows that noise is not correlated exactly with spont FR.
% So add a variation to the scaling. This is applied when noise is first added
% during sde_EM_evolution_v1.m:
noiseScalingMu = 1; % 3;
noiseScalingSig = 0;  %  assumed to be the same for all neuron types.

% give the multiplier used when calculating FRs. the FR probs created by
% the sim are multiplied by this, so that the FR added up over a given
% window are higher.
probMultiplier = 1;  % 1 => raw

%-------------------------------------------------------------------------------------------

% save these as a struct. 
% Lines below here should not be modified except when parameters are added or removed.

modelParams.nS = nS;
modelParams.nG = nG;
modelParams.nR = nR;
modelParams.nP = nP;
modelParams.nPI = nPI;
modelParams.nK = nK;
modelParams.nE = nE;
modelParams.numPperG = numPperG;
modelParams.numPIperG = numPIperG;

modelParams.tauR = tauR;
modelParams.tauP = tauP;
modelParams.tauPI = tauPI;
modelParams.tauL = tauL;
modelParams.tauK = tauK;
modelParams.tauE = tauE;

modelParams.cR = cR;
modelParams.cP = cP;
modelParams.cPI = cPI;
modelParams.cL = cL;
modelParams.cK = cK;

modelParams.slopeParam = slopeParam; 
modelParams.symmetricAboutZero = symmetricAboutZero;

modelParams.spontRmu = spontRmu; 
modelParams.spontRsig = spontRsig;
modelParams.RspontAsInput = RspontAsInput;
modelParams.spontRbase = spontRbase;
modelParams.spontRdistFlag = spontRdistFlag;

modelParams.RperSFrMu = RperSFrMu;  
modelParams.S2Rmu = S2Rmu;
modelParams.S2Rsig = S2Rsig;
modelParams.makeStimsOrthogonalFlag = makeStimsOrthogonalFlag;

modelParams.R2Gmu = R2Gmu;
modelParams.R2Gsig = R2Gsig;
modelParams.R2Pmult = R2Pmult;  
modelParams.R2Psig = R2Psig; 
modelParams.R2PImult = R2PImult;
modelParams.R2PIsig = R2PIsig;
modelParams.R2Lmult = R2Lmult; 
modelParams.R2Lsig = R2Lsig;

modelParams.L2Gfr = L2Gfr; 
modelParams.L2Gmu = L2Gmu;
modelParams.L2Gsig = L2Gsig;
modelParams.killIntraStimLwts = killIntraStimLwts; 

modelParams.L2Rmult = L2Rmult;   
modelParams.L2Rsig = L2Rsig;

modelParams.L2Pmult = L2Pmult; 
modelParams.L2Psig = L2Psig;
modelParams.L2PImult = L2PImult;
modelParams.L2PIsig = L2PIsig;
modelParams.L2Lmult = L2Lmult;
modelParams.L2Lsig = L2Lsig;

modelParams.GsensMu = GsensMu;
modelParams.GsensSig = GsensSig;

modelParams.KperPfrMu = KperPfrMu; 
modelParams.KperPIfrMu = KperPIfrMu; 
modelParams.P2Kmu= P2Kmu;
modelParams.P2Ksig = P2Ksig;
modelParams.PI2Kmu = PI2Kmu; 
modelParams.PI2Ksig = PI2Ksig;

modelParams.kGlobalDampFactor = kGlobalDampFactor;
modelParams.kGlobalDampSig = kGlobalDampSig;
modelParams.sparsityTarget = sparsityTarget;
modelParams.octoSparsityTarget = octoSparsityTarget;

modelParams.KperEfrMu = KperEfrMu;
modelParams.K2Emu = K2Emu;
modelParams.K2Esig = K2Esig;

modelParams.octoNegDiscount = octoNegDiscount;
modelParams.octoSigFactor = octoSigFactor;
modelParams.octo2Gmu = octo2Gmu;
modelParams.octo2Gsig = octo2Gsig;
modelParams.octo2Pmult = octo2Pmult;
modelParams.octo2Psig = octo2Psig;
modelParams.octo2PImult = octo2PImult;
modelParams.octo2PIsig = octo2PIsig;
modelParams.octo2Lmult = octo2Lmult;
modelParams.octo2Lsig = octo2Lsig;
modelParams.octo2Rmult = octo2Rmult;
modelParams.octo2Rsig = octo2Rsig;
modelParams.octo2Kmu = octo2Kmu;
modelParams.octo2Ksig = octo2Ksig;
modelParams.octo2Emu = octo2Emu;
modelParams.octo2Esig = octo2Esig;

modelParams.hebTauPK = hebTauPK;
modelParams.hebTauPIK = hebTauPIK;
modelParams.hebTauKE = hebTauKE;
modelParams.hebMaxPK = hebMaxPK;
modelParams.hebMaxPIK = hebMaxPIK;
modelParams.hebMaxKE = hebMaxKE;
modelParams.hebPowerP = hebPowerP; 
modelParams.hebPowerPI = hebPowerPI;
modelParams.hebPowerKin = hebPowerKin;
modelParams.hebPowerKout = hebPowerKout;
modelParams.hebPowerE = hebPowerE; 

modelParams.dieBackTauPK = dieBackTauPK;
modelParams.dieBackTauPIK = dieBackTauPIK;
modelParams.dieBackTauKE = dieBackTauKE;

modelParams.noiseR = noiseR;
modelParams.RnoiseSig = RnoiseSig;
modelParams.noiseP = noiseP;
modelParams.PnoiseSig = PnoiseSig;  
modelParams.noisePI = noisePI;
modelParams.PInoiseSig = PInoiseSig;
modelParams.noiseL = noiseL;
modelParams.LnoiseSig = LnoiseSig;
modelParams.noiseK = noiseK;
modelParams.KnoiseSig = KnoiseSig;
modelParams.noiseE = noiseE;
modelParams.EnoiseSig = EnoiseSig;
modelParams.noiseScalingMu = noiseScalingMu;
modelParams.noiseScalingSig = noiseScalingSig;

modelParams.modelParamsFilename = modelParamsFilename;  % name of this struct, used to load it later
modelParams.probMultiplier = probMultiplier;

if saveFlag 
    save(modelParamsFilename, 'modelParams');
end

out = modelParams; % to be passed to model_MBAL_v5_Fn


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