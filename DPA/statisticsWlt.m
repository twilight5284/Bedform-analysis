function [wltSignif, globalSignif, globalWs,lag1,scaleAveSignif,scaleAvg] =...
    statisticsWlt...
    (wltParameter, signalData, sampleFreq, deltaFreq,...
    signifLevel,wltPower)
%% Function that carries out the statistical analysis of the wavelets
% see also the Bedforms-ATM (Gutierrez et al., 2013).
% Author: Dr. Ronald R. Gutierrez

%% scale
lowerScale = 0.001;
upperScale = 10000;

%% General data
% Estimate the lag1 coefficient for the Red Noise Power Spectrum        
% for 10 Lag values (-10 to 10) according to matlab reference it is good
% for getting the lag of the white noise

[c,lags] = xcorr(signalData,10,'coeff');
indexLag1 = find(lags==1);
lag1 = c(indexLag1);
lag1 = 0.99;

nSignalData    = length(signalData); 
signalVariance = std(signalData)^2;
s0 = 2*sampleFreq ;
largestScale = fix((log2(nSignalData*sampleFreq/s0)))/deltaFreq;
wltScale = s0*2.^(deltaFreq*(0:largestScale)); % puede pasarse como parametro
% Get the variance of the wavelet function

%% Significance of the power wavelet
% Significance levels: (variance=1 for the normalized SST)
% [signalSignif, redNoiseSpectr] = significance(1.0,dt,scale,0,lag1,SIGLVL,-1,mother);
[signalSignif] = significance...
    (wltParameter, 1, nSignalData, sampleFreq, deltaFreq,...
    signifLevel,0,lag1,wltScale,lowerScale, upperScale);
wltSignif = (signalSignif')*(ones(1,nSignalData));
wltSignif = wltPower./wltSignif;

%% Global spectrum and significance
globalWs = signalVariance^1*(sum(wltPower')/nSignalData);
% globalWsfps = (sum(wltPower)/nSignalData);
% globalWsfps = globalWsfps/max(globalWsfps);
% globalFftps = (fourierTSignal');         % time-average over all times 

% globalFftps = globalFftps/max(ftSignal);
% dof = nSignalData-wltScale; it is inside significance function (RONALD)***

globalSignif = significance...
    (wltParameter, signalVariance, nSignalData, sampleFreq,...
    deltaFreq, signifLevel,1,lag1,wltScale,lowerScale, upperScale);

%% Scale Average Significance
avg = find((wltScale >= lowerScale) & (wltScale < upperScale));

  Cdelta = 2;        
  if (wltParameter==6), Cdelta = 0.776; end
   

scaleAvg = (wltScale')*(ones(1,nSignalData));
scaleAvg = wltPower./scaleAvg;
scaleAvg = signalVariance*sampleFreq*sampleFreq/Cdelta*sum(scaleAvg(avg,:));

scaleAveSignif = significance...
    (wltParameter, signalVariance, nSignalData, sampleFreq,...
    deltaFreq, signifLevel,2,lag1,wltScale,lowerScale, upperScale);


function [signalSignif] = significance...
    (wltParameter, signalVariance, nSignalData, sampleFreq, deltaFreq,...
    signifLevel,analysisMethod,lag1,wltScale,lowerScale, upperScale)
%     (wltParameter, signalData, sampleFreq, deltaFreq,...
%     signifLevel,analysisMethod,lag1)
%% Estimation of the Significance
% Set parameters of the mother wavlets
  degreOfFrdm = 2;        
  fourierFactor = (4*pi)/(wltParameter+sqrt(2+wltParameter^2)); % Scale-->Fourier [Sec.3h]
  empFactors = [2,-1,-1,-1];   
  if (wltParameter==6), empFactors(2:4) = [0.776,2.32,0.60]; end
% % Characteristics of the Signal
% nSignalData = length(signalData);    % estimate the length of the signal

% % Scale of the Wavelet 
% s0 = 2*sampleFreq;
% largestScale = fix((log2(nSignalData*sampleFreq/s0)))/deltaFreq;
% wltScale = s0*2.^(deltaFreq*(0:largestScale));

% Calculator
nWltScale = length(wltScale)-1;
% signalVariance = std(signalData)^2;
dj = log(wltScale(2)/wltScale(1))/log(2);
period = wltScale.*fourierFactor;
freqIndex = sampleFreq./period;   % normalized frequency
minDOF = empFactors(1);           % Degrees of freedom with no smoothing
recFactor = empFactors(2);        % reconstruction factor
gammaFactor = empFactors(3);      % time-decorrelation factor
dj0 =  empFactors(4);             % scale-decorrelation factor
redNoiseSpectr = (1-lag1^2)./(1-2*lag1*cos(freqIndex*2*pi)+lag1^2);  % from Eq(16)
redNoiseSpectr = signalVariance*redNoiseSpectr;  % include time-series variance
signalSignif = redNoiseSpectr;

% if (degreOfFrdm==1), degreOfFrdm = minDOF; end

switch(analysisMethod)
    case 0,                      % Just do a regular chi-square test Eqn (18)
        degreOfFrdm = minDOF;
        X = chi2inv(signifLevel,degreOfFrdm);
        chiSq = X/degreOfFrdm;
        signalSignif = redNoiseSpectr*chiSq;        
    case 1,                      % "time-average" test Eqn (23)
        degreOfFrdm = nSignalData-wltScale; % rearanged from the original code
        if (length(degreOfFrdm)==1),           
        degreOfFrdm = zeros(1,nWltScale+1)+degreOfFrdm; end
        trunc = find(degreOfFrdm<1);
        degreOfFrdm(trunc) = ones(size(trunc));
        degreOfFrdm = minDOF*sqrt(1+(degreOfFrdm*sampleFreq/gammaFactor./wltScale).^2); %Eq (23)
        trunc = find(degreOfFrdm<minDOF);
        degreOfFrdm(trunc) = minDOF*ones(size(trunc));
        for iterm = 1:nWltScale+1
            Y = chi2inv(signifLevel,degreOfFrdm(iterm));
            chiSq = Y/degreOfFrdm(iterm);
            signalSignif(iterm) = redNoiseSpectr(iterm)*chiSq;
        end       
    case 2,                     % do a "scale-average" test, Eqns (25)-(28)
        s1 = lowerScale;
        s2 = upperScale;
        avrg = find((wltScale>=s1) & (wltScale<=s2));
        navg = length(avrg);
        if (navg==0)
            error(['No valid scales between ',num2str(s1),' and ',num2str(s2)])
        end
%         Savg = (sum(1./WltScale(avrg)))^(-1);
        Savg = 1./sum(1./wltScale(avrg));
        Smid = exp((log(s1)+log(s2))/2);
        redNoiseSpectr = Savg*sum(redNoiseSpectr(avrg)./wltScale(avrg));
        degreOfFrdm = (minDOF*navg*Savg/Smid)*sqrt(1+(navg*dj/dj0)^2);
        Z = chi2inv(signifLevel,degreOfFrdm);
        chiSq = Z/degreOfFrdm;
        signalSignif = (dj*sampleFreq/recFactor/Savg)*redNoiseSpectr*chiSq;
    otherwise disp('No available type of statistical analysis');
end


