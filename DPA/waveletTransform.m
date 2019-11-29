function [wltPower, fourierPeriod, coneOfInfluence] = ...
    waveletTransform ( wltParameter, signalData, sampleFreq, deltaFreq)
%% wavelet transform
% see also the Bedforms-ATM (Gutierrez et al., 2013).
% Author: Dr. Ronald R. Gutierrez

%%

% Signal data
nsignalData = length(signalData);    % estimate the length of the signal

s0           = 2*sampleFreq;
largestScale = fix((log2(nsignalData*sampleFreq/s0)))/deltaFreq;
wltScale     = s0*2.^(deltaFreq*(0:largestScale));

% Normalizing with respect to the mean
padedSignal = signalData-mean(signalData);

% Padding with zeros
nSignalBase2    = nextpow2(nsignalData);
padedSignal     = [padedSignal, zeros(1,2^(nSignalBase2)-nsignalData)];  % Signal padded with zeros
nPaded          = length(padedSignal); % Length of the signal padded with zeros 

% Estimate angular frequency - equation (5)

% Compo's angular frequency
fourierFreq = [1:fix(nPaded/2)];
fourierFreq = fourierFreq.*((2.*pi)/(nPaded*sampleFreq));
fourierFreq = [0., fourierFreq, -fourierFreq(fix((nPaded-1)/2):-1:1)];

% Fourier transform of the paded signal
ftPaded_Signal = fft(padedSignal);

%% Wavelet Transform

% Set the wlt matrix and then make it complex
wlt0 = zeros(largestScale+1,nPaded);                 
wlt0 = wlt0+1i*wlt0; 

% Estimate the Parseval product
for iScale = 1:largestScale+1
    [daughterWavelet, fourierFactor, coneOfInfl] = ...
    continuousWltBasis ( wltScale(iScale), fourierFreq, wltParameter);
    wlt0(iScale,:) = ifft(ftPaded_Signal.*daughterWavelet); 
end

% timeEnergy = ifft(dautherWavelet);

% Take out the zeros of the paded wavelet transform
wltTransform = wlt0(:,1:nsignalData);                 

% Power of the wavelet
wltPower =(abs(wltTransform)).^2;

% Fourier Period of each scale
fourierPeriod = fourierFactor.*wltScale;

% Global Cone of influence
coneOfInfluence = coneOfInfl*sampleFreq*[1E-5,1:((nsignalData+1)/2-1),...
fliplr((1:(nsignalData/2-1))),1E-5];  % COI [Sec.3g]




