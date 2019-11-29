function [daughterWavelet, fourierFactor, coneOfInfl] = ...
    continuousWltBasis (waveletScale, fourierFreq, wltParameter)
% This function estimates the daughter of the mother wlt based on the
% parameters stated in the paper entitled A Practical Guide to Wavelet Analysis
% Compo & Torrence (1998)

% INPUTS
%     1. motherWavelet: the chosen mother wavelet,
%     2. waveletScale: the wavelet scale,
%     3. fourierFreq: a vector containing the Fourier Frequencies,  
%     4. wltParameter: the order of the derivative of the Gaussian (DOG) or
%        the central frequency of the Morlet wavelet
% 
% OUTPUTS
%     1. normalizationFactor: the factor estimated by appliying equation (6),
%     2. waveletScale: the scale of the wavelet transform estimated by applying eq. (9) and (10),
%     3. daughterWavelet: the conjugate of the Fourier transform of the mother Wavelet
%
% see also the Bedforms-ATM (Gutierrez et al., 2013).
% Author: Dr. Ronald R. Gutierrez
%%    
centralFreq = wltParameter;
nFourierFreq=length(fourierFreq);
normalizFactor = sqrt(waveletScale*fourierFreq(2))*(pi^(-0.25))*sqrt(nFourierFreq);
expArgument = -(waveletScale.*fourierFreq-centralFreq).^2/2.*(fourierFreq>0);
daughterWavelet = normalizFactor*exp(expArgument);
daughterMorlet = daughterWavelet.*(fourierFreq>0);
fourierFMorlet = (4*pi)/(centralFreq+sqrt(2+centralFreq^2));
coneMorlet = fourierFMorlet/sqrt(2); 

daughterWavelet = daughterMorlet;
fourierFactor = fourierFMorlet;
coneOfInfl = coneMorlet;

