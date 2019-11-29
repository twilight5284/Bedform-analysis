function [peaksCross,xPeakspos,yPeakspos] = getWltpeaks(globalWs,globalSignif,...
                                          fourierPeriod)                                      
% This function calculates the value and coordinates of the peaks
% in the global wavelet significance of the wavelet transform
% The peaks that are greater than the global significance

% INPUTS
%     1. globalWs: The global wavelet significance,
%     2. globalSignif: the global significance,
%     3. fourierPeriod: a vector containing the Fourier Periodes 
%        for each point in globalWs,  
% This values are outputs of the statisticsWlt function
% 
% OUTPUTS
%     1. peaksCross: value of the peak,
%     2. xPeakspos: x coordinate to print the peak,
%     3. yPeakspos: y coordinate to print the peak 
%
% see also the Bedforms-ATM (Gutierrez et al., 2013).
% Author: Dr. Ronald R. Gutierrez
%%

lengthWS = length(globalWs);
xposFM = 0.0;
PrintFPeak = 0;
peaksCross =[];
xPeakspos =[];
yPeakspos =[];

for rr=1:lengthWS
    if (globalWs(rr)>globalSignif(rr))
        if (globalWs(rr)>xposFM) 
            PrintFPeak = 1;
            xposFM = globalWs(rr);
            yposFM = log2(fourierPeriod(rr));
        end
    else
        if (PrintFPeak==1)
            peaksCross (end+1) =2^yposFM;
            xPeakspos(end+1) =  xposFM;
            yPeakspos(end+1) =  yposFM;
            PrintFPeak=0;
        end
    end
end
