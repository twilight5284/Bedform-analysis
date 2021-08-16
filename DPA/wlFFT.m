function [dirArrayM,lambdaArrayM,phiM,stdDir,lambdaM,stdLambda] = wlFFT(x,y,z,L,resolution)
%% Regional dominant pattern calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The dominant regional dune orientation and wavelength of a surface are 
% calculated by 2D DFT 
% 
% INPUTS:
%     x,y,z - grided matrixes of input surface;
%     L - wavelength of interest;
%
% OUTPUTS:
%     dirArrayM, lambdaArrayM - all picked results, orientation and wavelength;                        
%     phiM, lambdaM - the average of orientation and wavelength;
%     stdDir, stdLambda - the standard deviation of orientation and wavelength;
% 
% Several technical improvements have been applied to this calculation, based on
% QBA (Quantitative Bedform Analysis, Cazenave et al., 2013).
%
% See also README.txt
%
% Author:  Li Wang
% Email:   twilight528400@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the threshold for circular filtering based on the wavelength of interest
filtThresh = 3*L;
if L > 100
    filtThresh = 2*L;
end    

% the resolution of meshs
xInc = abs(x(1,2))-abs(x(1,1));
yInc = abs(y(2,1))-abs(y(1,1));

% Grids number of the long edges of the rectangular matrix
NX = size(z,2);
NY = size(z,1);

if NX > NY
    Nxy = NX;
else
    Nxy = NY;
end    

%  fimd the symmetric center of the frequency domain. They are different 
%  when the number of meshes is odd or even,
nn = mod(Nxy,2);
if nn == 1
    NN = (Nxy-1)/(2*Nxy*abs(xInc));
else
    NN = 1/(2*abs(xInc));
end    

% Frequency domain coordinates after displacement
kx2 = 0:1/(Nxy*abs(xInc)):(Nxy-1)/(Nxy*abs(xInc));
ky2 = 0:1/(Nxy*abs(yInc)):(Nxy-1)/(Nxy*abs(yInc));
kx1 = kx2 - NN;
ky1 = ky2 - NN;
[KX,KY] = meshgrid(kx1,ky1);

% Cutoff frequency for the circular filter
cutoff = 1/filtThresh;

%  Distance from each point to the symmetric center
filterDistance = sqrt(KX.^2+KY.^2);


% Multiple data filtering are applied before Fourier analysis
zz = SmDetrend(z,L,1,resolution);
zzz = SmDetrend(zz,L,1,resolution);
zzzz = SmDetrend(zzz,L,1,resolution);

% Apply the FFT with filling zero.    
S = fft2(zzz,max(size(z)),max(size(z)));
if L > 200
    S = fft2(zzzz,max(size(z)),max(size(z)));
end    

SF = S;
SF(:,[1:2,round(size(SF,2)*0.5):size(SF,2)]) = nan;
SHalf = fftshift(SF);

% Apply the circular high-pass filter to the Fourier result.
indexC = find (filterDistance <= cutoff);
SHalf(indexC) = nan;

% Find values that are greater than 0.85 times the maximum.
SHalf = abs(SHalf);
indexSMax =find(abs(SHalf) >= (0.85*max(abs(SHalf(:)))));

% Get the wavelength and orientations.
dirArrayM =(atan2(KX(indexSMax),KY(indexSMax))*(180/pi));
lambdaArrayM = 1./sqrt(KX(indexSMax).^2+KY(indexSMax).^2);
powerArrayM = abs(SHalf(indexSMax));

indexD = find(lambdaArrayM >= 2*min(lambdaArrayM));
dirArrayM(indexD) = [];
lambdaArrayM(indexD) = [];
powerArrayM(indexD) = [];

stdDir = std(dirArrayM);
stdLambda = std(lambdaArrayM);

phiM = mean(dirArrayM);
lambdaM = mean(lambdaArrayM);
