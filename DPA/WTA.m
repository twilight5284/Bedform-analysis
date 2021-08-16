function [WT] = WTA(X,Y,Z,filepath,name,resolution,L)
%% Wavelet analysis for bedform profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs the wavelet analysis for profiles.
% INPUTS:
%     X,Y,Z - grided matrixes of X, Y, and Z;
%     filepath - project folder path;
%     projectname - project folder name;    
%     L - the wavelength of interest;
%     resolution - the resolution of the input data.
% OUTPUTS
%     is saved in a mat file in the project folder.
%
% See also the function runWltAnalysis in Bedforms-ATM (Gutierrez et al., 2013)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TitleFig = name;

[raw,col]=size(Z);
lag = 1;
sampleFreq = resolution;

WT = {};

%% Wvelet transform
for j=1:fix(raw/lag)
  i = lag*j;
  temp = [];
  temp(:,1)=X(i,1:col);
  temp(:,2)=Y(i,1:col);
  temp(:,3)=Z(i,1:col);
  indexN = find(~isnan(temp(:,3)));
  tep = temp(indexN,:);
  n_3D_yJ=tep;
  signalData=tep(:,3)';
  signalDataX=tep(:,1);

% wavwlet transform
  [wltPower, fourierPeriod, coneOfInfluence] = ...
      waveletTransform (6, signalData, sampleFreq, 0.05);
% statistical analysis
  [wltSignif, globalSignif, globalWs,lag1,scaleAveSignif,scaleAvg] =...
      statisticsWlt...
      (6, signalData, sampleFreq, 0.05,...
      0.95,wltPower);
% getting peaks  
  [peaksCross,xPeakspos,yPeakspos] = getWltpeaks(globalWs,globalSignif,...
                                          fourierPeriod);
  
% save the result of wavelet transform    
    WT(j).i=i;
    WT(j).name=TitleFig;
    WT(j).signalDataX=signalDataX;
    WT(j).signalData=signalData;
    WT(j).globalWs=globalWs;
    WT(j).globalSignif=globalSignif;
    WT(j).coneOfInfluence=coneOfInfluence;
    WT(j).fourierPeriod=fourierPeriod;
    WT(j).wltPower=wltPower;
    WT(j).lag1=lag1;
    WT(j).wltSignif=wltSignif;
    WT(j).scaleAveSignif=scaleAveSignif;
    WT(j).scaleAvg=scaleAvg;
    WT(j).peaksCross=peaksCross;
    WT(j).xPeakspos=xPeakspos;
    WT(j).yPeakspos=yPeakspos;
    WT(j).sampleFreq=sampleFreq;
    WT(j).n_3D_yJ=n_3D_yJ;
end
%% Creating the *.mat file for print and statistical analysis.
     
nameOfMATFile = [filepath,'WT'];
save ([nameOfMATFile '.mat'], 'WT');
    