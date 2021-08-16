function [bedformProfile] = ScaleDiscr(filedir,WT,L,resolution)
%% bedform Discrimination
% Three scales of bedforms are separated by wavelet analyses and robust spline filters.
%
% INPUTS
%     filedir - project folder;
%     WT - outputs of Wavelet analysis;
%     L - the wavelength of interest;
%
% OUTPUTS
%     bedformProfile - separated bedform profiles
%
% Several technical improvements have been applied to this calculation, based on
% the Bedforms-ATM (Gutierrez et al., 2013).
% 
% See also README.txt

sParamPower = [-10:-1 0:2:50];
bedformProfile={};
  %%Discrimination
for j=1:length(WT)

    i = WT(j).i;
    TitleFig = WT(j).name;
    signalDataX = WT(j).signalDataX;
    signalData = WT(j).signalData;
    globalWs = WT(j).globalWs;
    globalSignif = WT(j).globalSignif;
    fourierPeriod = WT(j).fourierPeriod;
    sampleFreq = WT(j).sampleFreq;
    n_3D_yJ = WT(j).n_3D_yJ;
    
    % Profile info
    bedformProfile(j).name         = [TitleFig ' ' num2str(i)];
    bedformProfile(j).signalAbcise = signalDataX;
    bedformProfile(j).n            = signalData;  
    
    % Smooth versions of the signal
    smoothData = zeros(length(signalData),length(sParamPower));
    
    NsignalData = (signalData-mean(signalData))/std(signalData);
    
    % Actual value 12.5% of max amplitude and using normalized data
    param = 0.125 * ( max( NsignalData )-min(NsignalData));
    [~, mintab] = peakdet(NsignalData, param);
    
    troughSignal = signalData(mintab(:,1))';
   
    % 
    for nSParam = 1:length(sParamPower)

        Z = SMOOTHN(signalData,2^(sParamPower(nSParam)))'; % Base of the S Parameter
        delta = 0;
        troughSmooth = [];
        NsignalData = (Z-mean(Z))/std(Z);
        param = 0.125 * ( max( NsignalData )-min(NsignalData));
        [~, mintab] = peakdet(NsignalData, param);
        
        switch isempty(mintab)
            case 0
                troughSmooth = Z(mintab(:,1));
        end
        
        switch (length(troughSignal) - length(troughSmooth))
            case 0
                delta = mean(abs(troughSmooth-troughSignal));
        end
        
        Z = Z-delta; 
        smoothData(:,nSParam) = Z;
    end
    
    % L1 & H1, L2 & H2
    [L1,L2] = Ltarget(signalDataX,signalData);
     
    % period for target
    mintab = [];
    globalSW = log10(globalSignif)-log10(globalWs);
    H = 0.001*(max(globalSW)-min(globalSW));
    [~, mintab] = peakdet(globalSW, H);
    mintab(:,1)=fourierPeriod(1,mintab(:,1));
    
    % target
    [target1,target2,target3,scaleData] = TargetL(L1,L2,L,mintab,smoothData,signalData);
        
  %% Getting the 3rd level of the BFP 
    
  
    signal3 = [];    
    % border
    smoothWave = 2.^sParamPower;
    Border = min(find(smoothWave > L));
    [signal3] = TargetProf(Border,length(sParamPower),target3,scaleData,sampleFreq,globalSignif);
    
    % Getting the signal at level 3
    position3 = min(find(signal3(:,2)==min(signal3(:,2))));
    level3Position = signal3(position3,1)-1;
    % 3rd Level of th BFP
    bedformProfile(j).n33 = smoothData(:,level3Position);
%% Getting the 2nd level of the BFP 
    

    signal2 = [];
    forTarget1 = signalData'*ones(1,36)-smoothData;
    [signal2] = TargetProf(1,18,target1,forTarget1,sampleFreq,globalSignif);
        
    % Getting the signal at level 2
    position2 = min( find(signal2(:,2)==min(signal2(:,2))) );
    level2Position = signal2(position2,1);

    bedformProfile(j).n23trended = smoothData(:,level2Position);

    bedformProfile(j).n23 = bedformProfile(j).n23trended-...
                                 bedformProfile(j).n33;
    
    bedformProfile(j).n13 = bedformProfile(j).n'-...
                                 bedformProfile(j).n23trended;
     %%
    bedformProfile(j).n33 = bedformProfile(j).n33+...
                                 mean(bedformProfile(j).n23)+...
                                 mean(bedformProfile(j).n13);

    bedformProfile(j).n23 = bedformProfile(j).n23-mean(bedformProfile(j).n23);
    bedformProfile(j).n13 = bedformProfile(j).n13-mean(bedformProfile(j).n13);
    
              
    P1  = bedformProfile(j).n13;
    P2  = bedformProfile(j).n23;
    P3  = bedformProfile(j).n33;
    
   %% Smooth each Level
    P2S = smoothts(P2,'g',fix(0.25*L/resolution));
    P1S = P1+P2-P2S;
    bedformProfile(j).n13s = P1S;
    bedformProfile(j).n23s = P2S;
    bedformProfile(j).Y = n_3D_yJ(:,2);
    
    %% target 
    bedformProfile(j).target = roundn([target1 target2 target3],-2);
end
nameOfMATFile = [filedir,'bedformProfile'];
save ([nameOfMATFile '.mat'], 'bedformProfile');

function [signal] = TargetProf(a,b,targetL,profData,sampleFreq,globalSignif)

   levelIndex = 1;
   signal = [];
   for rrg = a:1:b
        % Wavelet transform
            [wltPower, fourierPeriod, ~] =...
            waveletTransform (6, profData(:,rrg)',...
                               sampleFreq   ,0.05 );

        % Wavelet significance
            [~, ~,globalWs,~,~,~] =...
            statisticsWlt( 6, profData(:,rrg)',sampleFreq,...
            0.05, 0.95, wltPower);
        
        % Estimate the closer peak to the target
        ScaledSig = log10(globalSignif) - log10(globalWs);
        
        ScaledSig = (ScaledSig-mean(ScaledSig))/std(ScaledSig);
                       
        param = 0.001*(max(ScaledSig)-min(ScaledSig));
        
        [~,mintab] = peakdet( ScaledSig, param );
        
        
        switch isempty(mintab)
            case 0
                mintab(:,1)=fourierPeriod(1,mintab(:,1));
                mintab(:,2)=(mintab(:,2)-min(mintab(:,2)))/abs(min(mintab(:,2)));
                maxp = mean(mintab((find(mintab(:,2) <= 0.5)),1));
        end
       
        switch isempty(maxp)
            case 0
                signal(levelIndex,1) = rrg;
                signal(levelIndex,2) = abs( targetL-maxp );
                levelIndex = levelIndex+1;
        end
        
   end 
