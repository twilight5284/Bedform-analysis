function [target1,target2,target3,scaleData] = TargetL(L1,L2,L,mintab,smoothData,signalData)
%%
% this function is applied to find the target wavelengths of a bedform
% profile.
% 
% INPUTS
%     L1 - the wavelength of the smallest rhythmic feature
%     L2 - the wavelength of the other rhythmic feature
%     L - the wavelength of interest
%     mintab - output of peakdet
%     smoothData - smoothed profile
%     signalData - bedform profile
%
% OUTPUTS
%     target1 - the smallest target wavelength
%     target2 - the target wavelength close to the wavelength of interest
%     target3 - the largest target wavelength
%     scaleData - the smoothed profile that best represents the dunes with a wavelength of target2
%
% Author:  Li Wang
% Email:   twilight528400@hotmail.com
 mintab1=[];
 mintab2=[];
 minp=[];
 maxp=[];
 mintab3=[];
 mintab0 = mintab;
 Ldiff = max(L/L1,L1/L) - max(L/L2,L2/L);
 
 
 switch (Ldiff < 0)
     case 1  
         % target2
         a = min(mintab(:,2));
         a = a/abs(a);
         switch  a
             case -1
                 mintab2=mintab((find(mintab(:,2)<=0)),:);
                 [~,target2] = searchclosest(mintab2(:,1),L);
                 if max(target2/L,L/target2)>5
                     mintab((find(mintab(:,2)<=0)),:)=[];
                     mintab2=[];
                     mintab2=mintab;
                     mintab2(:,2)=(mintab2(:,2)-min(mintab2(:,2)))/abs(min(mintab2(:,2)));
                     maxp = mintab2((find(mintab2(:,2) <= 0.5)),1);
                     [~,target2] = searchclosest(maxp,L);  
                 end  
             otherwise
                 mintab2=mintab(find((mintab(:,1) < 3*L) & (mintab(:,1) > L/5)),:);
                 mintab2(:,2)=(mintab2(:,2)-min(mintab2(:,2)))/abs(min(mintab2(:,2)));
                 maxp = mintab2((find(mintab2(:,2) <= 0.5)),1);
                 [~,target2] = searchclosest(maxp,L);  
         end
         % target1
         mintab1 = mintab((find(mintab(:,1)<(min(target2,L)/4))),:);
         b = size(mintab1,1);
         switch b
             case 0
                 target1=min(L1,min(target2,L)/4);
             otherwise
                 mintab1(:,2)=(mintab1(:,2)-min(mintab1(:,2)))/abs(min(mintab1(:,2)));
                 minp = mintab1((find(mintab1(:,2) <= 0.5)),1);
                 [~,target1] = searchclosest(minp,min(L1,min(target2,L)/4));
         end
         % target3
         mintab3 = mintab((find(mintab(:,1)>(3*max(target2,L1)))),:);
         c = size(mintab3,1);
         switch c
             case 0
                 target3=max(3*max(target2,L1),3.5*L);
             otherwise
                 target3=mintab3((find(mintab3(:,2)==(min(mintab3(:,2))))),1);
                 target3=max(target3, 3.5*L);
         end                 

     case 0
         % target2
         a = min(mintab(:,2));
         a = a/abs(a);
         switch  a
             case -1
                 mintab2=mintab((find(mintab(:,2)<=0)),:);
                 [~,target2] = searchclosest(mintab2(:,1),L);
                 if max(target2/L,L/target2)>5
                     mintab((find(mintab(:,2)<=0)),:)=[];
                     mintab2=[];
                     mintab2=mintab;
                     mintab2(:,2)=(mintab2(:,2)-min(mintab2(:,2)))/abs(min(mintab2(:,2)));
                     maxp = mintab2((find(mintab2(:,2) <= 0.5)),1);
                     [~,target2] = searchclosest(maxp,L);  
                 end    
             otherwise    
                 mintab2=mintab;
                 mintab2(:,2)=(mintab2(:,2)-min(mintab2(:,2)))/abs(min(mintab2(:,2)));
                 maxp = mintab2((find(mintab2(:,2) <= 0.5)),1);
                 [~,target2] = searchclosest(maxp,L);  
         end
         
         % target1
         mintab1 = mintab((find(mintab(:,1)<(min(target2,L)/4))),:);
         b = size(mintab1,1);
         switch b
             case 0
                 target1 = L1;
             otherwise
                 mintab1(:,2)=(mintab1(:,2)-min(mintab1(:,2)))/abs(min(mintab1(:,2)));
                 minp = mintab1((find(mintab1(:,2) <= 0.5)),1);
                 [~,target1] = searchclosest(minp,L1);
         end    
         % target3
         mintab3 = mintab((find(mintab(:,1)>(2*max([target2,L2,L])))),:);
         c = size(mintab3,1);
         switch c
             case 0
                 target3=max(2*max([target2,L2,L]),3.5*L);
             otherwise
                 target3=mintab3((find(mintab3(:,2)==(min(mintab3(:,2))))),1);
                 target3=max(target3, 3.5*L);
         end        
 end
    
 % target = L3 
 mintab0(:,2)=(mintab0(:,2)-min(mintab0(:,2)))/abs(min(mintab0(:,2)));
 target0 = mean(mintab((find(mintab0(:,2) <= 0.5)),1));
 scaleData=[];
 d = abs(target0-target2) - abs(target0-target3);
 d = d/abs(d);
 switch d
     case 1
         target3 = target2;
         scaleData = signalData';
         scaleData = scaleData(:,ones(1,36));
         scaleData = scaleData-smoothData;         
     otherwise
         scaleData = smoothData;
 end