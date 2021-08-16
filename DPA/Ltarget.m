function [L1,L2] = Ltarget(x,y)
%% 
% this function is applied to find the wavelength and wave height of the
% rhythmic features.
% 
% INPUTS
%     x - profile length
%     y - profile depth
% 
% OUTPUTS
%     L1 - the wavelength of the smallest rhythmic feature
%     L2 - the wavelength of the other rhythmic feature
%
% Author:  Li Wang
% Email:   twilight528400@hotmail.com

 H1=[];
 L1=[];
 H0=[];
 H2=[];
 L2=[];
    
 try
     maxtab = [];
     mintab = [];
     [maxtab, mintab] = peakdet(y, 0.1);
     H1=mean(maxtab(:,2))-mean(mintab(:,2));       
 catch
     try
        maxtab = [];
        mintab = [];
        [maxtab, mintab] = peakdet(y, 0.05);
        H1=mean(maxtab(:,2))-mean(mintab(:,2));
     catch
        maxtab = [];
        mintab = [];
        [maxtab, mintab] = peakdet(y, 0.01);
        H1=mean(maxtab(:,2))-mean(mintab(:,2));
     end
 end    
 L1=(max(x)-min(x))/(max(size(maxtab,1),size(mintab,1)));
    
 tb = (size(maxtab,1) - size(mintab,1));
 switch tb 
     case 0
         tb = tb*1;
     case tb < 0
         mintab(size(mintab,1),:)=[];
     case tb > 0
         maxtab(size(maxtab,1),:)=[];
 end    
    
 H0 = mean(abs(maxtab(:,2)-mintab(:,2)));
 Hdiff = H0 - H1;
 switch Hdiff
     case Hdiff >= 0
         H2 = H1+H0;
     otherwise
         H2 = 2*H1;
 end        
 
 [maxtab, mintab] = peakdet(y, H2);
 L2 = (max(x)-min(x))/max(size(maxtab,1),size(mintab,1));