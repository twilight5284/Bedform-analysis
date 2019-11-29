function [BQP] = bqp(B,L,Data,pos,resolution)
%% calculate the dune geometry parameters with dune crests and troughs   
% INPUTS:
%     B - locations of crests and troughs
%     L - wavelength of interest
%     Data - depth of profile data
%     pos - length of profile data
%     resolution - the resolution of profile data
% OUTPUTs:
%     BQP - the dune geometry parameters
%      columns:
%       1 - x of dune crest
%       2 - y of dune crest
%       3 - wavelength
%       4 - left wavelength
%       5 - right wavelength
%       6 - wave height
%       7 - dune asymmetry
%       8 - dune stepness
%       9 - average lee slope angle
%       10 - average lee slope angle between knick points
%       11 - maximum lee slope angle
%       12 - depth of dune crest
%
% Author:  Li Wang
% Email:   twilight528400@hotmail.com

%%
rawB = size(B,1);
BQP = zeros(abs((rawB-1)/2),11);
m=0;
for l = 1:2:(rawB-2)
    BQP(m+1,1) = B(l+1,1);             %  x of dune crest
    BQP(m+1,2) = B(l+1,2);             %  y of dune crest
    BQP(m+1,12) = B(l+1,3);            %  depth of dune crest
    BQP(m+1,3)=sqrt((B(l+2,1)-B(l,1)).^2+(B(l+2,3)-B(l,3)).^2);   %  wavelength
    a=B(l+1,:)-B(l,:);
    b=B(l+2,:)-B(l,:);
    c=B(l+1,:)-B(l+2,:);
    a(:,2)=[];
    b(:,2)=[];
    c(:,2)=[];
    d=-b;
    BQP(m+1,4)=(dot(a,b))/(sqrt(b(1,1).^2+b(1,2).^2));     
    BQP(m+1,5)=(dot(c,d))/(sqrt(d(1,1).^2+d(1,2).^2));     
    BQP(m+1,6)=abs(det([b;a]))/norm(b);              %  wave height
    BQP(m+1,7)=(BQP(m+1,4)-BQP(m+1,5))/BQP(m+1,3);   %  asymmetry
    BQP(m+1,8)=BQP(m+1,6)/BQP(m+1,3);                %  stepness
   
    S = atand((B(l+2,3)-B(l,3))/(B(l+2,1)-B(l,1)));  
    LAC=[];
    if BQP(m+1,7)>=0                                    %  find the lee slope
        LAC(:,1) = pos(find(pos == B(l+1,1)):find(pos == B(l+2,1)));
        LAC(:,2) = Data(find(pos == B(l+1,1)):find(pos == B(l+2,1)));
        S = S;
    else
        LAC(:,1) = pos(find(pos == B(l,1)):find(pos == B(l+1,1)));
        LAC(:,2) = Data(find(pos == B(l,1)):find(pos == B(l+1,1)));
        S = -S;
    end
    
    [LAmean,LSAmean,LAmax] = LeeAngle(LAC,resolution,S);  % calculate the lee slope angle 
    BQP(m+1,9) = LAmean;
    BQP(m+1,10) = LSAmean;
    BQP(m+1,11) = LAmax;
    
    if  (BQP(m+1,6) <= 0.1)|(BQP(m+1,3) >= 3*L)|(BQP(m+1,9) <= 0)
        BQP(m+1,:) = NaN;
        BQP(m+1,:) = NaN;
    end
    m=m+1;
end
BQP=BQP(find(~isnan(BQP(:,1))),:);