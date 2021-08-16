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
%       3 - wavelength (straight line of two troughs£©
%       4 - left wavelength
%       5 - right wavelength
%       6 - wave height (perpendicular)
%       7 - wavelength (horizontal)
%       8 - wave height (vertical, crest to the average depth of two troughs)
%       9 - wave height (vertical, crest to the line of two troughs)
%       10 - dune asymmetry
%       11 - dune stepness
%       12 - average lee slope angle
%       13 - average lee slope angle between knick points
%       14 - maximum lee slope angle
%       15 - depth of dune crest
%       16 - slope of the straight line of two troughs
% Author:  Li Wang
% Email:   twilight528400@hotmail.com

%%
rawB = size(B,1);
m=0;
begin=1;
if B(1,3)>B(2,3)
    begin=2;
end    
Num=floor((rawB-begin)/2);
BQP = zeros(Num,15);
for l = begin:2:(2*Num-2+begin)
    BQP(m+1,1) = B(l+1,1);             %  x of dune crest
    BQP(m+1,2) = B(l+1,2);             %  y of dune crest
    BQP(m+1,15) = B(l+1,3);            %  depth of dune crest
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
    BQP(m+1,7)=B(l+2,1)-B(l,1);                      
    BQP(m+1,8)=B(l+1,3)-(B(l,3)+B(l+2,3))/2; 
    BQP(m+1,9)=B(l+1,3)-B(l,3)-(B(l+1,1)-B(l,1))*(B(l+2,3)-B(l,3))/(B(l+2,1)-B(l,1)); 
    BQP(m+1,10)=(BQP(m+1,4)-BQP(m+1,5))/BQP(m+1,3);   %  asymmetry
    BQP(m+1,11)=BQP(m+1,6)/BQP(m+1,3);                %  stepness
   
    S = atand((B(l+2,3)-B(l,3))/(B(l+2,1)-B(l,1)));  
    LAC=[];
    if BQP(m+1,10)>=0                                    %  find the lee slope
        LAC(:,1) = pos(find(pos == B(l+1,1)):find(pos == B(l+2,1)));
        LAC(:,2) = Data(find(pos == B(l+1,1)):find(pos == B(l+2,1))); 
        S = S;
        Or = -1;
    else
        LAC(:,1) = pos(find(pos == B(l,1)):find(pos == B(l+1,1)));
        LAC(:,2) = Data(find(pos == B(l,1)):find(pos == B(l+1,1)));   
        S = -S;
        Or = 1;
    end
    BQP(m+1,16)=S;
    
    [LAmean,LSAmean,LAmax] = LeeAngle(LAC,resolution,Or);  % calculate the lee slope angle 
    TT = (LAmean > 30 | LSAmean > 35 | isnan(LSAmean));
    JJ = 0;
    while (TT) & (JJ<4)
        [LAmean,LSAmean,LAmax] = FilterLA(LAC,resolution,Or);
        JJ=JJ+1;
    end
    BQP(m+1,12) = LAmean;
    BQP(m+1,13) = LSAmean;
    BQP(m+1,14) = LAmax;
    
    if  (BQP(m+1,6) <= 0.075)|(BQP(m+1,3) >= 3.5*L)|(BQP(m+1,12) <= 0)
        BQP(m+1,:) = NaN;
        BQP(m+1,:) = NaN;
    end
    m=m+1;
end
BQP=BQP(find(~isnan(BQP(:,1))),:);

function [LAmean,LSAmean,LAmax] = FilterLA(LAC,resolution,Or)
    LAC(:,2) = smooth(LAC(:,2),floor(size(LAC,1)/3)); % ²»Æ½»¬¼ÆËãlee slope
    [LAmean,LSAmean,LAmax] = LeeAngle(LAC,resolution,Or);