function [LAmean,LSAmean,LAmax] = LeeAngle(LAC,resolution,Or)
%% lee slope angle calculation
% INPUTS
%     LAC - lee slope curve
%     resolution - the resolution of the profile
%     S - the slope of the seabed topography
%
% OUTPUTS
%     LAmean: the average angle of the lee slope
%     LSAmean: the average angle of the lee slope between knee points
%     LAmax: the max angle of the lee slope
%
% Author:  Li Wang
% Email:   twilight528400@hotmail.com

%%
DL1 = [];
LA = [];
DL2 = [];
LAmax = [];
LAMax = [];
loc = [];

% LAC(:,2) = smooth(LAC(:,2),size(LAC,1)/3,'moving'); % ²»Æ½»¬¼ÆËãlee slope

DL1 = Or * diff(LAC(:,2))/resolution;

LA = atand(abs(DL1));
LA0 = LA.*(DL1./abs(DL1));
if size(LA,1) < 3
    LAMax = find(LA0 == max(LA0));
else
    LAMax = find(LA0 == max(LA0(2:end-1)));
end    
loc = mean(LAMax);
LAC(:,3) = smooth(LAC(:,2),max(1,floor(size(LAC,1)/2.5)));
DL2 = diff(LAC(:,3),2)/resolution;
if size(DL2,1) < 3 
    LAmean = mean(LA0(find(LA0>0)));
    LSAmean = mean(LA0(find(LA0>0)));
    LAmax = max(LA0);
else 
  if loc == size(LA,1)
    loca = loc-1;
  elseif loc ==  (size(LA,1)-1)
    loca = loc;
  else
    loca = loc+1;
  end    
  if loc == 1
    locb = loc;
  else
    locb = loc-1;
  end
  if loc > size(LA,1)/5 & loc < 4*size(LA,1)/5
      Lmax = find(DL2(locb:floor(4*size(DL2,1)/5)) == min(DL2(locb:floor(4*size(DL2,1)/5)))) + locb - 1;
      Lmin = find(DL2(ceil(size(LA,1)/5):loca,1) == max(DL2(ceil(size(LA,1)/5):loca)));
  else    
      Lmax = find(DL2(locb:size(DL2,1)) == min(DL2(locb:size(DL2,1)))) + locb - 1;
      Lmin = find(DL2(1:loca,1) == max(DL2(1:loca)));
  end    
  if (size(Lmin) == size(Lmax)) & (Lmin == Lmax)
    LAmean = mean(LA0(find(LA0>0)));
    LSAmean = mean(LA0(find(LA0>0)));
    LAmax = max(LA0);
  else
    LAmean = mean(LA0(find(LA0>0)));
    LSAmean = mean(LA0(find(LA0(max(Lmin):min(Lmax)+1)>0)));
    LAmax = max(LA0);
  end
end
if isnan(LSAmean)
    LSAmean = LAmean;
end
