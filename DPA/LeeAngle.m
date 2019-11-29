function [LAmean,LSAmean,LAmax] = LeeAngle(LAC,resolution,S)
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

LAC(:,2) = smooth(LAC(:,2),size(LAC,1)/3,'moving');

DL1 = diff(LAC(:,2))/resolution;

LA = atand(abs(DL1));
LA0 = LA +S;
LAMax = find(LA == max(LA));
loc = mean(LAMax);
DL2 = diff(LAC(:,2),2)/resolution;
if size(DL2,1) < 3 
    LAmean = mean(LA0);
    LSAmean = mean(LA0);
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
  if mean(DL1) > 0  
      Lmax = find(DL2(locb:size(DL2,1)) == min(DL2(locb:size(DL2,1)))) + locb - 1;
      Lmin = find(DL2(1:loca,1) == max(DL2(1:loca)));
  else
      Lmin = find(DL2(1:loca) == min(DL2(1:loca)));
      Lmax = find(DL2(locb:size(DL2,1)) == max(DL2(locb:size(DL2,1)))) + locb - 1;
  end
  if Lmin == Lmax
    LAmean = mean(LA0);
    LSAmean = mean(LA0);
    LAmax = max(LA0);
  else
    LAmean = mean(LA0);
    LSAmean = mean(LA0(max(Lmin):min(Lmax)+1));
    LAmax = max(LA0);
  end
end
if LSAmean < LAmean
    LSAmean = LAmean;
end
