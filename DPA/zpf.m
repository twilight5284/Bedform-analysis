function [tc0,pc0] = zpf(tc,pc,L)
%% crossing points filter
% if the distance betwwen two adjacent crossing points is larger than 0.3
% times of the wavelength of interest, they would be filtered out.
% 

  [tc1,pc1] = L3L1(tc,pc,L);
  [tc2,pc2] = L3L1(tc1,pc1,L);
  [tc3,pc3] = L3L2(tc2,pc2,L);
  [tc4,pc4] = L3L2(tc3,pc3,L);
  
  tc0=tc4;
  pc0=pc4;
  function [tc0,pc0] = L3L1(tc,pc,L) 
  zeroi=[];
  for zeroi = 1:2:(length(pc)-2)
    sig1 = pc(zeroi+2)-pc(zeroi);
    if sig1 < 0.3*L
      tc(zeroi) = NaN;
      tc(zeroi+1) = NaN;
      pc(zeroi) = NaN;
      pc(zeroi+1) = NaN;
    end
  end
  tc0 = tc(1,find(~isnan(tc)));
  pc0 = pc(1,find(~isnan(tc)));      
 
  function [tc0,pc0] = L3L2(tc,pc,L) 
  zerol=[];
  for zerol = 2:2:(length(pc)-2)
    sig1 = pc(zerol+2)-pc(zerol);
    if sig1 < 0.3*L
      tc(zerol) = NaN;
      tc(zerol+1) = NaN;
      pc(zerol) = NaN;
      pc(zerol+1) = NaN;
    end
  end
  tc0 = tc(find(~isnan(tc)));
  pc0 = pc(find(~isnan(tc)));