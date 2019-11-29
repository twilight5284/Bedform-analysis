function [tc0,pc0] = zpf(tc,pc,L)
%% crossing points filter
% if the distance betwwen two adjacent crossing points is larger than 0.3
% times of the wavelength of interest, they would be filtered out.
% 

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
  tc = tc(1,find(~isnan(tc)));
  pc = pc(1,find(~isnan(tc)));
  zeroj=[];
  for zeroj = 1:2:(length(pc)-2)
    sig1 = pc(zeroj+2)-pc(zeroj);
    if sig1 < 0.3*L
      tc(zeroj) = NaN;
      tc(zeroj+1) = NaN;
      pc(zeroj) = NaN;
      pc(zeroj+1) = NaN;
    end
  end
  tc = tc(find(~isnan(tc)));
  pc = pc(find(~isnan(tc)));
  
  zerok=[];
  for zerok = 2:2:(length(pc)-2)
    sig1 = pc(zerok+2)-pc(zerok);
    if sig1 < 0.3*L
      tc(zerok) = NaN;
      tc(zerok+1) = NaN;
      pc(zerok) = NaN;
      pc(zerok+1) = NaN;
    end
  end
  tc = tc(find(~isnan(tc)));
  pc = pc(find(~isnan(tc)));
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
  tc = tc(find(~isnan(tc)));
  pc = pc(find(~isnan(tc)));
  tc0=tc;
  pc0=pc;
  
  