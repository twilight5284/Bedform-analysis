function [z0] = SmDetrend(z,L,d,resolution)
%% surface and profile filter
%
% INPUTS:
%     z, surface matrix;
%     L, wavelength of interest;
%     d, different filters; 
%        1 - exponential smoothing with a window of 0.75*L;
%        2 - moving average smoothing with a window of 3*L in the direction of rows;
%        3 - moving average smoothing with a window of 3*L in the direction of columns;
%
% OUTPUTS:
%     z0, surface matrix after filtering
%
% Author:  Li Wang
% Email:   twilight528400@hotmail.com

Num = size(z);

% exponential smoothing with a window of 0.75*L
if d == 1
  zz=[];  
  for i = 1:Num(1)
    zz(i,:) = smoothts(z(i,:),'e',fix(0.75*L/resolution))'; 
  end
  zzz=[];
  for j = 1:Num(2)
    zzz(:,j) = smoothts(zz(:,j),'e',fix(0.75*L/resolution)); 
  end
  z0=z-zzz;
  
% moving average smoothing with a window of 3*L in the direction of rows;
elseif d == 2
  zz=[];  
  for i = 1:Num(1)
    zz(i,:) = smooth(z(i,:),3*fix(L/resolution),'moving')'; 
  end
  z0=z-zz;
  
% moving average smoothing with a window of 3*L in the direction of columns;
else
  zzz=[];
  for j = 1:Num(2)
    zzz(:,j) = smooth(z(:,j),3*fix(L/resolution),'moving'); 
  end
  z0=z-zzz;
end
