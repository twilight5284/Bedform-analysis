function [X,Y,Z,Z0,Xb,Yb] = Rot_Grid(x,y,z,phiM,L,resolution)
%% surface rotation and regrid
%
% INPUTS:
%     x,y,z, surface matrixes
%     phiM, rotation angle
%     L, wavelength of interest
%     resolution, the data resolution
%
% OUTPUTS:
%     X,Y,Z,Z0, surface matrixes after rotation
%     Xb,Yb, the coordinate of the lower left corner of the original surface
% 
% Author:  Li Wang
% Email:   twilight528400@hotmail.com

% find the coordinates in the lower left corner
Xb = min(min(x));
Yb = min(min(y));
x0 = x-Xb;
y0 = y-Yb;
% rotate the coordinates
if phiM ~= 90 
  X0=(x0)*cosd(phiM-90)-(y0)*sind(phiM-90);
  Y0=(x0)*sind(phiM-90)+(y0)*cosd(phiM-90);
else
  X0=x0;
  Y0=y0;
end
[X,Y] = meshgrid(floor(min(min(X0))):resolution:ceil(max(max(X0))),floor(min(min(Y0))):resolution:ceil(max(max(Y0))));

% profile filtering 
if 45 < phiM < 135
  z0 = SmDetrend(z,L,2,resolution);
else
  z0 = SmDetrend(z,L,3,resolution);
end

% regird the surface depth
Z = griddata(X0,Y0,z,X,Y,'linear');
Z0 = griddata(X0,Y0,z0,X,Y,'linear');
indexL = find(resolution*sum((~isnan(Z)),2) < 3*L);
X(indexL,:)=[];
Y(indexL,:)=[];
Z(indexL,:)=[];
Z0(indexL,:)=[];




