function [BX,BY,x,y,z,X,Y,Z]=Partition(input,inputtype,L,resolution)
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bathymetry matrix partition
%
% INPUTS:
%     input - bathymetry matrix file;
%     inputtype - the type of input file;
%     L - the wavelength of interest;
%     resolution - the resolution of the orignal data.
%
% OUTPUTS:
%     BX, BY - the boundaries of subsets;
%     x, y, z - results of partition;
%     X, Y, Z - the grided matrixes of orignal data.
% 
% Author:  Li Wang
% Email:   twilight528400@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

load(input);

if inputtype == 2
    X = x;
    Y = y;
    Z = z;
else
    % data grid 
    [X,Y] = meshgrid(floor(min(min(data(:,1)))):resolution:ceil(max(max(data(:,1)))),floor(min(min(data(:,2)))):resolution:ceil(max(max(data(:,2)))));
    Z = griddata(data(:,1),data(:,2),data(:,3),X,Y,'linear');
    Z0 = griddata(data(:,1),data(:,2),data(:,3),X,Y,'nearest');
    Z(find(isnan(Z)))=Z0(find(isnan(Z)));
end

minx = min(X(1,:));
maxx = max(X(1,:));
miny = min(Y(:,1));
maxy = max(Y(:,1));
LengX = maxx-minx;
LengY = maxy-miny;

% the number of subsets
M1 = mod(LengX,10*L);
Q1 = fix(LengX/(10*L));
Mx = fix(M1/Q1);
M2 = mod(LengY,10*L);
Q2 = fix(LengY/(10*L));
My = fix(M2/Q2);

% calculate the boundaries of subsets
NumReso=1/resolution;
BorderX=[];
BorderY=[];
if Q1 <= 1
    BorderX(1,1) = 1;
    BorderX(1,2) = size(X,2);
    Q1 = 1;
else
    BorderX(1,1) = 1;
    BorderX(1,2*Q1) = size(X,2);
    for i = 2:Q1
        BorderX(1,(2*i-2)) = floor(NumReso*(10*(i-1)*L+Mx)+2*NumReso*L+1);
        BorderX(1,(2*i-1)) = ceil(NumReso*(10*(i-1)*L+Mx)-2*NumReso*L+1);
    end
end

if Q2 <= 1
    BorderY(1,1) = 1;
    BorderY(1,2) = size(X,1);
    Q2 = 1;
else
    BorderY(1,1) = 1;
    BorderY(1,2*Q2) = size(X,1);
    for j = 2:Q2
        BorderY(1,(2*j-2)) = floor(NumReso*(10*(j-1)*L+My)+2*NumReso*L+1);
        BorderY(1,(2*j-1)) = floor(NumReso*(10*(j-1)*L+My)-2*NumReso*L+1);
    end
end

% set buffer zones
BX=[];
BY=[];
BX=X(1,BorderX);
BY=Y(BorderY,1);
if size(BX,2)>2
  for k1=2:2:(size(BX,2)-1)
      BX(1,k1)=BX(1,k1)-2*L;
  end
  for k1=3:2:(size(BX,2))
      BX(1,k1)=BX(1,k1)+2*L;
  end
end
if size(BY,1)>2
  for k1=2:2:(size(BY,1)-1)
      BY(k1,1)=BY(k1,1)-2*L;
  end
  for k1=3:2:(size(BY,1))
      BY(k1,1)=BY(k1,1)+2*L;
  end
end

% the results of partition
x={};
y={};
z={};  
 for i=1:Q1
    for j=1:Q2
        Xa = BorderX(1,(2*i-1));
        Xb = BorderX(1,2*i);
        Ya = BorderY(1,(2*j-1));
        Yb = BorderY(1,2*j);
        x{j,i} = X(Ya:Yb,Xa:Xb);
        y{j,i} = Y(Ya:Yb,Xa:Xb);
        z{j,i} = Z(Ya:Yb,Xa:Xb);
    end
end