function [phiM,lambdaM,L_M,H_M,BPf]=Calculation(filepath,projectname,x,y,z,L,resolution)
%% Dune Parameters Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dune regional pattern are calculated by two dimensional Fourier analysis;
% Individual dune parameters are calculated by wavelet transform and
% zero-crossing analysis.
%
% INPUTS
%     filepath - project folder path;
%     projectname - project folder name;
%     x,y,z - grided matrixes of x, y, and z;
%     L - the wavelength of interest;
%     resolution - the resolution of the input data.
%
% OUTPUTS
%     phiM - the regional orientation of dunes;
%     lambdaM - the regional wavelength of dunes;
%     BPf - the parameters of all individual dunes;
%     L_M, H_M - the average wavelength amd wave height of all dunes;
% All outputs are saved as 'output.mat' in the project folder.
%
% Author:  Li Wang
% Email:   twilight528400@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate the dune regional pattern by 2-D DFT
[dirArrayM,lambdaArrayM,phiM,stdDir,lambdaM,stdLambda] = wlFFT(x,y,z,L);

%% rotate the surface according to the dune orientation 
[X,Y,Z,Z0,Xb,Yb] = Rot_Grid(x,y,z,phiM,L,resolution);

%% wavelet transform
[WT]=WTA(X,Y,Z0,filepath,projectname,resolution,L);

%% bedform discrimination 
[bedformProfile] = ScaleDiscr(filepath,WT,L);


%% extract the crest and trough points, calculate the  dune morphological parameters
[CT,BP] = CT_BP(bedformProfile,Y,Z,Z0,L,resolution);

%% the original coordinate
CT_O = [];
if (~isempty(CT))
   CT_O(:,1) = ((CT(:,1))*cosd(90-phiM)-(CT(:,2))*sind(90-phiM))+Xb;
   CT_O(:,2) = ((CT(:,1))*sind(90-phiM)+(CT(:,2))*cosd(90-phiM))+Yb;
   CT_O(:,3) = CT(:,3);
end
BP0 = [];
if (~isempty(BP))
  BP0(:,1) = ((BP(:,1))*cosd(90-phiM)-(BP(:,2))*sind(90-phiM))+Xb;
  BP0(:,2) = ((BP(:,1))*sind(90-phiM)+(BP(:,2))*cosd(90-phiM))+Yb;
  BP0(:,3:12) = BP(:,3:12);
end


%% calculate the average L and H
if (~isempty(BP))
  L_M = mean(BP(:,3));
  std_L_M = std(BP(:,3));

  std_H_M = std(BP(:,6));
  H_M = mean(BP(:,6));
end
%% delete the dunes which are too small 
for l=1:1:size(BP0,1)
   if  (BP0(l,6) <= (H_M/5))|(BP0(l,3) <= (L_M/5))
       BP0(l,:) = NaN;
   end
end   
BP0=BP0(find(~isnan(BP0(:,1))),:);

%% delete the isolated points
[BPf] = Distfilter(BP0,L);

%% give the dune moving direction 
if mean(BP0(:,7)) < 1
    Phi = 360-phiM;
elseif mean(BP0(:,7)) == 1
    Phi = 'symmetrical';
else
    Phi = phiM;
end 

%% save all outputs
save([filepath,'output.mat'],...
    'dirArrayM','Phi','phiM','stdDir','lambdaArrayM','lambdaM','stdLambda',...
    'x','y','z','X','Y','Z','Z0','Xb','Yb',...
    'bedformProfile',...
    'CT','BP','CT_O','BP0','BPf',...
    'L_M','std_L_M','H_M','std_H_M'); 

function [BPf] = Distfilter(BP0,L)
DistLag = max(2,floor(L/50));
DuneCrest = BP0(:,1:2);
DistNum = [];
for i = 1:size(BP0,1)
    Dist = [];
    DC0 = DuneCrest(i,:);
    Dist = sqrt((DuneCrest(:,1)-DC0(1)).^2+(DuneCrest(:,2)-DC0(2)).^2);
    DistNum(i,1) = length(find(Dist < 10*DistLag));
end  
DistFil = find(DistNum < 10);
BPf=BP0;
BPf(DistFil,:)=[];
