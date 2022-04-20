function [phiM,lambdaM,L_M,H_M,BPf,CT0]=Calculation(filepath,projectname,x,y,z,L,resolution)
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
[dirArrayM,lambdaArrayM,phiM0,stdDir,lambdaM,stdLambda] = wlFFT(x,y,z,L,resolution);

phiM=phiM0;
%% rotate the surface according to the dune orientation 
[X,Y,Z,Z0,Xb,Yb] = Rot_Grid(x,y,z,phiM,L,resolution);

%% wavelet transform
[WT]=WTA(X,Y,Z0,filepath,projectname,resolution,L);

%% bedform discrimination 
[bedformProfile] = ScaleDiscr(filepath,WT,L,resolution);


%% extract the crest and trough points, calculate the  dune morphological parameters
[CT,BP] = CT_BP(bedformProfile,Y,Z,Z0,L,resolution);

%% the original coordinate
CT0 = CT;
if (~isempty(CT))
   CT0(:,1) = ((CT(:,1))*cosd(90-phiM)-(CT(:,2))*sind(90-phiM))+Xb;
   CT0(:,2) = ((CT(:,1))*sind(90-phiM)+(CT(:,2))*cosd(90-phiM))+Yb;
end
BP0 = BP;
if (~isempty(BP))
  BP0(:,1) = ((BP(:,1))*cosd(90-phiM)-(BP(:,2))*sind(90-phiM))+Xb;
  BP0(:,2) = ((BP(:,1))*sind(90-phiM)+(BP(:,2))*cosd(90-phiM))+Yb;
else 
    BP=zeros(10,12);
end
  
%% calculate the average L and H

 L_M = mean(BP(:,3));
 std_L_M = std(BP(:,3));
 std_H_M = std(BP(:,6));
 H_M = mean(BP(:,6));


%% delete the dunes which are too small 
if (~isempty(BP0))
  for l=1:1:size(BP0,1)
   if  (BP0(l,6) <= (H_M/5))|(BP0(l,3) <= (L_M/5))
       BP0(l,:) = NaN;
   end
  end   
  BP0=BP0(find(~isnan(BP0(:,1))),:);
end  


%% delete the isolated points
if (~isempty(BP0))
    [BPf] = Distfilter(BP0,L);
else
    BPf=[];
end

%% give the dune moving direction 
if (~isempty(BP0))    
  if mean(BP0(:,7)) < 1
    Phi = 360-phiM;
  elseif mean(BP0(:,7)) == 1
    Phi = 'symmetrical';
  else
    Phi = phiM;
  end
else
    Phi=NaN;
end    


%% save all outputs

save([filepath,'output.mat'],...
    'dirArrayM','Phi','phiM','stdDir','lambdaArrayM','lambdaM','stdLambda',...
    'x','y','z','X','Y','Z','Z0','Xb','Yb',...
    'bedformProfile',...
    'CT','BP','CT0','BP0','BPf',...
    'L_M','std_L_M','H_M','std_H_M'); 

