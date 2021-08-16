function DPA(filedir,projectname,input,inputtype,L,resolution)
%% Dune Parameters Analysis (DPA) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DPA (Dune Parameters Analysis) quantitatively measures dunes.
% The geometry parameters (crest orientation, crest location, depth, 
% wavelength, wave height, asymmetry, lee slope angele, etc.) of each dune 
% could be analysed automatically with the bathymetry matrix. 
% Outputs would be saved in the project folder.
%
% INPUTS
%     filedir - project directory, for example, 'H:\bedform\', 
%               the file path must already exist;
%     projectname - name of project folder, for example, '201803',
%                   a new folder would be created to save the outputs; 
%     input - a file in format of mat, for example, 'H:\bedform\input.mat', 
%         type 1: a matrix named 'data' with 3 columns, they are x, y, and z.  
%         type 2: three matrixes with same rows and columns£¬they are 
%                 grided x, y, and z matrixes.
%     inputtype - the type of input matrix, it could be 1 or 2;  
%     L - the wavelength of the dunes of interest. Using the LT function,
%         you can find the wavelength values of interest;
%     resolution - the resolution of the data. 
%
% OUTPUTS
%     x0, y0, z0 - matrixes of x, y, and z after grid;
%     PhiM, LambdaM - the regional pattern of dunes in ecah subset calculated 
%                   by two dimensional Fourier analysis, PhiM is the orientation, 
%                   LambdaM is the wavelength. 
%     LM, HM - average parameters of all dunes in each subset. LM is the wavelength, 
%                   HM is the wave height.
%     BPall, BPALL - dune parameters of each dune. BPALL is the union of each
%                   subset (BPall).
% EXAMPLE:
% DPA('H:\','test','H:\input.mat',2,20,0.5);
% 
% See also README.txt
%
% Version: 1.0, 24/05/2019
% Author:  Li Wang
% Email:   twilight528400@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% divide the bathymetry matrix into subsets
[BX,BY,X,Y,Z,x0,y0,z0]=Partition(input,inputtype,L,resolution);
NumX=size(X,1);
NumY=size(X,2);

%% create the project folder
mkdir(filedir,projectname);
filepath1 = [filedir,projectname,'\'];

PhiM=[];
LambdaM=[];
LM=[];
HM=[];
BPall={};

%% calculate the dune parameters of each subset
for i=1:NumX
    for j=1:NumY
        
        % read the x, y, and z of each subset
        x=X{i,j};
        y=Y{i,j};
        z=Z{i,j};
        
        % create project folder of each subset                                   
        filename = [num2str(i),'_',num2str(j)];        
        mkdir(filepath1,filename);
        
        % save input file of each subset
        filepath2 = [filepath1,filename,'\'];
        save([filepath2,'input.mat'],'x','y','z');
        
        % calculate the dune parameters of each subset
        [phiM,lambdaM,L_M,H_M,BPf,CT0]=Calculation(filepath2,filename,x,y,z,L,resolution);
        
        % save the output of each subset
        PhiM(i,j) = phiM;
        LambdaM(i,j) = lambdaM;
        LM(i,j) = L_M;
        HM(i,j) = H_M;
        BPall{i,j} = BPf;
        CTall{i,j} = CT0;
        NS = num2str((i-1)*NumY + j);
        NT = num2str(NumX*NumY);
        disp(['The subset ' num2str(i) '_' num2str(j) ' is finished. (' NS '/' NT ')']);
    end
end

%% crest orientation
PhiM = PhiM - 90;
PhiM(find(PhiM < 0)) = PhiM(find(PhiM < 0)) + 180;
%% delete the buffer overlaps
BPALL=[];
for i = 1:size(BPall,1)
    for j = 1:size(BPall,2)
        Temp = BPall{i,j};
        if isempty(Temp)
            continue;
        end    
        Xa=BX(1,(2*j-1));
        Xb=BX(1,(2*j));
        Ya=BY((2*i-1),1);
        Yb=BY((2*i),1);
        Temp(find(Temp(:,1)<Xa),:) = [];
        Temp(find(Temp(:,1)>Xb),:) = [];
        Temp(find(Temp(:,2)<Ya),:) = [];
        Temp(find(Temp(:,2)>Yb),:) = [];
        BPALL = cat(1,BPALL,Temp);
    end
end
CTALL=[];
for i = 1:size(CTall,1)
    for j = 1:size(CTall,2)
        Temp = CTall{i,j};
        if isempty(Temp)
            continue;
        end   
        Xa=BX(1,(2*j-1));
        Xb=BX(1,(2*j));
        Ya=BY((2*i-1),1);
        Yb=BY((2*i),1);
        Temp(find(Temp(:,1)<Xa),:) = [];
        Temp(find(Temp(:,1)>Xb),:) = [];
        Temp(find(Temp(:,2)<Ya),:) = [];
        Temp(find(Temp(:,2)>Yb),:) = [];
        CTALL = cat(1,CTALL,Temp);
    end
end
%% save results of all
mkdir(filepath1,'all');
filepath3 = [filepath1,'all','\'];
save([filepath3,'output.mat'],'x0','y0','z0','PhiM','LambdaM','LM','HM','BPall','BPALL','CTALL');