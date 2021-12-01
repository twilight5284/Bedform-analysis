%% Dome
%% the first step: find the wavelength(s) of interest

inputdata = 'testdata_type_2.mat';             % a file in format of mat, for example, 'H:\bedform\input.mat', 
inputtype = 2;                                 % the type of input matrix, it could be 1 or 2; 
                                               % type 1: a matrix named 'data' with 3 columns, they are x, y, and z. For example 'testdata_type_1.mat'  
                                               % type 2: three matrixes with same rows and columns£¬they are grided x, y, and z matrixes. For example 'testdata_type_2.mat'
resolution = 1;                                % resolution - the data resolution

[LT, T, PHID] = LT(inputdata,inputtype,resolution); % get the wavelength(s) of interest through 2D Fourier analysis

%% the second step: bedform analysis

filedir = 'H:\';                               % the directory of the folder where you want to save your results
projectname='Demo';                            % the projectname
L=208;                                         % the wavelength of interest, you can get it from the LT function above or by manual selection
resolution=1;                                  % the data resolution
input='testdata_type_2.mat';                   % the input data, a file in format of mat, for example, 'H:\bedform\input.mat'
inputtype=2;                                   % the type of input matrix, it could be 1 or 2; 
                                               % type 1: a matrix named 'data' with 3 columns, they are x, y, and z. For example 'testdata_type_1.mat'  
                                               % type 2: three matrixes with same rows and columns£¬they are grided x, y, and z matrixes. For example 'testdata_type_2.mat'   
DPA(filedir,projectname,input,inputtype,L,resolution); % bedform analysis 

