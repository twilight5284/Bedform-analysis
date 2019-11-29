function [LT,T,PHID] = LT(input,inputtype)
%% find the wavelength of interest
% A series of 2D DFT calculations are used to find out the wavelength of
% interest.
% The bathymetric map and a bedform profile are ploted to determine the 
% reliability of the results. 
% 
% INPUTS
%     input - a file in format of mat, for example, 'H:\bedform\input.mat', 
%         type 1: a matrix named 'data' with 3 columns, they are x, y, and z.  
%         type 2: three matrixes with same rows and columns£¬they are 
%                 grided x, y, and z matrixes.
%     inputtype - the type of input matrix, it could be 1 or 2; 
%
% OUTPUTS
%     LT - bedform wavelength of interest
%     T - bedform wavelength and crest orientation calculated by 2D DFT 
%     PHID - bedform crest orientation
%
% EXAMPLE:
% [LT,T,peaks,PHID] = LT('H:\matlab files\bedform\input.mat',2);
%
% Version: 1.0, 10/09/2019
% Author:  Li Wang
% Email:   twilight528400@hotmail.com

%%
load (input);
resolution =1;
if inputtype == 2
    X = x;
    Y = y;
    Z = z;
else
    % data grid 
    [X,Y] = meshgrid(floor(min(min(data(:,1)))):resolution:ceil(max(max(data(:,1)))),...
        floor(min(min(data(:,2)))):resolution:ceil(max(max(data(:,2)))));
    Z = griddata(data(:,1),data(:,2),data(:,3),X,Y,'linear');
    Z0 = griddata(data(:,1),data(:,2),data(:,3),X,Y,'nearest');
    Z(1,:) = Z0(1,:);
    N1 = size(Z,2);
    N2 = size(Z,1);
    Z(N2,:) = Z0(N2,:);
    Z(:,N1) = Z0(:,N1);
    Z(1,:) = Z0(1,:);
    Z(:,1) = Z0(:,1);
end

% plot the bathymetry
figure;
surf(X,Y,Z)
shading interp;
set(gca,'DataAspectRatio',[1 1 0.1]);

Xmin = min(min(X));
Xmax = max(max(X));
Ymin = min(min(Y));
Ymax = max(max(Y));
t = ceil(max(Xmax-Xmin,Ymax-Ymin)/3);  
T=[];
i=1;
a=[1.5:1:20];
b=1.5.^a;
b=ceil(b);    % Exponential sequence
while b(1,i) < min(t,1000)
    % Calculation of possible regional dune wavelengths with 2D DFT
    [~,~,phiM,~,lambdaM,~] = wlFFT(X,Y,Z,b(1,i));  
    T(i,1)=b(1,i);
    T(i,2)=phiM;
    T(i,3)=lambdaM;
    disp(i) 
    i=i+1;
end    

TT=(T(2:end,3)-T(1:end-1,3))./T(2:end,3);
% Those their relative difference does not exceed 1/10 are considered to be the same regional dune wavelength.
tt1=find(TT<0.1);      
tt2=tt1+1;
tt=[tt1;tt2];
tt0=unique(tt);
t0=find((tt0(2:end)-tt0(1:end-1))>1);   
LT=[];
switch length(t0)  
    case 0                %   There is only one scale of dunes.
        LT(1,:) = mean(T(tt0,3));    %   wavelength
        LT(2,:) = std(T(tt0,3));     %   standard deviation
        PHID = mean(T(tt0,2));       %   orientation
    case 1                %   There are two scales of dunes.
        LT1 = mean(T(tt0(1:t0),3));
        PHID1 = mean(T(tt0(1:t0),2));
        LTstd1 = std(T(tt0(1:t0),3));
        LT2 = mean(T(tt0(t0+1:end),3));
        PHID2 = mean(T(tt0(t0+1:end),2));
        LTstd2 = std(T(tt0(t0+1:end),3));
        LT(1,:) = [LT1 LT2];
        LT(2,:) = [LTstd1 LTstd2];
        PHID = (PHID1+PHID2)/2;
    case 2                %    There are three scales of dunes.
        LT1 = mean(T(tt0(1:t0(1)),3));
        PHID1 = mean(T(tt0(1:t0(1)),2));
        LTstd1 = std(T(tt0(1:t0(1)),3));
        LT2 = mean(T(tt0(t0(1)+1:t0(2)),3));
        PHID2 = mean(T(tt0(t0(1)+1:t0(2)),2));
        LTstd2 = std(T(tt0(t0(1)+1:t0(2)),3));
        LT3 = mean(T(tt0(t0(2)+1:end),3));
        PHID3 = mean(T(tt0(t0(2)+1:end),2));
        LTstd3 = std(T(tt0(t0(2)+1:end),3));
        LT(1,:) = [LT1 LT2 LT3];
        LT(2,:) = [LTstd1 LTstd2 LTstd3];
        PHID = (PHID1+PHID2+PHID3)/3;
end       

NL = size(LT,2);
WL = [num2str(LT(1,1),'%.0f') ' m'];
for i=2:NL
    WL = [WL ', ' num2str(LT(1,i),'%.0f') ' m'];
end    
disp(['There are ' num2str(NL) ' wavelength(s) of interest, they are ' WL]);

Kprof = 1/tand(PHID);     
if PHID <= 90
    Xprof = X(1,:);
    Yprof = Kprof*(Xprof-Xmin) + Ymin;
else 
    Yprof = Y(:,1)';
    Xprof = (Yprof-Ymax)/Kprof + Xmin;
end
%  One profile perpendicular to the dune crest is selected, according to the dune orientation.
Zprof = interp2(X,Y,Z,Xprof,Yprof);     
Xprof0 = Xprof(find(isnan(Zprof) ~= 1));
Yprof0 = Yprof(find(isnan(Zprof) ~= 1));
Zprof0 = Zprof(find(isnan(Zprof) ~= 1));
Xprof2 = Xprof0 - min(Xprof0);
Yprof2 = Yprof0 - min(Yprof0);
L = sqrt(Xprof2.^2 + Yprof2.^2);
D = Zprof0;
L0 = [0:1:floor(max(L))];
D0 = interp1(L,D,L0);
L0 = L0(find(isnan(D0) ~= 1));          %     the length of the profile
D0 = D0(find(isnan(D0) ~= 1));          %     the depth of the profile

%  Analyzing the periodic characteristics of the profile with Wavelet transform
[wltPower, fourierPeriod, coneOfInfluence] = ...
    waveletTransform (6, D0, 1, 0.05);
[wltSignif, globalSignif, globalWs,lag1,scaleAveSignif,scaleAvg] =...
    statisticsWlt(6, D0, 1, 0.05, 0.95,wltPower);
[peaksCross,xPeakspos,yPeakspos] = ...
    getWltpeaks(globalWs,globalSignif,fourierPeriod);

signalData = D0;
signalDataX = L0;
titleFont = 12;
labelFont = 12;
axisFont = 12;
sampleName = 'Length';
signalName = 'Depth';
sampleUnit = 'm';
signalUnit = 'm';
signalVariance = std(signalData)^2;
plotXRange = [min(signalDataX), max(signalDataX)];
Yticks = 2.^(fix(log2(min(fourierPeriod))):fix(log2(max(fourierPeriod))));
FP = interp1(fourierPeriod,globalWs,LT(1,:));

figure;                   %   plot the profile 
set(gcf,'unit','centimeters','position',[5 2 14 6]);
subplot('position',[0.15 0.23 0.8 0.75]);
plot(L0,D0,'b');
xlabel([sampleName '(' sampleUnit ')'],'fontsize',12);
ylabel('Depth (m)','fontsize',12);
XLimits =  [min(L0) max(L0)];
set(gca, 'XLim', XLimits);
grid on;
set(gca, 'fontsize',axisFont);
axis tight;
y = get(gca,'YTick');
set(gca,'YTickLabel', cellstr(num2str(y(:), '%.1f')) );
