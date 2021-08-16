function [CT,BP] = CT_BP(bedformProfile,Y,Z,Z0,L,resolution)
%% bedform location and geometries calculation
% bedform crests and troughs are find in the bedform profile, 
% and then bedform geometry parameters are calculated.
%
% INPUTS
%     bedformProfile - bedfrom discrimination outputs
%     Y - y coordinate of the profile
%     Z - the orignal depth
%     Z0 - the smoothed depth
%     L - the wavelength of interest
%     resolution - the resolution of the profile
%       
% OUTPUTS
%     CT - the location of the crests and troughs
%     BP - the geometry parameters of all dunes 
%
% Author:  Li Wang
% Email:   twilight528400@hotmail.com

%%
CT = [];
BP = [];
ZZ = Z-Z0;
for l = 1:length(bedformProfile)
  %load bedformdata
  bedformdata = [];
  bedformdata = bedformProfile(l);
  
  %load xyz 
  data = [];
  pos = [];
  loc = [];
  indexY=[];
  data = bedformdata.n23s;
  data0 = bedformdata.n23s + bedformdata.n33;
  pos = bedformdata.signalAbcise;
  loc = bedformdata.Y;
  indexY = find(Y(:,1) == loc(1,1));
  L=bedformdata.target(2);
  
  %find the zero point
  tc = [];
  pc = [];
  tc0=[];
  pc0=[];
  peaktroughs = [];
  locations = [];
  [tc0,pc0] = crossing(data,pos);
  [tc,pc] = zpf(tc0,pc0,L);
  
%   % use data for filter peaktroughs
%   peaktroughs = nan(numel(tc)-1,1);
%   locations = nan(numel(tc)-1,1);
%   for i = 1:length(tc);
%     if i<length(tc) && i<length(pc) 
%         win = data0(tc(i):tc(i+1));
%         ps = pos(tc(i):tc(i+1));
%         if win(2)> win(1) % get maxima
%             [peaktroughs(i),locations(i)] = max(win);
%             locations(i) = ps(locations(i));
%         else % get minima
%             [peaktroughs(i),locations(i)] = min(win);
%             locations(i) = ps(locations(i));
%         end
%     end
%   end
%   PT=[];
%   PT = [locations,(loc(1,1)*ones(size(locations))),peaktroughs];
%  %set wave trough as start, keep num of peaktroughs an odd num.
%   if size(PT,1)>2
%     if PT(1,3)>PT(2,3)                           % start with a minmum
%       PT(1,:)=[];
%       tc(1)=[];
%       pc(1)=[];
%     end
%   end
  
  
  %find peaks and troughs
  Data = [];
  DeData = [];
  DeData = find(~isnan(ZZ(indexY,:)));
  Data = bedformdata.n23+bedformdata.n33+ZZ(indexY,DeData)';
  Data0 = detrend(Data);
  Data1 = bedformdata.n33+ZZ(indexY,DeData)';
  
  peaktroughs = nan(numel(tc)-1,1);
  peaktroughs0 = nan(numel(tc)-1,1);
  locations = nan(numel(tc)-1,1);
  locations0 = nan(numel(tc)-1,1);
  for i = 1:length(tc);
    if i<length(tc) && i<length(pc) % so we don't run beyond the end of the data.
        win = Data0(tc(i):tc(i+1));
        win0 = Data(tc(i):tc(i+1));
        ps = pos(tc(i):tc(i+1));
        if win(2)> win(1) % get maxima
            [peaktroughs(i),locations(i)] = max(win);
            locations0(i) = ps(locations(i));
            peaktroughs0(i) = win0(locations(i));
        else % get minima
            [peaktroughs(i),locations(i)] = min(win);
            locations0(i) = ps(locations(i));
            peaktroughs0(i) = win0(locations(i));
        end
    end
  end
  
  %save xyz of peaktroughs
  PT0=[];
  PT0 = [locations0,(loc(1,1)*ones(size(locations))),peaktroughs0];
  DPT=[];
  DPT=diff(PT0(:,3));
  DPT=DPT./abs(DPT);
  dpt=[];
  dpt=diff(DPT);
  PT0(find(dpt==0)+1,:)=[];
  PTf=[];
  PTf(:,1)=diff(PT0(:,1));
  PTf(:,2)=abs(diff(PT0(:,3)));
  Numf=find(PTf(:,1)<mean(PTf(:,1))/5 & PTf(:,2)<mean(PTf(:,2))/5);
  PT0(Numf,:)=[];
  DPT=[];
  DPT=diff(PT0(:,3));
  DPT=DPT./abs(DPT);
  dpt=[];
  dpt=diff(DPT);
  PT0(find(dpt==0)+1,:)=[];
  %calculate dune morphological parameters.
  BPQ = [];
  if size(PT0,1) > 2
    [BPQ] = bqp(PT0,L,Data,pos,resolution);
  end
  %save peaktroughs and locations;
  if size(PT0,1) > 0
    CT(((size(CT,1)+1):(size(CT,1)+size(PT0,1))),:) = PT0;
  end
  %save dune locations and their morphological parameters.
  if size(BPQ,1) > 0
    BP((size(BP,1)+1):(size(BP,1)+size(BPQ,1)),:) = BPQ;
  end
end

   
    