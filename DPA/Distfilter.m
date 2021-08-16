function [BPf] = Distfilter(BP0,L)
 DistLag = max(2,floor(L));
 DuneCrest = BP0(:,1:2);
 DistNum = [];
 for i = 1:size(BP0,1)
    Dist = [];
    DC0 = DuneCrest(i,:);
    Dist = sqrt((DuneCrest(:,1)-DC0(1)).^2+(DuneCrest(:,2)-DC0(2)).^2);
    DistNum(i,1) = length(find(Dist < 5*DistLag));
 end  
 DistFil = find(DistNum < 10);
 BPf=BP0;
 BPf(DistFil,:)=[];