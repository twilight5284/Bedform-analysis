% % Search value 'v' in sorted vector 'x' and find index and value
% % with respect to vector x that is equal or closest to 'v'.
% % 
% % If more than one value is equal then anyone can be returned
% % (this is property of binary search).
% % If more than one value is closest then first occurred is returned
% % (this is property of linear search).
% % 
% %
% % INPUT:
% % x: vector of numeric values,
% % x should already be sorted in ascending order
% %    (e.g. 2,7,20,...120)
% % v: numeric value to be search in x
% %  
% % OUTPUT:
% % i: index of v with respect to x. 
% % cv: value that is equal or closest to v in x 
       
function [i,cv] = searchclosest(x,v)

i=[];
lo=1;
hi=length(x);

% % Phase 1: Binary Search
while hi-lo>1
    mid = round((hi + lo)/2);    
    diff = x(mid)-v;
    if diff==0
        i=mid;
        cv=v;
        return
    elseif diff<0     % x(mid) < v
        lo=mid;
    else              % x(mid) > v
        hi=mid;			
    end
end

if ( abs(x(lo)-v) < abs(x(hi)-v) )
    i = lo;
    cv = x(lo);
else
    i = hi;
    cv = x(hi);
end
