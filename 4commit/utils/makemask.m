function MASK=makemask(iy,ix,X,method)
%              MASK=makemask(iy,ix,X)
%Create a mask:
%-if X is NOT given
%   size(MASK)=[length(iy) length(ix)], filled with NaNs iy/ix rows/colums
%-if X is given, 2 possible methods (applied/inferred) [2,for ambiguous cases]
%   1 - if size(X)=[sum(iy) sum(ix)], iy/ix specify locations in a
%   equal/bigger MASK that will be filled with X 
%       (i.e. MASK iy/ix==1 will be X, iy/ix==0 locations will be NaNs)
%   2 - if size(X)=[length(iy) length(ix)], then iy/ix specify locations
%   within X that will be filled with NaNs 
%       (i.e. X=MASK, iy/ix==0 locations will be NaNs)



if ~nargin
    help makemask
    return
end

iy=logical(iy(:));
ix=logical(ix(:));
if exist('X','var') && ~isempty(X)  
    
    m1ok=all([size(X,1) size(X,2)]==[sum(iy) sum(ix)]);
    m2ok=all([size(X,1) size(X,2)]==[length(iy) length(ix)]);
    
    if ~exist('method','var')||isempty(method)
    if m1ok && m2ok
        display('** ambiguous ix/iy interpretation!  **')
        method=2;        
    elseif  m1ok && ~m2ok
        method=1;
    elseif ~m1ok &&  m2ok
        method=2;
    else
        method=0;
    end
    end
    
    if (method==1 && ~m1ok) || (method==2 && ~m2ok)
        method=0;
    end            
    
    if method==1
%         display('method 1')
    MASK=nan(length(iy),length(ix), size(X,3));
    MASK(iy,ix,:)=X;    
    elseif method==2
%         display('method 2')
    MASK=X;
    MASK(:,ix,:)=NaN;
    MASK(iy,:,:)=NaN;    
    else
        %mat2str[size(X,1) size(X,2) sum(iy) sum(ix) length(iy) length(ix)]
        error('ix/iy - X inconsistency: not sure what to do with X!')
        MASK=[];
    end
else
    MASK=ones(length(iy),length(ix));
    MASK(:,ix)=NaN;
    MASK(iy,:)=NaN;
end

    




