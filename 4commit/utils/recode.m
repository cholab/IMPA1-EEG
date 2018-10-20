function [out1,out2,out3,out4]=recode(x,y)
% It recodes any input 'x'  in numerical vector form 'X', or  recodes 'x'
% and 'y' with the same coding scheme, so that X and Y share the same
% mapping into a unique 'key':
%              [X,key]=recode(x,y)
%              [X,Y, key,CuXuY,Y]=recode(x,y)
%'CuXuY' labels elements in 'key' (and numerical items in X and Y) as
%follows:
%    1 for intersection
%    2 for elements unique to x 
%    3 for elements unique to y

if ~nargin
    help recode
    return
end


remap2=0; flip=0;
if exist('y','var') && ~isempty(y)
remap2=1;    
if ~strcmp(class(x), class(y))
    error('x-y class missmatch')
end
if isvector(x) && isvector(y)
    flip=1;
    x=x(:);
    y=y(:);
end 
if size(x,2)~=size(y,2)
if size(x,2)==size(y,1)
    warning('x-y size missmatch: y was flipped')                
    y=y';
else
    error('x-y size missmatch')                        
end
end
end
    
[X,uX]=grp2idx_all(x);
if ~remap2
out1=X;
if ischar(x) && iscell(uX)
    uX=char(uX);
end
out2=uX; %i.e. the key


else %x-y common remapping
[Y,uY]=grp2idx_all(y);
if isnumeric(uY)
    nancode=max([uX(:);uY(:)])+1;    
    uX(isnan(uX))=nancode;
    uY(isnan(uY))=nancode;        
end

if iscell(uX)
XY=intersect(uX,uY);
uuX=setdiff(uX,XY);
uuY=setdiff(uY,XY);
key=[XY; uuX; uuY];  CuXuY=[ones(size(XY,1),1); ones(size(uuX,1),1)*2; ones(size(uuY,1),1)*3];
[~,iX]=ismember(uX,key);
[~,iY]=ismember(uY,key);
else
XY=intersect(uX,uY,'rows');
uuX=setdiff(uX,XY,'rows');
uuY=setdiff(uY,XY,'rows');
key=[XY; uuX; uuY];  CuXuY=[ones(size(XY,1),1); ones(size(uuX,1),1)*2; ones(size(uuY,1),1)*3];
[~,iX]=ismember(uX,key,'rows');
[~,iY]=ismember(uY,key,'rows');    
end
X=iX(X);     Y=iY(Y);   
 
if ischar(x) && iscell(key)
    key=char(key);
elseif isnumeric(x)
    key(key==nancode)=NaN;
end
    if ~flip
    out1=X; out2=Y;
    else
    out1=X'; out2=Y';        
    end
    out3=key;
    out4=CuXuY;
end


function [NV,GL]=grp2idx_all(L)
%Note that this differs from gr2idx  as follows:
% -it handles columns  differently, not as multiple grouping variables 
%but as different specifiers of the same grouping strategy, 
%i.e. it is a grp2idx(L,'rows') which is not supported. 
% - it handles NaNs in numerical input assuming NaN=NaN, 
%i.e. unique groups  are formed for NaN occuring in the same column position
if isnumeric(L) || islogical(L)
     nancode=max(L(:))+1;    
     L(isnan(L))=nancode;
    [GL,~,NV]=unique(L,'rows');               
    GL(GL==nancode)=NaN;
else 
    if ischar(L) 
        if size(L,2)>1
        L=cellstr(L);
        end
    end
     [NV,GL]=grp2idx(L);    
end


