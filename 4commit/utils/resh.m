function X2=resh(X,dims,reshaping_i)
%        X2=resh(X,dims,reshaping_i)
if ~nargin
    help resh
   return
end
 
d=size(X);
if ~exist('reshaping_i','var') %NO reshaping indices
if ~exist('dims','var')||isempty(dims)|| numel(dims)==1 %assuming we want to move from/to 2D/3D
    
if numel(d) ==3 
    X2=reshape(X,[d(1) d(2)*d(3)]);    
elseif numel(d)==2
    X2=reshape(X,[d(1) dims d(2)/dims]);    %one single dims value is taken to code  d2 in the reshaped 3D
end  

else
if numel(d)==3 && numel(dims)==2 %if 3D to 2D
if dims(1)==d(1)*d(3) %if 3D is aligned along 1st
X2=zeros(dims);
for n2=1:d(2)    
X2(:,n2)=reshape(squeeze(X(:,n2,:)),[d(1)*d(3) 1]);
end
elseif dims(1)==d1 && dims(2)==d(2)*d(3)
    X2=reshape(X,dims);
end

elseif numel(d)==2 && numel(dims)==3 %if 2D to 3D
if d(1)==dims(1)*dims(3) %if 3D is aligned along 1st
X2=zeros(dims);    
for n2=1:dims(3)    
X2(:,:,n2)=X((n2-1)*dims(1)+(1:dims(1)),:);
end
elseif dims(1)==d(1) && dims(2)*dims(3)==d(2)
    X2=reshape(X,dims);
end
end
end


else %reshaping indices
    if numel(d)==3 %from 3D to 2D for selected timepoints
        X2=X(:,reshaping_i);  
    elseif numel(d)==2 %from 2D to 3D 
        X2=X(:,reshaping_i);     X2=reshape(X2,d(1), dims, []);        
    end
end


