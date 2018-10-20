function [Y,y,y12]=reshaper(X,directed)
%Matrix-to-pairs/pairs-to-matrix converter:
%           [Y,y3,y12]=reshaper(X,directed)
% If X:
%- is recognized as a matrix (i.e. it is NxN), it is converted into a
% list of pairs indexing the relative value. 
%- is recognized as a list of pairs (i.e. it is Nx3) it is
%converted into a matrix. 
%- is a vector  (i.e. N x 1) or a single number, than a list of pairs is created.
%If X cannot be unambiguosly recognized, an error is produced.  
%Unless 'directed' is 1, pairs are returned in the most economic form (i.e.
%each pair is taken once).
% Optionally also y and y12 are returned. They are meaningful only for
% matrix2pairs conversion: Y is split into y, the vector of values, and
% y12, the corresponding indeces, i.e. Y=[y y12].


if ~nargin
    help reshaper
else

if ~exist('directed','var') || isempty(directed)
    directed=0;
end

for nd=1:size(X,3)
x=X(:,:,nd);    
D=size(x);



if ~diff(D) && D(1)~=1 && D(1)~=3
  %display('MATRIX2PAIRS')
   L=D(1)^2-D(1);
   y(1:L,1:3)=NaN;
   s=0;
   for n=1:D(1) 
   for m=1:D(1) 
       if n~=m
       s=s+1;
       y(s,:)=[n m x(n,m)];
       end
   end
   end
%    y=y(~isnan(y(:,3)),:);
   if ~directed
   [~,i] = unique(sort(y(:,1:2),2), 'rows');
   y=y(i,:);
   end
    
    
elseif sum(D==3)==1 || D(1)==3
      %display('PAIRS2MATRIX ')   
    if D(1)<D(2) && D(1)>1
        x=x';
    end
    nn=max([x(:,1); x(:,2)]);
    y=zeros(nn,nn);
    x(~(x(:,1)-x(:,2)),:)=[];
    b = unique(sort(x(:,1:2),2), 'rows');
    if ~sum(size(b,1)==D)
        for n=1:size(x,1)
            y(x(n,1),x(n,2))=x(n,3);
        end
    else
        for n=1:size(x,1)
            y(x(n,1),x(n,2))=x(n,3);
            y(x(n,2),x(n,1))=x(n,3);
        end
    end
    y(logical(eye(size(y))))=NaN;
    
elseif sum(D==1)==1 && directed>=max(x)   
   D=max(D);
   y=zeros(directed);
   for n=1:D
   for m=1:D       
       y(x(n),x(m))=1;
   end
   end
    for ii=1:directed
        y(ii,ii)=NaN;
    end
elseif sum(D==1)>0
%display('just pairs')  
   if numel(x)==1
       x=1:x;
   end
   y(1:numel(x)^2,1:2)=NaN;
   s=0;
   for n=1:numel(x)
   for m=1:numel(x)
       s=s+1;
       y(s,:)=[n m];
   end
   end
   y((diff(y,1,2)==0),:)=[];
   if ~directed
   [b,i] = unique(sort(y(:,1:2),2), 'rows');
   y=y(i,:);
   end
   y=x(y);
    
else
    help reshaper
    error('Input cannot be recognized!')
end

if (isempty(y) && nd>1) || nd>1 && length(y)~= size(Y,1)
    Y(1:nn,1:nn,nd)=NaN;
else
    Y(:,:,nd)=y;
end

end



if diff(size(Y(:,:,1))) && size(y,2)==3
    y=squeeze(Y(:,3,:)); y12=Y(:,1:2,1);
else
    y=[]; y12=[];
end

end