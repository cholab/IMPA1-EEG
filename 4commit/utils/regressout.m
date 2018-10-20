function X=regressout(X,M)
%Compute residuals from multiple linear regression of M on X, i.e. it
%"regresses out" M from X:
%          Xr=regressout(X,M)
%
if ~nargin
    help regressout
    return
end
   
    d3=size(X,3);    
    X=X(:,:);
    M=M(:,:);  
    M=M(~any(isnan(M),2),:);
    for n=1:size(X,1)
    [~,~,X(n,:)]=regress(X(n,:)',M');
    end
    
if d3>1
   X=reshape(X,n,[],d3);
end
    