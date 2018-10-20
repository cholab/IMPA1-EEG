function [V,R,i]=explainedvar(A,iW,dim)
%               [V,R,i]=explainedvar(A,iW,dim)
% Given component strengths 'A' (comp x tp x tr) and mixing matrix 'iW' (channels x comp), it returns variance
% accounted by each component 'V' (and ratio of total variance 'R'). 
% V is computed according to 'dim'[0]:
%   dim = 0 total variance (comp x1)
%   dim = 2 time-variance trajectory (comp x tp x tr)
%   dim = 3 trialwise (comp x tr)
% 'i' is an index that order components in decreasing order of explained
% variability, with the same dimensions as V.

if ~nargin
    help explainedvar
    return
end
    
if ~exist('dim','var') || isempty(dim)
    dim=0;
end

SW=sum(iW.^2,1);
[c,tp,tr,d4]=size(A);
switch dim
    case 0
        df=(c*tp*tr)-1;
        V = (SW.*sum((A(:,:)').^2,1))';  
    case 2
        df=c-1; 
        V=zeros(c,tp*tr);
        for n=1:tp*tr*d4
        V(:,n) = SW.*sum((A(:,n)').^2,1);  
        end
        if tr*d4>1
            V=reshape(V,[c tp tr d4]);
        end
    case 3
        df=c*tp-1; 
        V=zeros(c,tr);
        for n=1:tr
        V(:,n) = (SW.*sum((A(:,:,n)').^2,1))';  
        end    
    case 4
        df=c-1; 
        V=zeros(c,tr,d4);
        for n=1:tr
        for m=1:d4    
        V(:,n,m) = (SW.*sum((A(:,:,n,m)').^2,1))';  
        end
        end        
end
V=V/df;


if nargout>1    
R=V./repmat(sum(V,1),[size(V,1) 1 1 1]);
end

if nargout>2   
    [~,i]=sort(V(:,:),'descend');
    i=reshape(i,size(V));
end







