function [sparse_el,peripherical_el,S]=findsparsechannels(el,chanlocs,clsize)
%  [sparse_el,peripherical_el,S]=findsparsechannels(el,chanlocs,clsize)

if ~nargin
    help findsparsechannels
else
if ~islogical(el)
el_t=el;
el=false(size(chanlocs,2),1);
el(el_t(~isnan(el_t)))=true;
else
    el=el(:);
end

th=6;
if ~exist('clsize','var') ||isempty(clsize)
    clsize=1;
else
    if clsize>6
        th=clsize;
    end
end

xyz=[struct2mat(chanlocs,'X')  struct2mat(chanlocs,'Y') struct2mat(chanlocs,'Z')];  
d=squareform(pdist(xyz));   
peripherical_el=zscore(xyz(:,3))<-1.2;
[dn,n]=sort(d);
    
%     [~,n]=sort(n);
%     dn=dn./repmat(dn(5,:),[size(dn,1) 1]);
%     dn=dn(n);
%     for j = 1:size(dn,1), d(:,j) = dn(n(:,j),j); end
%     dth=1.1;
%     d=d<dth;
%     neighbors=sum(d(el,el));
%     sparse_el=sparse_el(neighbors==1);
S=n(1:3,:);
s=n(4:th,:);
s(dn(4:th,:)>min(dn(th+1,:)))=NaN;
S=[S;s];  
sparse_el=el & sum(ismember(S,find(el)))'>0 &  sum(ismember(S,find(el)))'<=clsize;
peripherical_el=el & peripherical_el;
if ~nargout
display( ['Sparse :',num2str(sum(sparse_el)),'/',num2str(sum(el))] ) 
end

if ~islogical(el)
    sparse_el=find(sparse_el);  
    peripherical_el=find(peripherical_el);  
end
end
