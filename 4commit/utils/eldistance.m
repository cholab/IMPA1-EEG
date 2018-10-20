function [D,T2,L,PC]=eldistance(X,chanlocs)
%    [D,T2,L,PC]=eldistance(X,chanlocs)
if ~nargin
    help eldistance
    return
end

% fprintf('Running exploratory PCA...')
warning('off','stats:princomp:colRankDefX')
[L,PC,E,T2]=princomp(X(:,:)');
% fprintf(' done.\n')

E=E/sum(E);
Le=L.*repmat(E',[size(L,1) 1]); %note that it emphasizes high variability components
distL=squareform(pdist(Le)); %weighted loadings distance
if exist('chanlocs','var') && ~isempty(chanlocs)
fprintf('Computing  geometry normalized multivariate distances... ')
xyz=[struct2mat(chanlocs,'X')  struct2mat(chanlocs,'Y') struct2mat(chanlocs,'Z')];
distG=squareform(pdist(xyz)); %geometric distance
distG=distG/max(distG(:));
distL=distL/max(distL(:));
D=distL.*(1-distG); %normalized between electrodes distance
D(eye(size(D))==1)=NaN;
fprintf(' done.\n')
else
% fprintf('Computing  NOT geometry normalized multivariate distances... ')
D=distL/max(distL(:));
D(eye(size(D))==1)=NaN;    
% fprintf(' done.\n')
end


if nargout>2
    PC=reshape(PC',[],size(X,2),size(X,3)); 
end
