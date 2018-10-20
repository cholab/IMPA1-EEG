function [Z,ZZ,mzZ,cZ,bel,btr]=checkeeg(X,sr,thresh)
% Normalized  distance of eeg segments X from an artifact free
% template:
%       [Z,ZZ,mzZ,cZ,bel,btr]=checkeeg(X,sr,thresh)
%   X(:,timepoints,:)  eeg segments
%   'sr'               sample rate
%Distance is assessed between distribution moments and spectral properties,
%    Z   whatever is higher
%and: 
%     ZZ   exports distances for each measure
%    mzZ   mean of within measure distances
%     cZ   distances 'centered' to best channels and typical best trials 
% 
%Segments further than threshold*sd from expected value are "bad".
%Accordingly:
% -badelectrodes:  frequencies of bad trials per channel
% -badtrials: frequencies of bad channel per trial
%      [Z,ZZ,mzZ,cZ,bel,btr]=checkeeg(X,sr,thresh)

if ~nargin
    help checkeeg
return
end
    
 warning('off','stats:pca:ColRankDefX'); 
 warning('off','stats:princomp:colRankDefX');

%distributions
fprintf('%s','Computing distance from eeg parameters...')
fprintf('%s',' distribution...')
[Vv,~,Zv]=checksegments(X);
%spectral
if ~exist('sr','var') || isempty(sr)
    fprintf('%s',' assuming default sample rate: 250Hz...')    
    sr=250;
end 
fprintf('%s',' spectrum...')
[Vb,~,Zb]=checksegments(X,sr);
fprintf(' done. \n')
ZZ=cat(3,Zv,Zb);
%max distance
Z=max(ZZ,[],3);
%% ??
mzZ=ZZ-repmat(nanmean(reshape(ZZ,[],1,size(ZZ,3))), [size(ZZ,1) size(ZZ,2) 1]);
mzZ=meandim(mzZ./repmat(nanstd(reshape(ZZ,[],1,size(ZZ,3))), [size(ZZ,1) size(ZZ,2) 1]),3);

%Do we want "bad" segments?
if nargout>4  
fprintf('%s','Finding "bad" segments...')  
if ~exist('thresh','var') || isempty(thresh)
    fprintf('%s',' assuming default threshold: 3std...')    
    thresh=3;
end
Zth=Z>thresh;
btr=mean(Zth);
bel=mean(Zth,2);
fprintf(' done. \n')
end

if nargout>3  
   fprintf('%s','Centering to dataset...')        
   %trying to adjust for "dataset" specific differences
   iel=find(zscore(mean(squeeze(mean(zscore(ZZ),2)),2))<0); %most reliable electrodes 
   itr=find(zscore(mean(squeeze(mean(zscore(ZZ(iel,:,:),1,2),1)),2))<0); %most reliable trials
   
   [~,~,~,T2]=princomp(resh(X(iel,:,itr))');   
   T2=mean(squeeze(reshape(T2(:)',size(X,2),[])));
   i=T2<quantile(T2,.75); %~typical trials
   itr=itr(i);
    
   V=cat(3,Vv,Vb);   
   V=reshape(V(iel,itr,:),[],size(V,3));
   m=mean(V);
   s=std(V);
   cZ=cat(3,Vv,Vb);   
   for n=1:11
       cZ(:,:,n)=abs(cZ(:,:,n)-m(n))./s(n);
   end   
   fprintf(' done. \n')
   cZ=max(cZ,[],3)/2;
end



 
% [~,~,~, bel1]=pca(X.mztrans);
% bel2=max(zscore(abs(zscore(X.mztrans))),[],2);
% bel=max(zscore([bel1 bel2]),[],2);
% 
% [~,~,~, bel1]=pca(X.mztrans);
% bel2=max(zscore(abs(zscore(X.mztrans))),[],2);
% bel_f=max(zscore([bel1 bel2]),[],2);
% bel=max([bel bel_f],[],2);
% bel=bel>thresh;
 
