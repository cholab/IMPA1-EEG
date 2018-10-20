function [M,SD]=meandim(M,d,squeeze_sw,allbut,param)
%      [M,SD]=meandim(M,d,squeeze_sw,allbut,param)
% 'M'  input matrix
% 'd'  dimensions that will be collapsed by averaging discarding NaNs, but:
%      -if 'allbut' is 1 these are the dimensions that are NOT collapsed [0]
%      -if 'param' is 0 nanmedian is used instead [1]
% Output matrix will be squeezed unless squeeze_sw is ~=`1 [1]
%     squeeze_sw == 0 same dimensions as M but d=1
%     squeeze_sw == 2 same dimensions as M,i.e. it is tiled back to original size
%
%WARNINGS/2DOes:
%check for numel(d)>1: differentiate the mean of the means from the whole
%population mean in presence of missing values.
%SD support added, BUT:
% -not adjusted for numel(d)>1 
% -interquantile range used instead of std for non parametric approach
%      [M,SD]=meandim(M,d,squeeze_sw,allbut,param)

if nargin==0
    help meandim
    return
end


dim=size(M);
if ~exist('d','var')
    d=0;
end
if numel(dim)==numel(d) && ~any(d>1)   
   d=find(d);
end
if numel(dim) < max(d)    
    error('target dimensions exceeds target matrix!')
end          
if ~exist('squeeze_sw','var') || isempty(squeeze_sw)
    squeeze_sw=1;
end
if ~exist('allbut','var') || isempty(allbut)
    allbut=0;
end
if ~exist('param','var') || isempty(param)
    param=true;
end
if allbut
   d2=1:numel(dim);
   d2(d)=[];
   d=d2;
end
try
M(isinf(M))=NaN;
end

if nargout<2
if param
if sum(d)==0 || (numel(dim)==numel(d) && sum(d)==sum(1:max(d))) %collapse all
  M=nanmean(M(:));
elseif numel(d)==1 %collapse 1 dimension
  M=nanmean(M,d);      
else
for n=1:numel(d)
    M=nanmean(M,d(n));
end
end

else %use nanmedian
if sum(d)==0 || (numel(dim)==numel(d) && sum(d)==sum(1:max(d))) 
  M=nanmedian(M(:));
elseif numel(d)==1
  M=nanmedian(M,d);      
else
for n=1:numel(d)
    M=nanmedian(M,d(n));
end
end    
end

if squeeze_sw==1
    M=squeeze(M);
elseif squeeze_sw==2
    M=repmat(M,[dim./size(M)]);
end


else %adding std support    
if param
if sum(d)==0 || (numel(dim)==numel(d) && sum(d)==sum(1:max(d))) 
  SD=nanstd(M(:));    
  M=nanmean(M(:));
elseif numel(d)==1
  SD=nanstd(M,1,d);          
  M=nanmean(M,d);      
else
for n=1:numel(d)
    SD=nanstd(SD,1,d(n));        
    M=nanmean(M,d(n));
end
end
else %use nanmedian
if sum(d)==0 || (numel(dim)==numel(d) && sum(d)==sum(1:max(d))) 
  SD=mean(quantile(M(:),[.25 .75]));    
  M=nanmedian(M(:));
elseif numel(d)==1
  SD=mean(quantile(M(:),[.25 .75]),d);        
  M=nanmedian(M,d);      
else
for n=1:numel(d)
    SD=mean(quantile(M,[.25 .75]),d(n));        
    M=nanmedian(M,d(n));    
end
end    
end

if squeeze_sw==1
    M=squeeze(M);
    SD=squeeze(SD);    
elseif squeeze_sw==2
    M=repmat(M,[dim./size(M)]);
    SD=repmat(SD,[dim./size(SD)]);    
end
end


