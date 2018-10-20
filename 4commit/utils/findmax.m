function [locs,v]=findmax(X, dim, list, interval,method)
%            [locs,v]=findmax(X, dim, list, interval,method)
%It sorts observations in X in descending order over  dimension(s) 'dim',
%and returns the full list ('list'=1) or just the first, i.e. the maximum ('list=0').    [list=0]
%Default is dim=2, unless size(X,2)=1, then the longest dimension is taken.
%If dim=[] than sorting is applied over all observations/absolute maximum is returned.
%Sorting is done along  the specified dim and run indipendently for each
%other dimenions,  unless  numel(dim)>1.
%If numel(dim)>1 the maximum along dim(1:end-1) is iteratively taken and
%then sorting is applied along dim(end). For dim(1) a different 'method'
%can be applied if provided (e.g. mean,std,kurtosis).
% 'Interval' narrows observations to the specified subset of samples in the smallest
%size(X)==numel(interval).
%           [locs,v]=findmax(X, dim, list, interval,method)

if ~nargin
    help findmax
    return
end

dims=[size(X) 1 1];

if ~exist('list','var') || isempty(list) 
    list=0;
end

if ~exist('dim','var')
    if dims(2)==1
        [~,dim]=max(dims);
    else
        dim=2;
    end
end


%get rid of NaN related issues
absmin=floor(min(X(:))-2);
X(isnan(X))=absmin;
% and extend it to observations outside the interval
if exist('interval','var') && ~isempty(interval) 
    interval=~interval;
    target=find(dims==numel(interval),1);
    eval([' X('  repmat( ':,',[1 target-1])   'interval'  repmat(',:',[1 numel(dims)-target])   ')=absmin;']);
end

absmax=0;
if isempty(dim)
    X=X(:);
    dim=1;
    absmax=1;
end

if numel(dim)>1
if exist('method','var') && ~isempty(method) 
    try
            eval(['X='   method  '(X,[],dim(1));']);
    catch
            eval(['X='   method  '(X,dim(1));'])
    end
    dim(1)=[];
end    
    for n=1:numel(dim)-1
        X=max(X,[],dim(n));
    end
    dim=dim(end);
end



if ~list
    [v,locs]=max(X,[],dim);    
else
    [v,locs]=sort(X,dim,'descend');
end



locs(v==absmin)=NaN;
v(v==absmin)=NaN;
if absmax
    [dim1,dim2,dim3,dim4,dim5]=ind2sub(dims,locs);
    locs=[dim1 dim2 dim3 dim4 dim5];
    locs=locs(:,1:find(dims>1,1,'last'));
end

locs=squeeze(locs);
v=squeeze(v);


