function   [OUT,LBL,centralFq]=f2b(IN,fq,dim,peakvalue)
%                [OUT,LBL,centralFq]=f2b(IN,fq,dim,peakvalue)
%It collapses frequency-wise values  in IN, within the following band of interests:
%
%           Delta,Theta,Alfa,Beta1,Beta2,Gamma1,Gamma2*  <  4,8,13,18,30,48,"100"
%
%'fq'   is the vector of input frequencies [1:size(IN,2), as in the standard wavelet output]
%'dim'   specifies the dimension in IN that matches 'fq' for ambiguous cases  [first matching dimension]
%'IN'   can have ANY number of dimensions, ANY dimension can be the fq dimension
%'peakvalue'   calls for max value within each band instead of mean value [0]
%
%'OUT' has the same dimensions as IN but the fq>bands dimension 
%'LBL' provides labels for the used bands
%          Use:            LBL=f2b      to get labels at any time
%'centralFq' is the central value ofr the used bands (for numerical labelling)
%
% *   Note that 57<fq<63 is not included and upper value is actually 80Hz ,
% in order to avoid issues related to line noise and 100Hz lowpass filter
% response.





LBL{1}='DELTA <4';
LBL{2}='THETA <8';
LBL{3}='ALFA <13';
LBL{4}='BETA1 <18';
LBL{5}='BETA2 <30';
LBL{6}='GAMMA1 <48';
LBL{7}='GAMMA2';

if ~nargin
if ~nargout
    help f2b
    return
end
    OUT=LBL;
else

    
    
dims=size(IN);    

if ~exist('fq','var')||isempty(fq)
    fq=1:dims(1);
end
nfq=numel(fq);


bandup=[4 8 14 18 30 48 100]; %bands upper values

bandup=[0 bandup 1000];
blim_min=find(bandup<min(fq), 1,'last');
if isempty(blim_min)
   blim_min=1; 
end
if max(fq)>bandup(end)
    blim_max=9;
else
    blim_max=find(bandup>max(fq), 1,'first');
end
bandup=bandup(blim_min:blim_max);

[~,B]=histc(fq,bandup);
B=B+1;

B(fq>57 & fq<63)=NaN; B(fq>80)=NaN;
band=unique(B(B>0));
LBL=LBL(blim_min-2+band);

mB=min(B(B>0));
B=B-mB+1;
mB=max(B);
B(find(B==mB,1,'last')+1:end)=NaN;


if ~exist('dim','var') || isempty(dim)
    dim=size(IN);
    dim=find(dim==nfq,1);
end
if ~exist('peakvalue','var') || isempty(peakvalue)
    peakvalue=false;
end



if numel(dim)>1
    error('dimensions interpration is ambigous! (add ''dim'' argument)')    
elseif isempty(dim)
    error('fq vectors does not match input dimensions!')
else
   
dimsB=dims;
dimsB(dim)=mB;
OUT=zeros(dimsB);

str=['OUT(' repmat(':,',1,dim-1)   'n'   repmat(',:',1,numel(dims)-dim ) ')='];
if peakvalue
str=[str 'max(IN(' repmat(':,',1,dim-1)   'B==n'   repmat(',:',1,numel(dims)-dim ) '),[],'  num2str(dim) ');'];
else
str=[str 'nanmean(IN(' repmat(':,',1,dim-1)   'B==n'   repmat(',:',1,numel(dims)-dim ) '),'  num2str(dim) ');'];
end

for n=1:mB
    eval(str)
end
if nargout>2
    bandup2=[1 bandup(1:end-1)];
    centralFq=mean([bandup2;bandup]);
    centralFq=centralFq(band);
end
end
end


