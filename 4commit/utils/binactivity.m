function bA=binactivity(A,wlength,threshold)
%          bA=binactivity(A,wlength,threshold)
if ~nargin
  help binactivity
  return
end


tp=size(A,2);
if size(A,3)==1
tr=size(A,1);
na=1;
else
na=size(A,1);    
tr=size(A,3);
end

if ~exist('wlength','var') || isempty(wlength)
%default: window of 100ms (assuming sr=250Hz)
wlength=.1*250;
end
if ~exist('threshold','var') || isempty(threshold)
    thresholdmode=0;    
else
    thresholdmode=1;
end




if wlength==tp
rshp=0;
nw=1;
else
rshp=1;
wcenters=linspace(0,tp,ceil(tp/wlength));
wcenters=wcenters(1:end-1)+(wcenters(2)-wcenters(1))/2;
wcenters=[ceil(wcenters(1)) round(wcenters(2:end-1)) fix(wcenters(end))];

nw=numel(wcenters);
w=1:fix(wlength/2); w=[fliplr(w) 0 w];
i=repmat(wcenters,[tr 1])+repmat(tp*(0:tr-1)',[1 nw]);
i=reshape(i',numel(i),[]);
i=repmat(i,[1 numel(w)])+repmat(w,[numel(i) 1]);
end


% bA=zeros(tr,nw,na);
bA=zeros(na,nw*tr);

fprintf('%s','Binning activities ')
if na<10
    na2=1:na;
else
    na2=round(linspace(1,na,10));
end
for n=na2;     fprintf('%s', '.'); end
for n=1:na 
if any(n==na2); fprintf('\b\b%s', ' '); end    
if na>1
a=squeeze(A(n,:,:))';
else
a=A;
end
 
 
if ~thresholdmode %transform data in some meanigful way to detect activity levels
      a=(detrend(a')); 
%    a=wfilter(a(:),5);
      a=zscore(detrend(a(:)));
%    [p,l]=findpeaks(abs(a),'minpeakheight',1)
%    a(abs(a)<1)=0;
   a=reshape(a,tp,[])';
end
 

if rshp %needs to be reshaped
%make it continuous
a=a';
a=a(:);
%apply i
a=a(i);
end


%NOW a is (nw,wlength)
%-whatever the number of input activities
%-whatever window length 

if thresholdmode  
%for th let's mark if any occurance is present inside the window
a=any(a>th,2);
a=reshape(a,nw,[])';
else
%are we interested in std or peaks?    
a=[a zeros(size(a))];%this trick should help to detect both...
a=zscore(std(a,1,2));
a=reshape(a,nw,[])';
end

bA(n,:)=a(:)';
end
bA=reshape(bA,na,nw,tr);
fprintf(': done.\n')

