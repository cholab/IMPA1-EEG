function [btr,bel,maxZ,mzZ,centZ]=findbad(X,chanlocs,method,TH,sr)
%      [btr,bel,Z,mzZ,cZ]=findbad(X,chanlocs,method,TH,sr)
%It finds bad channels/trials in data X according to "method":
%  -Z    distance from template approach
%  -PCA  correlation based approach
%   and FULL runs both and finds agreement

if~nargin
    help findbad
    return
end

if ~exist('TH','var')||isempty(TH)
    TH=3;
end
if ~exist('sr','var')||isempty(sr)
    sr=250;
end
if ~exist('method','var')||isempty(method)
    method='FULL';
end

method=upper(method);



if method(1)=='Z'||method(1)=='F'
%"ABSOLUTE DISTANCE" from template...
%So that differences are captured in the direction of "worse" eeg quality
%% Rereferencing to improve comparison with template
M=eye(size(X,1))-ones(size(X,1))*1/size(X,1);
X=reshape(M*X(:,:),size(X,1),size(X,2),[]);

       [maxZ,~,mzZ,centZ]=checkeeg(X,sr);
       
if TH<0
   TH=abs(TH);
   Zloc=centZ; %centered version
else
   Zloc=maxZ;
end
Zth=Zloc>TH;

%WORSE trials/channels approach
%trials during which more than 30% channels got worse (+ Z value check)
TR=mean(zscore(Zloc,1,2)>TH  & Zth,1)<=.3; 
%channels that are worse than the others for over half trials (+ Z value check)
CH=mean(zscore(Zloc(:,TR,:),1,1)>TH & Zth(:,TR) ,2)<=.5; 
%now deviations should be more local in time/space...

%BAD trials/channels approach
%bad trials (i.e. where the number of bad channel is greater than some expected value)
btr=mean(Zth(CH,TR),1); %mean bad segments per trial
btr(btr==1)=1-eps;
btr=zscore(atanh(btr))>TH; %the number of bad channels per trial is greater than TH
TR(TR)=~btr;
%bad electrodes
bel=mean(Zth(CH,TR),2); %mean bad segments per channel
bel=bel-min(bel); %the best channel is all we have, this centers the observations in really bad cases
[max_el,worst_el]=max(bel);
if max_el>.2  %finding channels bad for more than 1/5 trials
%I use some data segmentation to maximize dataset specific distribution
try
warning('off','stats:kmeans:EmptyCluster');
warning('off','stats:kmeans:EmptyClusterInBatchUpdate');
cl=kmeans(bel,4,'replicates',100);
bel=cl==cl(worst_el);
CH(CH)=~bel;
end
end
btr1=~TR;
bel1=~CH;
end







if method(1)=='P'||method(1)=='F'
V=sum(squeeze(std(X,[],2)));  %exclude biasing trials
ref=quantile(V,[.05 .95]); %robust zscoring
i=V>ref(1) & V < ref(2);
V=(V-mean(V(i)))/std(V(i));
i=V<4;
[D,T2,iW,A]=eldistance(X(:,:,i),chanlocs);
%Bad channels are uncorrelated 
bel2=(zscore(nanmean(D)))'>abs(TH);
%Bad trials are distant observations (in the multivariate sense)
T2=squeeze(reshape(T2',size(X,2),[]));
T2=max(T2);
btr2a=zscore(T2)>abs(TH);
%Searching for local principal component activities
[V,R]=explainedvar(A,iW,3);
[~,E]=explainedvar(A,iW);

tAV=zscore(V,1,2);
tAR=zscore(R,1,2);
c=squareform(pdist(R'));
%bad local activity criteria...  
%only main components
ii=cumsum(E)<.9; 
if sum(ii)==1
btr2b1=tAV(ii,:)>10;    
else
btr2b1=any(tAV(ii,:)>10);
end
%all components
btr2b2=mean(tAV>5)>0.05;
btr2c=any(tAR>20) | mean(tAR>5)>0.05;
btr2d=zscore(mean(c))>4;

% btr2=mean(maxZ(~bCH2,~bTR2)>5)>0.02; %in presence of bad segments.
btr2=btr2a | btr2b1 | btr2b2 | btr2c | btr2d;

i0=i;
i( i0)=btr2;
i(~i0)=1;

btr2=i;
end







if method(1)=='Z'
    btr=btr1; bel=bel1;
elseif method(1)=='P'
    btr=btr2; bel=bel2;
    maxZ=[];,mzZ=[];centZ=[];
elseif method(1)=='F'
    btr=btr1 & btr2;
    bel=bel1 & bel2;    
end






