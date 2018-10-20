function [V,X,Z]=checksegments(M,sr,verbose,freq)
%  [V,X,Z]=checksegments(M)     -> quick distribution check            (1)
%  [V,X,Z]=checksegments(M,sr) -> quick relative spectrum check  (2)
% (1)   STD  RANGE SKEWNESS KURTOSIS
% (2)   DELTA  THETA ALFA BETA1 BETA2 GAMMA1 GAMMA2
%Outputs:
% V                  raw values
% X .mztrans   average absolute between-electrodes distance
%    .mzlong    average absolute between-trials distance
% Z                  z-score from expected values from eeg channels

% EXPECTED VALUES  
%relative band power
%      DELTA   THETA     ALPHA   BETA1   BETA2   GAMMA1   GAMMA2
%b(1,:)=[0.32    0.197     0.15    0.11    0.08    0.069   0.054]; %mean
b(1,:)=[0.38    0.19      0.21    0.10    0.06    0.03     0.02]; %mean
b(2,:)=[0.08     0.05     0.05    0.03    0.02     0.02     0.02];  %std
%distributions
STD=[4.6   1.67];
RN=[25.2  9.2];
SK=[0   0.28];
KU=[2.9 0.52];
d=[STD' RN' SK' KU'];

if ~nargin
    help checksegments
else
  
[e,~,tr]=size(M);

if ~exist('verbose','var') || isempty(verbose)
verbose=0;
end
if ~exist('freq','var') || isempty(freq)
freq=0;
end
if ~exist('sr','var') || isempty(sr)
%examine distributions...    
if verbose
fprintf('Check distributions...') 
end
V=zeros(e,tr,4);
for n=1:tr
V(:,n,1)=std(M(:,:,n),1,2);
V(:,n,2)=range(M(:,:,n),2);
V(:,n,3)=skewness(M(:,:,n),1,2);
V(:,n,4)=kurtosis(M(:,:,n),1,2);
end
i=all(V(:,n,2)==0,2); %flat channel
if verbose
fprintf(' done.') 
fprintf('\n') 
fprintf('Absolute and relative distribution distances...') 
end
%distance from expected distributions...
Z=V;
for n=1:4
Z(:,:,n)=abs(Z(:,:,n)-d(1,n))./d(1,2);
end
zVtrans=zscore(V(~i,:,:),1,1); %electrodes distances (for each trial)
zVlong=zscore(V(~i,:,:),1,2);  %trials distances (for each electrode)
X.mztrans=nan(e,4);
X.mztrans(~i,:)=squeeze(mean(abs(zVtrans),2)); %average between-electrodes distance
X.mzlong=squeeze(mean(abs(zVlong),1)); %average between-trials distance
% 
% zVtrans=zscore(Z(~i,:,:),1,1); %electrodes distances (for each trial)
% zVlong=zscore(Z(~i,:,:),1,2);  %trials distances (for each electrode)
% Xz.mztrans=nan(e,4);
% Xz.mztrans(~i,:)=squeeze(mean(abs(zVtrans),2)); %average between-electrodes distance
% Xz.mzlong=squeeze(mean(abs(zVlong),1)); %average between-trials distance
if verbose
fprintf(' done.\n') 
end



else %power analysis   
if verbose    
fprintf('Check power...')     
end
[~,f,P2]=PSD_FFT(M,sr,2,2);
if ~freq
P2=f2b(P2,f,2);
%make it relative 
P2=P2./repmat(sum(P2,2),[1 size(P2,2) 1]);
V=zeros([size(P2,1) size(P2,3) size(P2,2) ]);
for n=1:size(P2,2)
    V(:,:,n)=P2(:,n,:);
end
if verbose
fprintf(' done.\n') 
fprintf('Absolute and relative band power distances...') 
end
%distance from expected distributions...
Z=V;
for n=1:7
Z(:,:,n)=abs(Z(:,:,n)-b(1,n))./b(1,2);
end
zVtrans=zscore(V,1,1); %electrodes distances (for each trial)
zVlong=zscore(V,1,2);  %trials distances (for each electrode)
X.mztrans=squeeze(mean(abs(zVtrans),2)); %average between-electrodes distance
X.mzlong=squeeze(mean(abs(zVlong),1)); %average between-trials distance
else
V=P2;
if verbose
fprintf(' done.\n')     
fprintf(['Absolute frequency-wise spectrum output [fq=linspace(' num2str(f(1)) ',' num2str(f(end)) ',' num2str(numel(f)) ')]...']) 
end
end
if verbose
fprintf(' done.\n') 
end
end
end






