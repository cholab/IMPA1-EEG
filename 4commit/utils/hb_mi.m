function [MI,JP]=hb_mi(TS,nbin)
%
%        [MI,JP]=hb_mi(TS,nbin,lag)
%MI   Histogram based mutual information.
%JP   Jaccard "proximity' (1- Jaccard Distance)

if ~nargin
     help hb_mi
else
    
if ~exist('nbin','var')||isempty(nbin)
    nbin=10;
end

misigma=false;
if nargout>1
    misigma=true;
end

if ~iscell(TS)
    [r c]=size(TS);
    c=ones(1,c);
    TS=mat2cell(TS, r, c);
end

nts=size(TS,2); 
%avoid unnecessary comparisons
pairs=reshaper(nts);

MI = zeros(size(pairs,1),1); 
JP  = zeros(size(pairs,1),1); 
% if misigma
% MI_sigma = zeros(size(pairs,1));
% end



for np=1:size(pairs,1)    
     P2=h2loc(TS{pairs(np,1)},TS{pairs(np,2)},nbin);    % 2D distribution          
     P1x=sum(P2,2);    %marginal distribution for x
     P1y=sum(P2);       %marginal distribution for y               
         % compute the mutual information
      I1=[-sum((P1x(P1x~=0)).*log(P1x(P1x~=0))), -sum((P1y(P1y~=0)).*log(P1y(P1y~=0)))];  % entropies of Px and Py
      I2=-sum(P2(P2~=0).*log(P2(P2~=0)));  % entropy of joint distribution Pxy
      I2_syserr=(length(P1x(P1x~=0))+length(P1y(P1y~=0))-length(P2(P2~=0).*log(P2(P2~=0)))-1)/(2*length(TS{pairs(np,1)})); % standard error for estimation of entropy
      MI(np)=I1(1)+I1(2)-I2+I2_syserr; %MI=Hx+Hy-Hxy     
      JP(np)=MI(np)/I2; %1-Jaccard distance
      
%      if misigma     % compute the standard errors
%      P2xy=P1x*P1y;
%      i=(P2xy~=0 & P2~=0);
%      MI_sigma(np)=sqrt( sum( (log(P2xy(i)./P2(i))+MI(np)).^2 .*(P2(i).*(1-P2(i)))) /length(TS{pairs(np,1)}));
%      end    
end

 if isnan(MI(1,1))
     warning('Mutual information is NaN. Use smaller bin size.')
 end
     if np>1
     MI=reshaper([pairs MI]);
     JP=reshaper([pairs JP]);
     JP(eye(size(JP))==1)=1;
     end
end

% NORMALIZED MI?
%       v = sqrt((MI/Hl)*(MI/Hr)) ; 

% MI can also be understood as the expectation of the KL divergence of the
% univariate distribution p(x)  from the conditional distribution p(x|y):
% the more different the distributions p(x|y) and p(x), the greater the
% information gain.





function P=h2loc(x,y,nbin,~) %lag=0 version
% x=x(:); y=y(:);
% normalise the value range to [0 1)
x1 = (x - min(x)) / max(x - min(x)) - eps;
y1 = (y - min(y)) / max(y - min(y)) - eps;

% if length(x)<=lag, lag=length(x)-1; warning(['Lag is too large. Using ',num2str(lag),' instead.']), end

% this is the main trick: put the first data to values [1:1:nbin] and the second data to [nbin:nbin:nbin^2]
temp = fix(x1 * nbin) + (fix(y1 * nbin)) * nbin;
% call histc function (faster than hist)
P=fliplr(reshape(histc(temp,0:(nbin^2-1))',nbin,nbin)/length(x));
% function P=h2loc(x,y,nbin,lag)
% x=x(:); y=y(:);
% % normalise the value range to [0 1)
% x1 = (x - min(x)) / max(x - min(x)) - eps;
% y1 = (y - min(y)) / max(y - min(y)) - eps;
% 
% if length(x)<=lag, lag=length(x)-1; warning(['Lag is too large. Using ',num2str(lag),' instead.']), end
% 
% % this is the main trick: put the first data to values [1:1:nbin] and the second data to [nbin:nbin:nbin^2]
% temp = fix(x1(1:(end-lag)) * nbin) + (fix(y1((lag+1):end) * nbin)) * nbin;
% % call histc function (faster than hist)
% P=fliplr(reshape(histc(temp,0:(nbin^2-1))',nbin,nbin)/length(x));




