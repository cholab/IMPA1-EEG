function R=cfilter(S,W,c,m)
%              R=cfilter(S,W,c,m)
% It accepts scores (S) and a weight matrix (W), such as Loadings from PCA
% or Weights from ICA, and  reconstruct a signal (R) by retaining components
% 'c' (or removing 'c' if negative):
%                   R=cfilter(S,W,c,m)
%  m:   mean from original signal (it is added back if provided)
%Note: the number of the original sources/signals (or its rank) is expected to
% be equal to d2 in d1Xd2 S(i.e. time is vertical) and W matrices
% Optional call for pca/ica (mode):
%                             R=cfilter(X,mode,c)
%                             R=cfilter(X)  [pca] [1]     (robust averager)

if ~nargin
    help cfilter
else
    flip=false;
    if nargin==1
        W='pca';
        c=1;
    end
    if ischar(W)
    if size(S,1)<size(S,2)
       flip=true;
       S=S';
    end
    if strcmpi(W(1),'p')
        m=mean(S);  
        %princomp expects time to be vertical
        [W,S]=princomp(S);
    elseif strcmpi(W(1),'i')
        %fastica expects time to be horizontal      
%     [S,separatingM,mixingM]=fastica(S');        
        [S,W,~]=fastica(S');
        %runica expects time to be horizontal
%         [W,SPH,bpcompvars,bias,signs,lrates,ACT]=runica(S);
%         unmixingM = W*SPH {if sphering off -> eye(chans)}
%        winv = inv(W*SPH);
        S=S';
    end
    end
dS=size(S);
dW=size(W);
dimok=true;
if dS(2)~=dW(2)
   if dS(1)==dW(2) 
       S=S';
       dS=size(S);
   else
       dimok=false;
   end
end

if ~dimok
    display('S/W dimensions missmatch: check input size and orientation!' )
    display(' Note: the number of the original sources/signal is taken to be equal to d2 in d1Xd2 S and W matrices')
    R=[];
else
if  exist('m','var') && ~isempty(m) && length(m)~=dS(2) 
    display('mean-W dimensions missmatch: check provided mean!' )
    R=[];
else    
    
if any(c<0) %remove components
    c=abs(c);
    C=true(1,dS(2));
    C(c)=false;
else %retain components
    C=false(1,dS(2));
    C(c)=true;
end

% Inner matrix dimensions must agree:
% reconstructed = repmat(mean(x,1),n,1) + score(:,1:ndim)*coeff(:,1:ndim)';
% residuals = x - reconstructed;
                R=S(:,C)*W(:,C)';    


if exist('m','var') && ~isempty(m)
    R=R+repmat(m(:)',[dS(1) 1]);
end
end
if flip
    R=R';
end
end
end


