function [P,iW] = act2proj(A,W,i,S,rM)
%             [P,iW] = act2proj(A,W,i,S,rM)
%Dimensional check and basic algebra to reconstruct from weights and
%sphere(W,S) or from the inverse of the unmixing matrix (iW=inv(WS)) the contribution of the
%activations A of components 'i' in native space, i.e. their projection:
%             P = act2proj(A,W,i,S) 
%             P = act2proj(A,W,i,[]) [eye(ncomps))]
%             P = act2proj(A,iW,i)
%Basics:
%      A=WSX
%      iWA=IX,  where iW=inv(WS);
%      iWA=X, and for single projections...
%      Pi=iW(:,i)* A(i,:)
%Here:
%   A: (comp,tp,tr)/(comp,tp) -> P : (ch,tp,tr)/(ch,tp)
%   W: (comp,ch)
%   S: (ch,ch)   [eye(size(W,2))]
%   U=WS;
%   i can be  linear index(es) of component(s), if negative components are removed!
%BUT:
%     if i=0 maxmode is run (see "maxmode*")
%     if size(i)==size(A), i is a temporal weighting matrix:  P = iW* A.*i;
%The rows of P are the time courses of the chosen component at
%the respective scalp sensors, in the original input data units and polarity.
%Add rM to have original offsets.
%* maxmode: it projects all components at the maximum weights, note that in
%this case activations NOT channel activity is produced!
%
%             [P,iW] = act2proj(A,W,i,S,rM)

if nargin < 1
help act2proj
return
end

dims0=size(A);
[comp,tp,tr]=size(A);
[compw,ch] = size(W);
if tr>1
A=reshape(A, comp, tr*tp);
end
if tp*tr<comp %fewer tp then components means data are misoriented
    A=A';
   [ch,tp,tr]=size(A); 
end
if compw~=comp && compw==ch %W matrix should be transposed, note that for square W this cannot be determined!
    W=W';
    ch=compw;
end
if ~exist('i','var')||isempty(i)
    i=1:comp; %reproject all
end
if exist('S','var')
    if isempty(S)
        S=eye(size(W,2));
    end
    W=W*S; 
end

if size(W,1) == size(W,2)
  iW = inv(W);% inverse W 
else %overcomplete case
  iW = pinv(W);% pseudo-inverse W 
end

if size(i,1)==comp && size(i,2)==tp  %i is a temporal weighting matrix
    A=A.*i;
    i=1:comp;
end


if numel(i)==1 && i==0
       fprintf('Max projected activations...')
       P=A;
       for n=1:comp
        P(n,:) = max(abs(iW(:,n)))* A(n,:);
       end
       fprintf('done.\n')
else 
if any(i<0) %remove components
    i=abs(i);
    C=true(1,comp);
    C(i)=false;    
else %retain components
    C=false(1,comp);
    C(i)=true;
end
        %%%A=UX > X=iWA%%%%%
                   P = iW(:,C)* A(C,:);
        %%%%%%%%%%%%%%%%
end
        
if exist('rM','var')
   rM=rM(:);
   P=P+repmat(rM,[1 tp*tr]); 
end
 

if tr>1
    P=reshape(P, [size(P,1) dims0(2:end)]);
end   
   
   







