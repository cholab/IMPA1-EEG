function [A,rM] = data2act(X,W,S)
%           A = data2act(X,W,S)
%           A = data2act(X,W,[])
%           A = data2act(X,U)
%Dimensional check and algebra to get components timecourses
%(activations) from weights and full projections (data):
%             A = WS(X - rM);
%   X: (ch,tp,tr)/(ch,tp) -> A : (comp,tp,tr)/(comp,tp)
%   W: (comp,ch)
%   S: (ch,ch)   [eye(size(W,2))]
%   U=WS;
%  rM: (ch,1) (row means to center data)
% Note: 
% A is not scaled to the input data units
% A polarities are not meaningful
% rM if not centered reconstruction of X is needed
%           A = data2act(X,W,S)

if nargin < 1
help data2act
return
end
  
  

[ch,tp,tr]=size(X);
[comp,chw] = size(W);

if tp<comp %fewer tp then components means data are misoriented
    X=X';
   [ch,tp,tr]=size(X); 
end
dims=size(X);

if tr>1
X=reshape(X, ch,[]);
end

nnans=~any(isnan(X));
X=X(:,nnans);
if chw~=ch && comp==ch %W matrix should be transposed, note that for square W this cannot be determined
    W=W';
    comp=chw;
end
if exist('S','var')
    if isempty(S)
        S=eye(size(W,2));
    end
    W=W*S;
end
rM=mean(X,2);%row means, it might matters for signal reconstruction
X=X-repmat(rM,[1 size(X,2)]); 

        %%%% A=WSX %%%%
            A=W*X;
        %%%%%%%%%%%%%%%
        
A2=nan([size(A,1) dims(2:end)]);
A2(:,nnans)=A;
% A=reshape(A,[size(A,1) dims(2:end)]);
A=A2;









