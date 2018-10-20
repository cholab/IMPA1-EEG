function tW=maketemporalW(A,W,i)
%Localize weak component activity, both sub-trial (kurtosis) and trialwise
%(ratio of explained variability)
%      tW=maketemporalW(A,W,i)

if ~nargin
    help maketemporalW
    return
end
%sub-trial activity indexed by kurtosis
K=squeeze(kurtosis(A,[],2));
K=K(i,:);

%trialwise activity indexed explained variability
[~,R]=explainedvar(A,pinv(W),3);
R=R(i,:);

%keep component when strength is weak
tW= K<5 | R<.5 | zscore(R,1,2)<2;
tW=reshape(tW,[size(tW,1) 1 size(tW,2)]);
tW=repmat(tW,[1 size(A,2) 1]); 



