function [kA,kTR]=findkurtotic(tA)
% [kA,kTR]=findkurtotic(tA)

if ~nargin
    help findkurtotic
    return
end

[mv,mtr]=max(tA');
[f,u]=occ(mtr);
f=f/sum(f);

%kurtotic trials
kTR1=find(mean(tA>2)>.5);
kTR2=u(f>0.25); 
kTR3=find(mean(tA)>2);
ktr=unique([kTR1(:)' kTR2(:)' kTR3(:)']);
kTR=false(1,size(tA,2));
kTR(ktr)=true;

kA1=kurtosis(tA')>80;
kA2=ismember(mtr,kTR1);
kA3=mv>8;
kA= (kA1 | kA2) & kA3; %kurtotic activations