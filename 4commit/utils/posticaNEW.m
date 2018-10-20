function [i,salvage]=posticaNEW(X,W,MASK,sr,epochs,chanlocs,EYE,ieye,keepall)
% [i,salvage]=postica(X,W,MASK,sr,epochs,chanlocs,EYE,ieye)

if ~nargin
    help postica
    return
end


btr=all(isnan(MASK));
bch=all(isnan(MASK),2);
if size(X,1)==numel(bch)-1
    X=ref(X);
end
[~,~,~,u2TR]=epochs2vect(epochs, [],~btr); %get reshaping indices
ntp=median(epochs(:,3)-epochs(:,2))+1;
X=resh(X,ntp,u2TR);

if ~exist('EYE','var')||isempty(EYE) || size(EYE,2)*size(EYE,3)~= size(X,2)*size(X,3)
EYE=findocularsignals(X,sr,chanlocs,find(~bch));
end
if size(EYE,3)==1 
    EYE=resh(EYE,ntp);
end
X=X(~bch,:,:);
if ~exist('ieye','var') ||isempty(ieye)
    ieye=0;
end

fprintf('%s','Unmixing activations... ')
if ieye
    if size(X,1)+numel(ieye)==size(W,2)+1 && all(X(end,1:10)==0) %ICA on data without reference in
     X(end,:,:)=[];
    end
    A=data2act(cat(1,X,EYE(ieye,:,:)),W);
else
    if size(X,1)==size(W,2)+1 && all(X(end,1:10)==0) %ICA on data without reference in
     X(end,:,:)=[];
    end
    %keyboard
    A=data2act(X,W);    
end
[comp,tp,tr]=size(A);
fprintf(' done.\n')


fprintf('%s','Projecting max activations... ')
[P,iW] = act2proj(A,W,0);

 
display('ICA activations') 
bA=binactivity(A);
display('ICA activations, trialwise') 
[tA,~]=explainedvar(A,iW,3);
tA=zscore(tA,1,2);

display('ICA activations, kurtosis') 
[kA,kTR]=findkurtotic(tA);

typA=zeros(comp,tr);%typical trials
for n=1:comp
    j=false(comp,1); j(n)=true;
    typical=tA(j,:)-nanmean(tA(~j,:));
    typical(:,kTR)=0;
    [~,typA(n,:)]=sort(typical,'descend');
end




fprintf('\n%s','C<->C distances...')
% CHECK DECOMPOSITION this kind of works,  not in a sensible way
[~,mA]=sort(tA,2,'descend');
mA=mA(:,1:2);
umA=unique(mA(:));
umA(ismember(umA,find(kTR)))=[];
% umA=unique(typA(:,1));
MI=hb_mi(reshape(A(:,:,umA),size(A,1),[])');
MI=transf(MI,'Z',2);
MI=MI>5;
[r1,c1]=find(MI);
MI(sub2ind(size(MI),c1,r1))=1;

CCb=zeros(size(bA,1));
CCb=repmat(CCb,[1 1 numel(umA)]);
for t=1:numel(umA)
   CCb(:,:,t)=corr(bA(:,:,umA(t))'); %corr(bA) in each trial
end
CCb=abs(CCb);
CCb=median(CCb,3)>.6;
CCb(eye(size(CCb))==1)=0;
% [r2,c2]=find(CCb);
TWINS=MI & CCb;
fprintf(':done.\n')





fprintf('%s','Finding high range activations... ')
Vr=checksegments(P);
Vr=Vr(:,:,2); %projected ranges
v=zeros(comp,10);
for n=1:comp
    v(n,:)=Vr(n,typA(n,1:10));
end
Awr=all(v>100,2);
fprintf(' done.\n')

fprintf('%s','Finding spectrally deviant activations... ')
Vsp=checksegments(P,sr);
Vsp=meandim(Vsp,2); %mean relative power
b(1,:)=[0.38    0.19      0.21    0.10    0.06    0.03     0.02]; %mean
b(2,:)=[0.08     0.05     0.05    0.03    0.02     0.02     0.02];  %std
Z=(Vsp-repmat(b(1,:),[size(Vsp,1) 1]))./repmat(b(2,:),[size(Vsp,1) 1]);
Asp=any(Z>5,2) & meandim(Vr,2)>50;
fprintf(' done.\n')

 
try
    
fprintf('\n%s\n','Z trialwise distances from target eeg parameters... ') 
X2=regressout(X,EYE);  %Remove blinks and ocular stuff...
[~,ZZ]=checkeeg(ref(X2,0),sr); %0:centered on dataset
ZZ=reshape(ZZ(:,:)',size(X,3),[])'; %this is [el1[1:11] el2[1:11]] ... x tr
btr=find(any(ZZ>2));
btr(ismember(btr,kTR))=[];
ZZ(isnan(ZZ))=0;
[~,Zpc,E]=princomp(zscore(ZZ(:,btr))');
if any(E/sum(E)>.1)
 Zpc=Zpc(:,E/sum(E)>.1)';
else
 Zpc=nanmean(Zpc(:,1),2)';
end
fprintf(' done.\n')

fprintf('%s','Z<->C')
C_Z=abs(corr([Zpc' tA(:,btr)']));
C_Z(eye(size(C_Z,1))==1)=0;
C_Z=C_Z(1:size(Zpc,1),size(Zpc,1)+1:end);
C_Z=squeeze(max(C_Z,[],1));
if ~isempty(C_Z)
    Az=C_Z(:)>.8;
else
    Az=false(comp,1);
end

catch
        Az=false(comp,1);
end
fprintf(':done.\n')






fprintf('\n%s','C<->EYE ')
if ~ieye
C_EYE=zeros(comp,4);
C_EYE=repmat(C_EYE,[1 1 tr]);
EYE=reshape(EYE,4,[],tr);
tr2=round(linspace(1,tr,10));
for t=tr2    fprintf('%s', '.'); end
for t=1:tr
       if any(t==tr2); fprintf('\b\b%s', ' '); end  
   C_EYE(:,:,t)=corr(A(:,:,t)',EYE(:,:,t)');   %corr(A,EYE) in each trial
end
C_EYE=abs(C_EYE); 
fprintf(':done.\n')
Aeye=false(comp,1);
[v,aeye]=max(abs(corr((A(:,:))',(EYE(:,:))')));

Aeye(aeye(v>.5))=true;
if v(3)>.2
   Aeye(aeye(3))=true; 
end
Aeye= Aeye | any(mean(C_EYE>.5,3)>.4,2);
else    
    [~,Aeye2]=max(abs(iW(end-(numel(ieye)-1):end,:)),[],2);
    Aeye=false(size(Awr));
    Aeye(Aeye2)=1;
end



fprintf('\n%s'  ,['WR components: ' num2str(find(Awr)')])
fprintf('\n%s'  ,['SP components: ' num2str(find(Asp)')])
fprintf('\n%s'  ,['Z components: ' num2str(find(Az)')])
fprintf('\n%s\n',['EYE components: ' num2str(find(Aeye)')])
i= Aeye(:) | Awr(:) | Az(:) | Asp(:);


if any(i)
twin=any(TWINS(i,:));
twin=twin(:);
fprintf('\n%s\n',['Found twins: ' num2str(find(twin)')])
if any(twin)
i=i|twin;
end

fprintf('\n%s\n',['Found kurtotic: ' num2str(find(i' & kA))])    
if any(i(kA))
i=i & ~kA(:);
end
end


salvage=false;
if ~exist('keepall','var') ||isempty(keepall)
if mean(i)>0.08 %like 10 for 129 components
  if sum(Aeye)<4;
  display('*** Too many components were found, resetting to Aeye ****')
  i=Aeye;
  else
  display('*** Too many components were found, use salvage mode? ****')
  salvage=true;
  end
end
end




if ~salvage
fprintf('\n%s','Iterative trialwise spectral quality pursuit ')
newC=1;
while any(newC)
fprintf('\n%s','Loop Start')  
if any(i)
 P2= act2proj(A,W,-(find(i)));
else
 P2= act2proj(A,W);    
end
if ieye
P2=P2(1:end-numel(ieye),:,:);
end
fprintf('\n%s','Estimating spectral content...')
% [~,ZZ]=checkeeg(X2,sr,0); %0:centered on dataset
[~,~,Z]=checksegments(reref(P2),sr);
Z=reshape(Z,size(Z,1),[]);
Z=reshape(Z',size(X,3),[])';
btr=any(Z>3);
if mean(btr)>.05
btr= (btr | all(Z<1)) & ~kTR; %add some variability, but not kurtotic trials
fprintf('\n%s','Finding main directions of deviant power...')
[~,pc,E]=princomp(zscore(Z(:,btr))');
Zpc=pc(:,E/sum(E)>.1)';
if isempty(Zpc)
   Zpc=mean(Z(:,btr));
end
C_Z=abs(corr([Zpc' tA(:,btr)']));
C_Z(eye(size(C_Z,1))==1)=0;
C_Z=C_Z(1:size(Zpc,1),size(Zpc,1)+1:end);
C_Z=squeeze(max(C_Z,[],1));


C_Z(i)=0;
Az=zeros(comp,1);
[v,az]=max(C_Z);
Az(az)=v;
Az=Az>.5;
Az(kA)=false;
fprintf('\n%s',['Related new component: ' num2str(find(Az)')])
newC=Az;
i=i|Az;
if sum(i)>9
    newC=0;
end
else
fprintf('\n%s\n','No bad segments found')        
newC=0;
end
end
end



if sum(i & ~Aeye)>10
  display('*** Too many components were found, use salvage mode? ****')
  salvage=true;
end




