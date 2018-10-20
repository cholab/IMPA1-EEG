function [X,bTR,bCH,MASK]=finaldatacleanup(P,sr,chanlocs,bch,kTR,TH)
% [X,bTR,bCH,MASK]=finaldatacleanup(P,sr,chanlocs,bch,kTR,TH)

if ~nargin
    help finaldatacleanup
    return
end
display('final data cleanup')


if ~exist('TH','var')||isempty(TH)
    TH=3;
end
X=nan(length(bch),size(P,2),size(P,3));
if sum(~bch)==size(P,1)
X(~bch,:,:)=reref(P);  
else
X(~bch,:,:)=reref(ref(P));  
end



%% CLEANUP1
[tr,ch,MASK]=cleaneeg(X(~bch,:,:),chanlocs(~bch),'1',TH,sr);
if exist('kTR','var') && ~isempty(kTR)
    tr(kTR)=1;
end
ch1=bch; ch1(~ch1)=ch;
%% CLEANUP2
[tr1,ch1]=cleaneeg(X(~ch1,:,~tr),chanlocs(~ch1),'2',TH-.5,sr);
tr(~tr)=tr1; ch(~ch)=ch1; 
MASK=makemask(ch,tr,MASK(:,:,1),2); 
MASK=makemask(~bch,ones(size(MASK,2),1),MASK,1);


Zth=isnan(MASK)|MASK>3;
bTR=mean(Zth)>.5;
bCH=mean(Zth,2)>.5;


if sum(bCH)>sum(bch)
fprintf('Rereference after new bad channels were found... ' )    
X=nan(length(bch),size(P,2),size(P,3));
if sum(~bch)==size(P,1)
X(~bch,:,:)=reref(P);  
else
X(~bch,:,:)=reref(ref(P));  
end

X(~bCH,:,:)=reref(X(~bCH,:,:));
fprintf(' done. \n')
end

fprintf('Removing bad  %s trials... ' ,num2str(sum(bTR)))    
X=X(:,:,~bTR);
Zth=Zth(:,~bTR);
fprintf(' done. \n')



%handling bad segments when bad segments=bad channels
badsegments=any(Zth(~bCH,:));
[sel,pel]=findsparsechannels(bCH,chanlocs);
sel=sel&~pel;
if any(sel)
    fprintf('Interpolating  %s sparse bad channel(s)... ' ,num2str(sum(sel)))
    ibch=bCH*1;
    ibch(sel)=2;
    target=find(ibch(ibch~=1));
    ibch(sel)=0;  
    X(ibch==0,:,~badsegments)=repchan(X(ibch==0,:,~badsegments),target,chanlocs(ibch==0));
    fprintf(' done. \n')
end


Zth=isnan(MASK)|MASK>4;
Zth=Zth(:,~bTR);
badsegments=any(Zth(~bCH,:));
%handling residual bad segments
fprintf('Solving segment per segment...')
badsegments=find(badsegments);
for n=1:numel(badsegments)
ii=Zth(:,badsegments(n))==1;   
sbel=findsparsechannels(ii,chanlocs);
if any(sbel)
    ii=ii*1;
    ii(sbel)=2;
    target=find(ii(ii~=1));
    ii(sbel)=0;
    X(ii==0,:,badsegments(n))=repchan(X(ii==0,:,badsegments(n)),target,chanlocs(ii==0));
end
X(ii==1,:,badsegments(n))=NaN;
end
fprintf(' done. \n')

MASKIN=makemask(true(1,size(X,1)),~bTR,Zth,1);
MASKOUT=makemask(true(1,size(X,1)),~bTR,squeeze(any(isnan(X),2)),1);
MASKIN(isnan(MASKIN))=1;
MASKOUT(isnan(MASKOUT))=1;

fprintf('ALL DONE \n')









