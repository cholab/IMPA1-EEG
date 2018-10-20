function EYE=findocularsignals(X,sr,chanlocs,chans,advanced,justblinks)
%  EYE=findocularsignals(X,sr,chanlocs,chans,advanced,justblinks)
% Note that X and chanlocs shoould include all electrodes used for
% recording, use 'chans' to define a subset of good electrodes.
%   where  EYE=[BL VEOG HEOG REOG] unless 'justblinks'

if ~nargin
    help findocularsignals
    return
end
[el,tp,tr]=size(X);
if tr>1
    X=reshape(X,[el tp*tr]);
end
if el==size(chanlocs,2)-1 %missing ref?
    X=ref(X);
    el=size(X,1);    
end
if ~exist('advanced','var') || isempty(advanced) 
    advanced=0;
end
if ~exist('chans','var') || isempty(chans) 
    chans=ones(el,1);
else
    chans=convertind(chans,el);
end


veog=findchannel(chanlocs,'VEOG');
heog=findchannel(chanlocs,'HEOG');

if ~isempty(heog)
    eog=[heog veog];
    chans(eog)=0;
    heog=X(heog,:,:);    
else
    heog=0;
end
if ~isempty(veog)
    veog=X(veog,:,:);
else
    veog=0;
end





if ~exist('justblinks','var') || isempty(justblinks)
    justblinks=0;    
end



%Let's select some key electrodes
x=zscore(cell2mat({chanlocs.X}')); %"anteriorness"
y=zscore(cell2mat({chanlocs.Y}')); %"rightness"
z=zscore(cell2mat({chanlocs.Z}')); %"highness"
if exist('chans','var') && ~isempty(chans)
    x=x(chans);
    y=y(chans);
    z=z(chans);    
end
%%%%
% X=reref(X);
%%%%
if isa(X,'single')
    X=double(X);
end

% Define subsets of channels of interest
% ant -->BL and all EOGs
ant=x>1 & z<0.5*max(z);
if ~heog %more channels if EOG is not available
ant=ant | (x>.4 & z<-.8);
end
%here I am interested in the best possible estimate of specific signals,
%so... check for local (anterior) not correlated channels
D=zscore(nanmean(eldistance(X(ant,:))))<3;
ant(ant)=D; %and not include them in the search
% and refBL -->BL
refBL=z<-1;
%I will refine the others after removing blinks.



warning('off','stats:regress:RankDefDesignMat');



%BLINK
fprintf('Creating BL vector...')
%blinks are anterior events...
T=X(ant,:);
%whos variability is increased by referencing to most inferior electrodes
newref=mean(X(refBL,:));
T=T-repmat(newref,[size(T,1) 1]);
if heog     %of course it is in the EOG channels too (especially veog?)
   T=[T;veog;heog];
end

%it introduces activity <12Hz...
% T=fqfilter(T,20,sr,'low'); %maybe not, especially for segmented data

%its dipole is seen differently across electrodes, 
%going to PC captures its direction
[L,BL]=princomp(T'); %on cov matrix, absolute differences matters
%what component is it?
%Let's exploit that we know it should look like something like this
bl=getablink(sr);   
Tc=T;
for n=1:size(T,1)
    Tc(n,:)=conv(T(n,:),bl,'same');
end
[~,BLc]=princomp(Tc'); %on convolved data, blinks are consistently the first components
BLc=BLc(:,1);  %the first component
%I now use this to find what the corresponding component on nonconvolved data is.
%This step is not redundant as convolution alters the real blinks profiles.
[~,i]=max(abs(corr(BLc,BL)),[],2);
if advanced
  t=ceil(sr/100);  
  W=decompose(BL(1:t:end,unique([i 1:3]))');   
  W=W*pinv(L(:,1:3)); %complete transformation from data (to conserve absolute values)
  A=data2act(T,W);
  BL=BL(:,i);
  BL(abs(zscore(BL))<1.5)=0;
  [~,ic]=max(abs(corr(BL,A')),[],2); %marking corresponding component
  iW=pinv(W);
  [~,ie]=max(abs(iW(:,ic))); %marking max projection on scalp
  BL=iW(ie,ic)* A(ic,:);
else
  [~,ie]=max(abs(L(:,i))); %marking max projection on scalp
  BL=L(ie,i)*BL(:,i)';  
  %BL=fqfilter(BL,20,sr,'low'); %temporarily
end
% BL2=improveblinks(BL,sr);
fprintf(' done.\n')

    EYE=nan(4,numel(BL));
    EYE(1,:)=BL;

if ~justblinks    
try

       
%Let's remove the blinks from target anterior electrodes
%Since blinks vectors are ~spectrally clean this will conserve any higher frequency content
T=reref(X);
T=T(ant,:);
% if heog ;     T=[T;veog;heog]; end
T=regressout(T,BL);
%Any outlier here?
D=zscore(nanmean(eldistance(T)))<3;
T=T(D,:);
ant(ant)=D; 
%get a blink free EEG for other outliers
X2=regressout(X,BL);

%now I can define other subsets:
%% ant_r/l-->HEOG
ant_r=ant&y>.4;
ant_l=ant&y<-.4;
%% sup/inf--> VEOG
% depending on montage, electrodes inferior to the eyes can be more or less available...
sup=ant&z>-.8 & z<.9;    
if heog
  inf=z<-1.4;
else
  inf=ant&z<-.8;
end
if sum(inf)>1
D=zscore(nanmean(eldistance(X2(inf,:))))<3;
inf(inf)=D;
end

%% orbital,intermediate,Pz-->  REOG
if heog
 orbital=sup;    
else
 orbital=sup|inf;
end
orbital=orbital & x > 1.2;

intermidiate= x>=0 & x<1.1 & abs(y)<1;
D=zscore(nanmean(eldistance(X2(intermidiate,:))))<3;
intermidiate(intermidiate)=D;

Pz=z>0 & abs(y)<.6 & x<-.3 & x>-1.3;
D=zscore(nanmean(eldistance(X2(Pz,:))))<3;
Pz(Pz)=D;
clear X2



%HEOG
fprintf('Creating surrogate HEOG...')
%Corneo-retinal dipole activity is in the 4-20Hz range
% R=fqfilter(R,[4 25],sr,'pass');
T=fqfilter(T,25,sr,'low');
%marking left and right channels
ant_lr=ant_l(ant)*1+ant_r(ant)*2;
%initial cleanup by removing inconsistent components
[L,C]=princomp(T');
LL=mean(L(ant_lr==1,:));
LR=mean(L(ant_lr==2,:));
LM=mean(L(ant_lr==0,:));
%left-right gradient
ok=(LL>LM & LM>LR) | (LL<LM & LM<LR);
L(:,~ok)=0;
T=(C*L')';

L=T(ant_lr==1,:); R=T(ant_lr==2,:); 
mL=mean(L);       mR=mean(R);
%referencing to?
L=L-repmat(mR,[size(L,1) 1]);
R=R-repmat(mL,[size(R,1) 1]);
[LL,CL]=princomp(L');
[LR,CR]=princomp(R');
c=corr(CL, CR);
[v,i]=sort(c(:));
i=unique([i(1); i(v<-.5)]);
[l r]=ind2sub(size(c),i);
L=mean((CL(:,l)*LL(:,l)'),2)';    
R=mean((CR(:,r)*LR(:,r)'),2)';    
HEOG=L-R;
HEOG=fqfilter(HEOG,22,sr,'low'); %redundant?
fprintf(' done.\n')



%VEOG
fprintf('Creating surrogate VEOG...')
i=ant|inf;
T=reref(X);
T=X(i,:);
T=regressout(T,[BL;HEOG]);
%Corneo-retinal dipole activity is in the 4-20Hz range
% R=fqfilter(R,[4 25],sr,'pass');
T=fqfilter(T,25,sr,'low');
%marking superior and inferior channels
supinf=sup(i)*1+inf(i)*2;
ant_lr=ant_l(i)*1+ant_r(i)*2;
%initial cleanup by removing inconsistent components
[L,C]=princomp(T');
LSL=mean(L(supinf==1 & ant_lr==1,:));
LSR=mean(L(supinf==1 & ant_lr==2,:));
LI=mean(L(supinf==2,:));
%sup-inf gradient, not lateralized
ok=(LI>LSL & LI>LSR) | (LI<LSL & LI<LSR);
L(:,~ok)=0;
T=(C*L')';
SR=T(supinf==1 & ant_lr==1,:); 
SL=T(supinf==1 & ant_lr==2,:); 
SM=T(supinf==1 & ant_lr==0,:);
I=T(supinf==2,:);
mS=mean([SR;SL;SM]); mI=mean(I);
%referencing to?
SR=SR-repmat(mI,[size(SR,1) 1]);
SL=SL-repmat(mI,[size(SL,1) 1]);
SM=SM-repmat(mI,[size(SM,1) 1]);
I=I-repmat(mS,[size(I,1) 1]);

[LSR,CSR]=princomp(SR');
[LSL,CSL]=princomp(SL');
[LSM,CSM]=princomp(SM');
[LI,CI]=princomp(I');

%find sup vs inf anti correlated ones
try
c=corr(CSR, CI);
[v,i]=sort(c(:));
i=unique([i(1); i(v<-.5)]);
[r,ic1]=ind2sub(size(c),i);
catch
 r=0;
 ic1=0:size(CI,2);
end

try
c=corr(CSL, CI);
[v,i]=sort(c(:));
i=unique([i(1); i(v<-.5)]);
[l,ic2]=ind2sub(size(c),i);
catch
 l=0;
 ic2=0:size(CI,2);
end

try
c=corr(CSM, CI);
[v,i]=sort(c(:));
i=unique([i(1); i(v<-.5)]);
[m,ic3]=ind2sub(size(c),i);
catch
 m=0;
 ic3=0:size(CI,2);
end

try
check1=corr(CSR(:,r), CSL(:,l));
[r_l,l_r]=ind2sub(size(check1),find(check1>.2));
catch
 r_l=0:size(CSR,2);
 l_r=0:size(CSL,2);  
end
try
check2=corr(CSR(:,r), CSM(:,m));
[r_m,m_r]=ind2sub(size(check2),find(check2>.2));
catch
 r_m=0:size(CSR,2);
 m_r=0:size(CSM,2);
end
try
check3=corr(CSL(:,l), CSM(:,m));
[l_m,m_l]=ind2sub(size(check3),find(check3>.2));
catch
 l_m=0:size(CSL,2);
 m_l=0:size(CSM,2);
end

r=r_l(ismember(r_l,r_m)); %the one positively correlated with both other subgroups
l=l_r(ismember(l_r,l_m));
m=m_l(ismember(m_l,m_r));
ii=ic1(ismember(ic1,ic2)); ii=ii(ismember(ii,ic3)); 

try
SR=(CSR(:,r)*LSR(:,r)');
catch
SR=nan(size(CSR,1),1);
end
try
SL=(CSL(:,l)*LSL(:,l)');
catch
SL=nan(size(CSR,1),1);
end
try
SM=(CSM(:,m)*LSM(:,m)');
catch
SM=nan(size(CSR,1),1);    
end
try
I=(CI(:,ii)*LI(:,ii)');
catch
I=nan(size(CSR,1),1);    
end
S=nanmean([SR SL SM],2)'; I=nanmean(I,2)';
VEOG=S-I;
VEOG=fqfilter(VEOG,22,sr,'low'); %redundant?
fprintf(' done.\n')





%REOG
fprintf('Creating surrogate REOG...')
%We focus on periorbital and Pz electrodes 
i=orbital|intermidiate|Pz;
T=reref(X);
T=T(i,:);
T=regressout(T,[BL;HEOG;VEOG]);
%Saccadic spike potentials are above 20Hz...
T=fqfilter(T,30,sr,'high');
%marking orbital, intermidiate and Pz channels
oip=orbital(i)*1+intermidiate(i)*2+Pz(i)*3;
%initial cleanup by removing inconsistent components
%referencing to?
mT=mean(T(oip==2,:));
T=T-repmat(mT,[size(T,1) 1]);

[L,C]=princomp(T');
LO=mean(L(oip==1,:));
LP=mean(L(oip==3,:));
LI=mean(L(oip==2,:));
%orb-pz gradient
ok=(LO>LI & LI>LP) | (LO<LI & LI<LP);
L(:,~ok)=0;
T=(C*L')';
mO=mean(T(oip==1,:));
mP=mean(T(oip==3,:));

REOG=mP-mO;
REOG=fqfilter(REOG,30,sr,'high'); %redundant?
w=zscore(abs(wavtransform(REOG,[35:5:80],sr)),1,2);
[l,c]=princomp(w');
c=c(:,find(mean(l>0)>.8,1))';
[~,c]=findpeaks(c,'minpeakheight',1,'minpeakdistance',round(sr/3)); %up to 3 per second
c=repmat(c',[1 2])+repmat(round([-sr*0.03 sr*0.03]),[length(c) 1]); %it lasts up to 25ms
i=epochs2vect(c);
i(i<1|i>numel(REOG))=[];
REOG2=zeros(size(REOG));
REOG2(i)=REOG(i);
REOG=REOG2;
fprintf(' done.\n')

    EYE(2:end,:)=[VEOG(:) HEOG(:) REOG(:)]';

end

if tr>1
    EYE=reshape(EYE,[4 tp tr]);
end
else
   if tr>1
    EYE=reshape(EYE(1,:,:),[1 tp tr]);
   end 
end



return


%refine EYE vectors for ICA
%% BL
function BL2=improveblinks(BL,sr,tp)
%How do I suppress non-blink peak activity?
%A simple threshold based approach might work
%Note that now BL has real absolute values so I can set it to be 50microV
th=50;
% BL2=BL-smooth(BL,sr*2)'; % removing slow drifts
%detecting blink peaks
[~,pk]=findpeaks(BL,'minpeakheight',th,'minpeakdistance',round(sr/5));
%checking  sign so that peaks are in the blink-up direction
%NOTE: given the referencing this should never happen!
% [~,pk_down]=findpeaks(-BL2,'minpeakheight',th'minpeakdistance',round(sr/5));
% if numel(pk)<numel(pk_down); BL2=-BL2; pk=pk_down; end 
%I will remove high frequency content from everything outside blink peaks
%Define blink peak boundaries

%Exploiting average blink
maxblinkduration=round(sr*1);
% REMEMBER BLINKS AT TRIAL BOUNDARIES!!!!
[~,ipk]=findpeaks(BL,'minpeakheight',th,'minpeakdistance',maxblinkduration*2);
pkt=zeros(size(BL));
pkt(ipk)=1;
pkt=reshape(pkt',tp,[])';

pkt(:,1:maxblinkduration)=0;
pkt(:,end-maxblinkduration:end)=0;
pkt=pkt';pkt=pkt(:);
pkt(end-maxblinkduration*2,end)=0;
pkt=find(pkt);

    i=epochs2vect([pkt-maxblinkduration pkt+3*maxblinkduration]);
    BLtr=reshape(BL(i)',[],numel(pkt));
    
    i=ttest(gradient(BLtr(1:maxblinkduration,:)'));
    on=find(~i,1,'last');
    BLtr=BLtr-repmat(mean(BLtr(1:on,:)),[size(BLtr,1) 1]);  
    i=ttest(BLtr(round(maxblinkduration*1.5):end,:)');    
    off=round(maxblinkduration*1.5)+find(~i,1,'first');

%Exploiting derivatives 
% D1=gradient(BL2); 
% D2=gradient(D1);
% D3=gradient(D2);
% I=true(size(BL2)); %non peak points
% I2=true(size(BL2)); %safe non peak points
% for n=1:numel(pk)  
%forward
% [~,i1]=findpeaks(abs(D1(pk(n):end)),'npeaks',1);
% [~,i2]=findpeaks(abs(D2(pk(n)+i1:end)),'npeaks',1);
% [~,i3]=findpeaks(abs(D3(pk(n)+i1+i2:end)),'npeaks',1);
% i4=find(abs(BL2(pk(n)+i1+i2+i3:end))<20,1,'first');
% upb=pk(n)+i1+i2+i3+i4;
% %backward
% [~,i1]=findpeaks(fliplr(abs(D1(1:pk(n)-1))),'npeaks',1);
% [~,i2]=findpeaks(fliplr(abs(D2(1:pk(n)-i1-1))),'npeaks',1);
% [~,i3]=findpeaks(fliplr(abs(D3(1:pk(n)-i1-i2-1))),'npeaks',1);
% i4=find(abs(fliplr(abs(BL2(1:pk(n)-i1-i2-i3-1))))<20,1,'first');
% lowb=pk(n)-i1-i2-i3-4;
% if lowb<1; lowb=1; end
% if upb>numel(I2); upb=numel(I2); end
% I(lowb:upb)=false;
% adjlowb=lowb-round(sr/10); if adjlowb<1; adjlowb=1; end
% adjupb=upb+round(sr/10);   if adjupb>numel(I2); adjupb=numel(I2); end
% I2(adjlowb:adjupb)=false;
% end
% I=~(~I&BL2>0);
% I2=~(~I2&BL2>0);
% BLs=0.5*(smooth(BL2,round(sr/10))'+BL);
% BL2(I)=BLs(I);
% BLs=fqfilter(BLs,10,sr,'low')/10;
% BL2(BL2<-1)=-1;
% BL2(I2)=BLs(I2);
% BL2=smooth(BL2,round(8*(sr/250)))';




% %%VEOG-HEOG
% HEOG=smooth(H3,round(sr/4));
% art=abs(HEOG)>200;
% lateral=abs(HEOG)>10;
% HEOG(~lateral|art)=0;
% 













% % % [L,PC]=princomp(Rf(ant_lr>0,:)');
% % % ant_lr=ant_lr(ant_lr>0);
% % 
% % % a way to index the degree of asimmetry,  i.e. left-right dipole:
% % asimmetry=[mean(L(ant_lr==1,:))' mean(L(ant_lr==2,:))'];
% % a=sum(asimmetry>0,2)==1; %opposite sign components only
% % asimmetry=abs(diff(asimmetry,1,2))';
% % asimmetry(~a)=0;
% % %not driven by isolate electrodes (70% in the same direction)
% % a=abs(mean(L(ant_lr==1,:)>0)-0.5)>.2 &  abs(mean(L(ant_lr==2,:)>0)-0.5)>.2;
% % asimmetry(~a)=0;
% % asimmetry=asimmetry/max(asimmetry);
% % if heog
% %     h=abs(L(:,end))';
% %     h=h/max(h);
% %     asimmetry=asimmetry.*h;
% %     asimmetry=asimmetry/max(asimmetry);
% % end
% % % asimmetry(zscore(asimmetry)<0)=0;
% % %and we reconstruct the signal weighting components accordingly 
% % % PC=repmat(asimmetry,[size(PC,1) 1]).*PC;
% % W=repmat(asimmetry,[size(L,1) 1]);
% % R=(PC*(L.*W)')';    
% % %this should emphasize lef-right differences and reduce SNR in the
% % %calculation of the surrogate EOG:
% % Lt=R(ant_lr==1,:);
% % Rt=R(ant_lr==2,:);
% % if heog
% %        %I assign real eog channels to one of the scalp subset
% %         v=R(numel(ant_lr)+1,:);
% %         h=R(numel(ant_lr)+2,:);
% %         [~,i]=max(corr(v', [mean(Lt); mean(Rt)]'));
% %         if i==1
% %             Lt=[Lt; v];
% %         else
% %             Rt=[Lt; v];    
% %         end
% %         [~,i]=max(corr(h', [mean(Lt); mean(Rt)]'));
% %         h=h*1.5; %overweighting h for heog
% %          if i==1
% %             Lt=[Lt; h];
% %         else
% %             Rt=[Lt; h];    
% %          end
% % end
% % HEOG=[mean(Lt) mean(Rt)];
% % HEOG=fqfilter(HEOG,[4 25],sr,'pass');
% % HEOG=fqfilter(HEOG,22,sr,'low');
% % 
% % fprintf(' done.\n')
% % 
% % 




% % %VEOG
% % fprintf('Creating surrogate VEOG...')
% % % Let's dampen BL and  HEOG   contribution from the residual signal of anterior electrodes
% % i=ant|inf;
% % R=X(i,:);
% % if heog   
% %     R=[R;veog;heog];
% % end
% % for n=1:size(R,1)
% % [~,~,R(n,:)]=regress(R(n,:)',BL');
% % end
% % %Corneo-retinal dipole activity is still in the 4-20Hz range
% % R=fqfilter(R,[3 25],sr,'pass');
% % for n=1:size(R,1)
% % [~,~,R(n,:)]=regress(R(n,:)',HEOG');
% % end
% % 
% % 
% % %marking superior and inferior channels
% % ant_infsup=inf(i)*1+sup(i)*2;
% % [L,PC]=princomp(zscore(R'));
% % %now indexing the degree of asimmetry highlights the inf-sup dipole
% % asimmetry=abs(mean(L(ant_infsup==1,:))-mean(L(ant_infsup==2,:)));
% % asimmetry=asimmetry/max(asimmetry);
% % if heog
% %     h=abs(L(:,end-1))';
% %     h=h/max(h);
% %     asimmetry=asimmetry.*h;
% %     asimmetry=asimmetry/max(asimmetry);
% % end
% % asimmetry(zscore(asimmetry)<0)=0;
% % %and we reconstruct the signal weighting components accordingly 
% % PC=repmat(asimmetry,[size(PC,1) 1]).*PC;
% % R=(PC*L')';    
% % i=find(ant_infsup);
% % if heog
% %     i=[i'  numel(ant_infsup)+(1:2) ];
% % end
% % VEOG=mean(R(i,:));
% % fprintf(' done.\n')
% % 
% % 
% % 
% % 

% % %marking anterior right, anterior left and Pz electrodes
% % lrp=ant_l(i)*1+ant_r(i)*2+Pz(i)*3;
% % if heog
% %     lrp=[lrp;0;0];
% % end
% % [L,PC]=princomp(zscore(Rf'));
% % % this provides a way to index the degree of asimmetry, i.e. ant-Pz dipole
% % asimmetry=abs(mean(L(lrp~=3,:))-mean(L(lrp==3,:)));
% % asimmetry=asimmetry/max(asimmetry);
% % %but also:
% % %-the left-right assimetry should be penalized...
% % lrasimmetry=abs(mean(L(lrp==1,:))-mean(L(lrp==2,:)));
% % lrasimmetry=lrasimmetry/max(lrasimmetry);
% % asimmetry=asimmetry.*(1-.5*lrasimmetry);
% % %-and the consistency in Pz rewarded
% % Pzconst=L(lrp==3,:);
% % Pzconst(abs(Pzconst)<=.1)=NaN;
% % Pzconst(Pzconst<-.1)=-1;
% % Pzconst(Pzconst>.1)=1;
% % Pzconst=abs(nanmean(Pzconst));
% % Pzconst(isnan(Pzconst))=1;
% % asimmetry=asimmetry.*Pzconst;
% % asimmetry=asimmetry/max(asimmetry);
% % asimmetry(zscore(asimmetry)<-1)=0;
% % %and we reconstruct the signal weighting components accordingly 
% % PC=repmat(asimmetry,[size(PC,1) 1]).*PC;
% % R=(PC*L');    
% % REOG=mean(R(:,lrp~=3 ),2)-mean(R(:,lrp==3),2);
