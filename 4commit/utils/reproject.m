function    P=reproject(X,W,MASK,sr,epochs,chanlocs,i,method,forceblink)
 %    P=reproject(X,W,MASK,sr,epochs,chanlocs,i,method,forceblink);
 %Selective component removal:
 %method 1 blink handling + focal reprojection     [1]
 %method 2 blink handling
 %method 3 component removal

if ~nargin
    help reproject
    return
end


if ~exist('method','var')||isempty(method)
    method=1;
end
if ~exist('forceblink','var')||isempty(forceblink)
    forceblink=1;
end

veog=findchannel(chanlocs,'VEOG');
heog=findchannel(chanlocs,'HEOG');

btr=all(isnan(MASK));
bch=all(isnan(MASK),2);
if size(X,1)==numel(bch)-1
    X=ref(X);
end

[~,~,v_u2tr,v_u2TR]=epochs2vect(epochs, [],~btr); %get reshaping indices
ntp=median(epochs(:,3)-epochs(:,2))+1;

if method~=3
Xtemp=X; %for EYE
end

X=X(~bch,:,:);
if size(X,1)==size(W,2)+1 && all(X(end,1:10)==0) %ICA on data without reference in
    X(end,:,:)=[];
    bch(end)=1;
end
if numel(i)==1 && isnan(i)
   i=zeros(size(W,1),1);
end

if any(i) || forceblink
display('Full length components strengths')  
A=data2act(X,W);
if method<3 %partial removal 
T=ones(size(A)); %this temporal filter reprojects all
T=resh(T,ntp,v_u2tr);
if ntp/sr<5 %special case for blinks in short trials
display('Short trials: peak blinks detection strategy')
%I will localize blinks and remove components in corresponding trials
Xtemp=resh(Xtemp,ntp,v_u2TR);
BL=findocularsignals(Xtemp,sr,chanlocs,find(~bch),1,1);
Xtemp=Xtemp(~bch,:,:);
Atemp=data2act(Xtemp,W);
c=corr(BL(:,:)',Atemp(:,:)');
[bl,v]=findmax(abs(c));
clear *temp
if v>.3 %minimum to define it the blink component 
BL=act2proj(A,W,bl);
x=zscore(cell2mat({chanlocs.X}')); %"anteriorness"
y=zscore(cell2mat({chanlocs.Y}')); %"rightness"
z=zscore(cell2mat({chanlocs.Z}')); %"highness"

ant=x>1 & z<0.5*max(z);
refBL=z<-1;

ant=ant(~bch);
refBL=refBL(~bch);
BL=mean(BL(ant,:))- mean(refBL(ant,:));

%detecting blinks 
th=15;
BL=fqfilter(BL,[.5 40],sr)'; 
% [~,BL]=findpeaks(BL,'minpeakheight',th, 'minpeakdistance',round(sr/10));
BL=BL>th;
try
% BL=convertind(BL,size(A,2));
BL=extendtwin(BL,round(sr/4)); 
i(bl)=0;
BL=resh(BL,ntp,v_u2tr)==1;

%all trials where blinks are present
BL=repmat(any(BL,2),[1 size(BL,2) 1])==0; 
T(bl,:,:)=BL; %keep blink component in selected trials
catch err
    warning(err.message)
end
else
    display('**No BL component was detected**')
end 
end

A=resh(A,ntp,v_u2tr);  

if method==1
if any(i)  
    display('retain weak activations')
    tW=maketemporalW(A,W,i);    
    T(i,:,:)=tW;  
end
fprintf('Focal reprojection...')
else
end

A=A.*T;
P=pinv(W)*A(:,:);
P=reshape(P,[],size(A,2),size(A,3));
fprintf(' done.\n')
bch=all(isnan(MASK),2);


elseif method==3 %complete removal
fprintf('Full component removal...')
A=resh(A,ntp,v_u2tr);
if islogical(i);
P=act2proj(A,W,-(find(i)));    
else
P=act2proj(A,W,-i);    
end
fprintf(' done.\n')
end


else
display('No component to be removed')    
P=X;
end

