function [bTR,bCH,MASK]=cleaneeg(X,chanlocs,method,TH,sr)
%        [bTR,bCH,MASK]=cleaneeg(X,chanlocs,method,TH,sr)

if ~nargin
    help cleaneeg
    return
end

if ~exist('TH','var')||isempty(TH)
    TH=3;
end
if ~exist('sr','var')||isempty(sr)
    sr=250;
end
if ~exist('method','var')||isempty(method)
    method='FULL';
end



if method(1)=='1' || method(1)=='F'
%% CLEANUP1
display('**CLEANUP ROUTINE 1**')
[bTR1,bCH1,maxZ,mzZ,centZ]=findbad(X,chanlocs,'Z',TH,sr); %this initializes bTR1/bCH1
btr=findbad(X(~bCH1,:,~bTR1),chanlocs(~bCH1),'PCA',TH,sr);
bTR1(~bTR1)=btr;
fprintf('    Bad trials: %s\n' ,num2str(sum(bTR1)))   
fprintf('    Bad channels: %s\n' ,num2str(sum(bCH1))) 
[~,bch]=findbad(X(:,:,~bTR1),chanlocs,'PCA',TH);
bCH1=bCH1 & bch;
fprintf('    Bad channels after retaining correlated ones: %s\n' ,num2str(sum(bCH1))) 
%RESIDUAL HIGH V segments
%conservative threshold (1mV, blinks are up to 1mV)
HV=squeeze(any(abs(X)>1000,2));
if any(any(HV(~bCH1,~bTR1)))
   fprintf('Handling detected high V segments...\n')    
   %they are bad trials unless it is about some specific channels:
   %- high range peripherical channels (likely to be instable)
   [~,bch1]=findsparsechannels(any(HV(:,~bTR1),2),chanlocs);   
   %-"frequent" (>10%) high range channels
   bch2=mean(HV,2)>.1; 
   if mean(bch2)>.5 %too many 
        bch2=mean(HV,2)>.25; 
   end
   %-"unique" high range channels (involved in at least 1/2 trials with high V)
   bch3=mean(HV(:,any(HV(:,~bTR1))),2)>=.5;
   %but only if trials are numerous 
   if mean(any(HV(~bCH1,~bTR1)))<0.01
      bch3(:)=0;
   end
   bch=(bch1|bch2|bch3) & ~bCH1;   
   fprintf('    assigned to channels: %s\n' ,num2str(sum(bch)))    
   bCH1=bCH1 | bch; %update bad channels
   %residual, trial specific cases
   btr=any(HV(~bCH1,:)) & ~bTR1;
   fprintf('    assigned to trials: %s\n' ,num2str(sum(btr)))   
   bTR1=bTR1 | btr; %update bad trials
   fprintf(' done.\n')   
else
  display('No residual high V segments  detected')    
end
Z=cat(3,maxZ,mzZ,centZ);
MASK1=makemask(bCH1,bTR1,Z);
end


if method(1)=='F'
    X=X(~bCH1,:,~bTR1);
    chanlocs=chanlocs(~bCH1);
end


if method(1)=='2' || method(1)=='F'
%% CLEANUP2 
display('**CLEANUP ROUTINE 2**')
[bTR2,bCH2,maxZ,mzZ,centZ]=findbad(X,chanlocs,'FULL',TH,sr);
fprintf('    Bad channels: %s\n' ,num2str(sum(bCH2))) 
fprintf('    Bad trials: %s\n', num2str(sum(bTR2)))   
fprintf('Searching for local principal component activities...')
btr=findbad(X(~bCH2,:,~bTR2),chanlocs(~bCH2),'PCA',TH,sr);
bTR2(~bTR2)=btr;
fprintf('    Summing up to residual bad trials: %s\n', num2str(sum(bTR2)))   

Z=cat(3,maxZ,mzZ,centZ);
MASK2=makemask(bCH2,bTR2,Z);
end




if method(1)=='F'
    bCH=bCH1;
    bTR=bTR1;
    bCH(~bCH)=bCH2;
    bTR(~bTR)=bTR2;
    MASK=makemask(~bCH1,~bTR1,MASK2); 
elseif method(1)=='1'
    bTR=bTR1;
    bCH=bCH1;
    MASK=MASK1;
elseif method(1)=='2'
    bTR=bTR2;
    bCH=bCH2;
    MASK=MASK2;    
else
    error('cleaneeg method could not be determined!')
end



