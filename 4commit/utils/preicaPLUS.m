function [MASK,mask]=preicaPLUS(X,sr,chanlocs,chans,TH,method,suppressref)
%Initial channels/trials quality assessment and removal of deviants that
%are unlikely to be recovered by/would affect decomposition:
%        MASK=preicaPLUS(X,sr,chanlocs,chans,TH,method,suppressref)
%'X' data (electrodes x time points x trials)
%'sr' sampling rate [250]
%'chans' subset of channels to process including reference [all]
%'chanlocs' montage structure
%'reference' index of reference electrode (it will be added if exceeds n electrodes)
%'TH' sd threshold for all distribution based selections [3]
%'method'it refers to levels of preprocessing info included in MASK 
%routines    1   2           1  2
%        1   1--(2)--ocular--3--4--5 mask=[1 2 3 4 5]   [1]
%        2   1-------ocular--3--4--5 mask=[1   3 4 5]
%        3   1---2                   mask=[1 2      ]
%        4   1                       mask=[1        ]
%'suppressref' it suppresses reference check/adjustment [0]
%       MASK=preicaPLUS(X,sr,chanlocs,chans,TH,method,suppressref)


if ~nargin
    help preicaPLUS
    return
end
fprintf('\n_______INITIAL DATA CLEANUP <preica.m_______\n')



if ~exist('TH','var')||isempty(TH)
    TH=3;
end
if ~exist('suppressref','var')||isempty(suppressref)
    suppressref=0;
end
if ~exist('method','var')||isempty(method)
    method=1;
end
switch method
    case 1
        ocularout=1;
        fullstep1=1;
        mask=1:5;
    case 2
        ocularout=1;
        fullstep1=0;
        mask=[1 3:5];        
    case 3
        ocularout=0;
        fullstep1=1;
        mask=1:2;        
    case 4
        ocularout=0;
        fullstep1=0;
        mask=1;           
end
if ~suppressref
[X,other]=ref(X);
else
other=0;    
end
[nel,ntp,ntr]=size(X);
fprintf('Data size: %s channels X %s time points X %s trials\n' ,num2str(nel),num2str(ntp),num2str(ntr))    
if exist('chans','var') && ~isempty(chans)
   X=X(chans,:,:);
   chanlocs=chanlocs(chans);
end


display('____WORKING ON RAW DATA________________________')
%% CLEANUP1 on RAW 
[bTRr,bCHr,MASKr]=cleaneeg(X,chanlocs,'1',-TH,sr);
if other
    bCHr(other)=1;
end
if fullstep1
%% CLEANUP2 on RAW 
%I add it for comparison with the final post eye removal cleanup, but I 
%will proceed with the bCH/bTR identified so far
[btr,bch]=cleaneeg(X(~bCHr,:,~bTRr),chanlocs(~bCHr),'2',TH,sr);
bTRr2=bTRr;  bTRr2(~bTRr2)=btr;
bCHr2=bCHr;  bCHr2(~bCHr2)=bch;
MASKr2=makemask(bCHr2,bTRr2,MASKr); 
end

if ocularout
%% GET bad ch/tr and EYE removed X 
%Now data should be clean enough to reliably extract ocular info
EYE=findocularsignals(X(:,:,~bTRr),sr,chanlocs,~bCHr);
%and remove it in order to improve detection of other sources of noise
fprintf('%s','Regressing out EYE variability... ') 
X=regressout(X(~bCHr,:,~bTRr),EYE);
fprintf(' done.\n')    

display('____WORKING ON OCULAR NOISE REMOVED DATA__________')
%% CLEANUP1 on EYEREMOVED
chanlocs=chanlocs(~bCHr);
[bTRe,bCHe,MASKe]=cleaneeg(X,chanlocs,'1',TH,sr);
%% CLEANUP2 on EYEREMOVED
[btr,bch]=cleaneeg(X(~bCHe,:,~bTRe),chanlocs(~bCHe),'2',TH,sr);
bTRe2=bTRe;  bTRe2(~bTRe2)=btr;
bCHe2=bCHe;  bCHe2(~bCHe2)=bch;
MASKe2=makemask(bCHe2,bTRe2,MASKe); 



display('____SUMMARY___________________________________')
%% final bad channels/trials count
btr=bTRr; bch=bCHr;
btr(~btr)=bTRe2;
bch(~bch)=bCHe2;
fprintf('    Total bad trials: %s\n' ,num2str(sum(btr)))   
fprintf('    Total bad channels: %s\n' ,num2str(sum(bch))) 

%% MASKS
%fill the rest with maxZ info
MASKe=makemask(~bCHr,~bTRr,MASKe); 
MASKe2=makemask(~bCHr,~bTRr,MASKe2); 
if fullstep1
MASK=cat(3,MASKr(:,:,1),MASKr2(:,:,1),MASKe(:,:,1),MASKe2(:,:,1));
else
MASK=cat(3,MASKr(:,:,1),MASKe(:,:,1),MASKe2(:,:,1));    
end
MASKre=MASKr(:,:,1);
MASKre(any(isnan(MASK),3))=NaN;
MASK=cat(3,MASK,MASKre);
else
if fullstep1   
MASK=cat(3,MASKr(:,:,1),MASKr2(:,:,1));   
else
MASK=MASKr;    
end
end
display('  ** MASKS CREATED**  ')
display('______________________________________________ ')






% % el=1:size(Z,1);
% % tr=1:size(Z,2);
% % el=repmat(el',[1 size(Z,2) ]);
% % tr=repmat(tr,[size(Z,1) 1]);
% % 
% % x=[el(:) tr(:)];
% % y=Z(:);
% % [b,~,s]=glmfit(x,y);
% % yfit = glmval(b, x,'identity');
% % 
% % y=Z(:)>5;
% % [~,~,s]=glmfit(x,y,'binomial');
% % b = glmfit([el(:) tr(:)],Z(:)>4,'binomial','link','probit');
% % yfit = glmval(b, x,'probit');
% % plot(x, y,'o',x,yfit,'-','LineWidth',2)




