function viscomp(U,X,chanlocs,icomp,AREA,scalpkind,nrows,realvalue,time,events,labels)
%    viscomp(U,X,chanlocs,icomp,AREA,scalpkind,nrows,realvalue,time,events,labels)
%  -'U' (signals,channels) is the unmixing matrix (i.e. inv(L),W), multiple
%  U can be nested within a single cell input and each decomposition will 
%  be plot as a different column
%  -'X' (channel,time points) is data or for unmatching first dimension
%  (signals,time points) it is intepreted as signals (i.e. PC,A)
%  -'chanlocs' is channel location
%  -'icomp' are the subset of components to display
%  -'AREA' is of display of regional topographies
%  -'scalpkind' control scalp map kind or suppress it for scalpkind=0
%  -'nrows' forces rows number: if multiple U are provided only up to nrows
%  components are plotted, if a single U is provided the single
%  decomposition is arranged across ceil(signals/nrows) columns.
%  -'realvalue' max projections are plotted instead of signals and plots
%  are rescaled within the [-realvalue realvalue] range
%Notes:
%-U is used to define the identity of X, therefore only the umixing
%matrix is accepted as an input, X can be Signals but it is recognized as
%such only for not square U so that the number of rows matches
%unambiguously the number of signals and not the number of channels

if ~nargin
    help viscomp
    return
end

onlytimecourse=0;
if ~exist('scalpkind','var')
    scalpkind=[];
elseif ~scalpkind
    onlytimecourse=1;
end
if ~exist('realvalue','var')||isempty(realvalue)
    realvalue=0;
end



if iscell(U) %decompositions comparison mode
    nw=length(U);
else
    nw=1;
    U={U};
end   
if ~exist('nrows','var') || isempty(nrows)
    nrows=0;
    if nw>1
    for n=1:nw    
    nrows=max([nrows size(U{n},1)]);
    end
    end
end

if nw==1 && nrows %split single decomposition mode
if nrows<1
      nrows=ceil(sqrt(size(U{1},1))); 
%       adj=1;
%     while nrows~=ceil(nrows)
%       nrows=sqrt(size(U{1},1)-adj); adj=adj+1;
%     end
end
nw=ceil(size(U{1},1)/nrows);  
icomp=[]; %override whatever other input
end


   
hsplit=linspace(0,1,nw+1);
unmix=1;

for nnww=1:nw
if unmix    
[nsources,nchannels]=size(U{nnww}); 
M=pinv(U{nnww});  %mixing matrix

dim=size(X);
if ismember(dim(2),size(U{nnww})) % 2D transposed X
    display('Transposing TimeSeries input')    
    X=X';
    dim=size(X);
end
if dim(1)==nchannels %input is raw data
    A=data2act(X,U{nnww});
elseif dim(1)==nsources %input is activation
    A=X;
else
   error('Unmixing matrix-TimeSeries input missmatch!') 
end
if realvalue
   A=act2proj(A,U{nnww},0); %max projection
end
if exist('chanlocs','var') && ~isempty(chanlocs)
    scalpmode=1;
    if exist('AREA','var') && ~isempty(AREA)
    areamode=1;  
    if size(AREA,2)~=nchannels
    error('Unmixing matrix-AREA locations missmatch')    
    end
    else
    areamode=0;          
    if size(chanlocs,2)~=nchannels
    error('Unmixing matrix-Channel locations missmatch')    
    end
    end
else
     scalpmode=0;
end


if exist('icomp','var') && ~isempty(icomp)
M=M(:,icomp);
A=A(icomp,:,:);
end
if length(U)==1 && nw>1 %split single decomposition across different columns
unmix=0;
AA=cell(nw,1);
MM=cell(nw,1);
ii=unique([0:nrows:nsources nsources]);
for nii=1:nw
    AA{nii}=A(ii(nii)+1:ii(nii+1),:,:,:,:);
    MM{nii}=M(:,ii(nii)+1:ii(nii+1));    
end
end
else
    A=AA{nnww};
    M=MM{nnww};
end

[ncomp,ntp,~]=size(A);
if ~exist('time','var') || isempty(time)
    time=1:ntp;
end
if ~exist('events','var') 
    events=[];
end


if nrows
    vsplit=linspace(0,1,nrows+1); 
    ncomp=min([nrows ncomp]);   
    vsplit=vsplit([nrows-ncomp+1 end]);
else
    vsplit=[0 1];
end

if onlytimecourse
h=subplotter(subplotter([ncomp 1],[hsplit(nnww) vsplit(1) hsplit(nnww+1) vsplit(2)],[],[]));        
else
spacer=3;
makeas1= repmat((0:ncomp-1)'*spacer,[1 spacer-1]) + repmat(2:spacer,[ncomp 1]);    
h=subplotter(subplotter([ncomp spacer],[hsplit(nnww) vsplit(1) hsplit(nnww+1) vsplit(2)],[],[],makeas1));    
end

set(gcf,'color','k')
axn=[];
for n=1:ncomp
if ~onlytimecourse    
%% W/L    
subplot(h ((n-1)*2+1))
if scalpmode
    if areamode
    a=el2a(M(:,n),AREA,1,-(max(size(chanlocs))));
    a(isnan(a))=0;
    scalp(a,chanlocs,scalpkind,max(abs(a)))        
    else
    scalp(M(:,n),chanlocs,scalpkind)     
    end
else
   bar(M(:,n))
   xlim([0 ncomp+1])
   mm=max(max(abs(M(:,1:ncomp))));
   ylim([-mm  mm])
end
end

%% ACT/PC
if ~onlytimecourse
subplot(h ((n-1)*2+2))
else
subplot(h(n))    
end

if nnww==1 && n==1 && exist('labels','var') && ~isempty(labels)
    plot(time,squeeze(A(n,:,:)));
    legend(labels)
end   


cla
if realvalue
    yr=realvalue;
else
    yr=1.1*max(abs(A(n,:)));
end
if ~isempty(events) 
%line([events;events], [-yr yr],'color',.2*[1 1 1],'LineStyle','-')  
 line([events;events], [-yr yr],'color',.4*[1 1 1],'LineStyle',':')      
end
 line([time(1);time(end)], [yr yr; -yr -yr]','color',.4*[1 1 1],'LineStyle',':')      

hold on
plot(time,squeeze(A(n,:,:)));
xlim([time(1) time(end)])
ylim([-yr yr])   



axn=[axn gca];
end
linkaxes(axn,'x')
end


