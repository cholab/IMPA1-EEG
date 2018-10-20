function [epochs,epoch_lookup]=events2epoch(lockingevent,in1,in2,in3,in4)
% [epochs,epoch_lookup]=events2epoch(lockingevent,eesgtructure,trange,sr)
% [epochs,epoch_lookup]=events2epoch(lockingevent,eesgtructure,trange)
% [epochs,epoch_lookup]=events2epoch(lockingevent,eesgtructure)
% [epochs,epoch_lookup]=events2epoch(type,value,latencies,trange,sr)
%Notes:
%latencies: are samples
%dtrange: if missing is 0, if empty opens dialog
%sr: default is 1, i.e trange=samples range

if ~nargin
    help events2epoch
    return
end

if isstruct(in1) %input1 is a structure
    if isfield(in1,'event') && ~isempty(in1.event) %like in the full eeglab structure
       in1=in1.event;
    end
    %only the event field
    if isfield(in1,'type') && ~isempty(in1(1).type)
       type={in1.type}';
    end  
    if isfield(in1,'value') && ~isempty(in1(1).value)
       value={in1.value}';
    end  
    if isfield(in1,'latency') && ~isempty(in1(1).latency)
       latencies={in1.latency}';
    end
    try
    trange=in2;    
    end
    try
    sr=in3;    
    end
else %all info is given as separated variables
    try
    type=lockingevent(:);
    end
    try
    value=in1(:);
    end
    try
    latencies=in2(:);
    end
    try
    trange=in3;    
    end
    try
    sr=in4;    
    end    
end

if (~exist('type','var')||isempty(type)) && (exist('latencies','var')||~isempty(latencies)) 
   latencies=latencies(:);
   lockingevent=ones(size(latencies));
   epoch_lookup=[];
else %if latencies are not input in isolation    
if ~exist('type','var')||isempty(type)
    error('event "type" is a required input')
end
if ~exist('latencies','var')||isempty(latencies)
    error('event "latencies" is a required input')
end

if iscell(latencies)
    latencies=cell2mat(latencies);
end
    %get unique types           
    ut=unique(type);
    T=zeros(size(type,1),size(ut,1));
    for nu=1:size(ut,1)
       T(:,nu)=strcmpi(type,ut(nu))*nu;
    end
    T=sum(T,2);
                   
    if exist('value','var') && ~isempty(value)
       %get unique values (optional)                   
        uv=unique(value);
        V=zeros(size(value,1),size(uv,1));
        for nu=1:size(uv,1)
           V(:,nu)=strcmpi(value,uv(nu))*nu;
        end
        V=sum(V,2);               
        TV=[T V];
        uTV=unique(TV,'rows');        
    else
        uv={''};        
        V=ones(size(T));               
        TV=[T V];
        uTV=unique(TV,'rows');           
    end
               
               
    utv=cell(size(uTV,1),1);
    for nu=1:size(uTV,1)
       utv{nu}=[ut{uTV(nu,1)} ' ' uv{uTV(nu,2)}];
    end
               
    if ~exist('lockingevent','var') || isempty(lockingevent) 
      c=listdlg('ListString',utv,'PromptString','SELECT LOCKING EVENT','Name','Missing event info');
    elseif strcmpi(lockingevent,'ALL')
      c=1:length(utv);
    else
       if iscell(lockingevent)
           lockingevent=lockingevent{1};                   
       end
       c=zeros(size(utv,1),1);
       for nu=1:size(utv)
           c(nu)=~isempty(strfind(upper(utv{nu}),upper(lockingevent)));
       end
       c=find(c);
    end
       epoch_lookup=utv(c);
       lockingevent=zeros(size(TV,1),numel(c));
       for nc=1:numel(c)
       lockingevent(:,nc)=ismember(TV,uTV(c(nc),:),'rows')*nc;               
       end      
       lockingevent=sum(lockingevent,2);
end
       
    
if exist('trange','var')
    if isempty(trange)
   prompt={'Enter epoch range in seconds (locking event + range)):'};
   name='Missing event info'; numlines=1;
   defaultanswer={'-.6 1'};
   trange=inputdlg(prompt,name,numlines,defaultanswer);
   trange=str2double(cell2mat(trange));
    end
else
    trange=0;
end
 
if size(lockingevent,2)==1 && size(trange,1)>1 %multiple tranges from the same locking event
   lockingevent=repmat(lockingevent,[1 size(trange,1)]);
end   



       
if ~exist('sr','var') || isempty(sr)
    sr=1;
end

epochs=cell(1,size(lockingevent,2));
for n=1:size(lockingevent,2)
anchors=latencies(lockingevent(:,n)>0);   
irange=round([trange(n,1)/(1/sr) trange(n,end)/(1/sr)]);
i=[lockingevent(lockingevent>0) repmat(irange,[numel(anchors) 1])+repmat(anchors,[1 numel(irange)])];
if ~trange
    i=i(:,1:2);
end
epochs{n}=i;
end



%check for consistency 
for n=1:size(epochs,2)
epochs_t=epochs{n};

i=any(epochs_t(:,2:end)<0,2); %if they start before acquisition? 
if any(i)
    display('Some epochs seem to start before recording...removed.')
    epochs_t(i,:)=[];
end

i=diff(epochs_t(:,2))==0; %if for some reason  some epochs are repeated 
if any(i)
    display('Some epochs were repeated... fixed.')
    epochs_t(i,:)=[];
end

if ~all(diff(epochs_t(:,2))>0)
    warning('Epochs times are not monotonically increasing... is it ok?')
end
epochs{n}=epochs_t;
end

%output it in vector format if cell is unnecessary
if n==1
    epochs=epochs{1};
end




