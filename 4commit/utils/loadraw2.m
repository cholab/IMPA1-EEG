function [X,sr,chanlocs,eegstr]=loadraw2(fn,sr,chanlocs)
       %[X,sr,chanlocs,eegstr]=loadraw2(fn,sr,chanlocs)



if ~nargin
    help loadraw2
    return
end


display(['Loading ' fn '...'])
if iscell(fn)
fn=fn{1};
end

if ischar(fn)
    try
        ext=fn(find(fn=='.',1,'last')+1:end);
        switch ext
            case 'set'
             X=pop_loadset(fn);
            case 'raw'
             X=pop_readegi(fn); %X=pop_readegi2(fn);  
            case 'nsf'
             X=pop_readegi2(fn);               
            case 'bdf'
%              X=pop_readbdf(fn); 
             X=pop_biosig(fn); 
            case 'cnt'
             X=pop_loadcnt(fn);     
            case 'vhdr'
              %fix path info, now it works only if fn is in the pwd 
             X=pop_loadbv([],fn);                
            otherwise
             X=pop_fileio(fn);
        end
    catch err
            display(err.message)
            X.data=[];
    end
else
    fprintf('%s\n','No need to load data file')
    X=fn;
end
    



if isstruct(X);
   eegstr=X;
   X=eegstr.data;
   eegstr.data=[];
   
    if ~exist('sr','var')||isempty(sr)
    if  isfield(eegstr,'srate') && ~isempty(eegstr.srate)
        sr=eegstr.srate;
    else
        sr=[];
    end
    end

    if ~exist('chanlocs','var')||isempty(chanlocs)
    if  isfield(eegstr,'chanlocs') && ~isempty(eegstr.chanlocs)
        chanlocs=eegstr.chanlocs;
    else
        chanlocs=[];
    end
    end
    
end


if ~exist('sr','var')||isempty(sr)
        sr=250; %set default
        warning('loadraw:missingSampleRate', ['Sample rate is set to default [' num2str(sr) 'Hz]'] )        
end
if ~exist('chanlocs','var')||isempty(chanlocs)
        chanlocs=[]; %set default?
end
display('..done.')

% 
% if ~isa(X,'double')
% X=double(X);
% end
% 

%     if size(X,3)==1 %continous data
% %         cfr  EEG=pop_epoch(EEG,{'clk+'},[-1.5 2]);
% [epochs,epoch_lookup] =events2epoch(lockingevent,eegstr,trange,sr);
%     end







