function [raw2tr,tr2u,u2tr,u2TR,TR2u,u2TR2u]=epochs2vect(epochs,NTP,TR)
%        [raw2tr,tr2u,u2tr,u2TR,TR2u,u2TR2u]=epochs2vect(epochs,NTP,TR)
% It returns index vectors to reshape data given a description of epochs
% timing "epochs", which can be:  epochs_on
%                                                                  [epochs_on epochs_off]
%                                                                  [type epochs_on epochs_off] 
%'NTP' is epoch length or boundaries in n of samples. If given it will overwrite (epochs_off-epochs_on)
%'TR' is a vector marking the trials to be retained, useful to extract
% only a subset of trials or to adjusts indices after trials removal.
%Outputs:
%'raw2tr' orders raw continuous data X into  3D trialswise format Xtr, as in:
%       Xtr=reshape(X(:,raw2tr),[size(X,1) NTP ntr]);
%'tr2u' orders trials format into unique occurences of interest Xu (i.e the
%most compact rapresentation of data):
%       Xu=Xtr(:,tr2u);  
%'u2tr' reverts previous transformation and orders  unique occurences back into trials format:
%       Xtr=reshape(Xu(:,u2tr),[size(Xtr,1) ntp ntr]);
%'u2TR', same but only for a the subset of trials TR
%'TR2u', makes the latter unique/continuous
%NOTES:
%-unique occurences are continous only in the sense of being sorted over
%time. They restore adjacency of originally continous data points only if
%epochs are fully overlapping.
%-Use of tr2u/u2tr can be more or less useful depending on the degree of
%epochs overlapping. If epochs overlap, Xu it is useful for storage, ICA or
%frequency/time-frequency analysis. Otherwise Xu=Xtr, tr2u/u2tr conversion
%is useless and for some form of analysis stitching might be more desirable.
%       [raw2tr,tr2u,u2tr,u2TR,TR2u]=epochs2vect(epochs,NTP,TR)

if ~nargin
    help epochs2vect
    return
end

if size(epochs,2)==3
    epochs=epochs(:,2:3);
end

if size(epochs,2)==2
        ntp=1+unique(epochs(:,2)-epochs(:,1));
        if numel(ntp)~=1
        warning('the epoch AZ vector does not code for trials of equal length')
        ntp=1+median(epochs(:,2)-epochs(:,1));
        display(['A length of ' num2str(ntp)  ' will be  assigned to all trials.'])
        end        
        epochs=epochs(:,1);
end

if ~exist('NTP','var') || isempty(NTP) 
    if exist('ntp','var') %trial length is implicit in epochs vectors
    NTP=ntp;
    else %ask for it
    prompt={'Enter trial length (n of samples)...','or  boundaries [(n of samples)  (n of samples)]:'};
    name='Missing epoch info'; numlines=1;        
    NTP=inputdlg(prompt,name,numlines);
    ntp1=NTP{1}; ntp2=NTP{2};    
    NTP=str2double(ntp1);
    if isempty(NTP)
       i=ntp2==' ';
       i(1:find(~i,1)-1)=0;
       NTP=[ str2double(ntp2(1:find(i,1))) str2double(ntp2(find(i,1):end))];
    end
    end
end
if numel(NTP)==2
   epochs=epochs-NTP(1);
   NTP=NTP(2)-NTP(1)+1;
end

raw2tr=repmat(epochs,[1 NTP]);
raw2tr=(raw2tr+repmat(0:NTP-1,[size(raw2tr,1) 1]))';
raw2tr=raw2tr(:);
[~,tr2u,u2tr]=unique(raw2tr(:));

if ~exist('TR','var') || isempty(TR) 
    TR=1:numel(epochs);
end


raw2TR=reshape(raw2tr,NTP,[]);
raw2TR=reshape(raw2TR(:,TR),[],1);  

u2TR=reshape(u2tr,NTP,[]);
u2TR=reshape(u2TR(:,TR),[],1);

[u2TR2u,TR2u]=unique(u2TR(:));




% Xtr=reshape(X(:,raw2tr),[size(X,1) NTP ntr]);
% Xc=resh(Xtr);
% Xc=Xc(:,tr2u);
% Xtr2=Xc;
% Xtr2=Xtr2(:,u2tr);
% Xtr2=resh(Xtr2,NTP);




