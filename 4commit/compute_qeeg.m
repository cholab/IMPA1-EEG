function [MeanPow,PeakFreq,DomFreqPow,DomFreq,DomFreqMn,DomFreqVar]=compute_qeeg(X,cond,sr,window_in_sec,whiten)

% Mfreq = mean freq power by band [1->6] [delta theta alpha beta gamma]
% PeakFreq = peak frequency across whole spectra (usually in the alpha range)
% DOMfreq = index of max power in each band by subject
% DOMfreqvar = trial to trial variance of max power 
%
%Bands:
%        delta < 4 Hz
%  4  <= theta < 8
%  8  <= alpha < 14
%  14 <= beta  < 31
%  31 <= gamma < 55
%  65 <= highG < 100 


if nargin<2; cond=[ones(24,1)*11; ones(24,1)*12]; end  %default conditions, assumes 5 sec trials
if nargin<3 || isempty(sr) ; sr=5000; end  %default sampling rate (5000 Hz)
if nargin<4 || isempty(window_in_sec) ; window_in_sec=10; end %default window (10 sec)
if nargin<5 || isempty(whiten) ; whiten=false; end %default window (10 sec)

band=[1e-20 4; 4 8; 8 14; 14 31; 31 55; 65 100];

window=sr*window_in_sec;
if whiten;  %Kleinfeld & Mitra (2014) -Spectral methods for functional brain imaging - Cold Spring Harbor Protocols
    X=diff(X,1,2); 
    X(:,end+1,:)=1e-10; 
end
dims=size(X);
ZEROS=X==0; X(ZEROS)=NaN;
COND=unique(cond);
MeanPow=NaN(size(X,1),size(band,1),numel(COND));
PeakFreq=NaN(size(X,1),numel(COND));
DomFreq=NaN(size(X,1),size(band,1),numel(COND));
DomFreqMn=NaN(size(X,1),size(band,1),numel(COND));
DomFreqVar=NaN(size(X,1),size(band,1),numel(COND));


for c=1:numel(COND);
%reshape 

    Xc=X(:,:,cond==COND(c)); Xc=reshape(Xc(:,:),[size(Xc,1) window size(Xc,2)*size(Xc,3)/window]); 


    ch=~any(any(isnan(Xc),3),2); bch=find(~ch);

    if ~isempty(bch) %if any trials are nans, run first pass with whole data....

        [~,f]=pwelch(Xc(find(ch,1,'first'),:)',window,0,[],sr,'psd');
           PSD=NaN(size(f,1),dims(1)); 
        [PSD(:,ch)]=pwelch(Xc(ch,:)',window,0,[],sr,'psd');
         
        for i=1:length(bch) %then run second pass for each channel with nans excluded
           filt=~isnan(Xc(bch(i),:));
           if sum(filt)>=(window*10);
               [PSD(:,bch(i))]=pwelch(Xc(bch(i),filt)',window,0,[],sr,'psd');
           end
        end

    else
        [PSD,f]=pwelch(Xc(:,:)',window,0,[],sr,'psd');
    end
    
    %compute relative power
    PSD=bsxfun(@rdivide,PSD,nansum(PSD(f>0&f<100,:),1));
        %PSD = spect X chan
    for b=1:size(band,1);
        MeanPow(:,b,c)=nanmean(PSD(f>=band(b,1)&f<band(b,2),:),1);
    end

    %find peak frequency
    f2=f(f>0&f<100);
    for chan=1:size(X,1);
        mPSD=PSD(f>0&f<100,chan);
        [TREND]=polyfit(log10(f2),log10(mPSD),3);
        Y=polyval(TREND,log10(f2));
        if ~isnan(TREND(1,1));
            PeakFreq(chan,c)=f2((log10(mPSD)-Y)==max(log10(mPSD)-Y));
        end
    end
    
    if nargout>2
       for ch=1:size(Xc,1);
           tr=squeeze(~isnan(Xc(ch,1,:)));
           if sum(tr)>=10;
           tmp=squeeze(Xc(ch,:,tr));
            [PSD]=pwelch(tmp,window,0,[],sr,'psd'); 
            PSD=bsxfun(@rdivide,PSD,nansum(PSD,1)); % get relative power
            dom=zeros(size(band,1),size(PSD,2));
            for b=1:size(band,1);
                idx=find(nanmean(PSD(f>=band(b,1)&f<band(b,2),:),2)==max(nanmean(PSD(f>=band(b,1)&f<band(b,2),:),2),[],1));
                ff=f(f>=band(b,1)&f<band(b,2));
                DomFreq(ch,b,c)=ff(idx);
                
                idx=find(ismember(PSD(f>=band(b,1)&f<band(b,2),:),max(PSD(f>=band(b,1)&f<band(b,2),:),[],1)));
                ff=repmat(ff,[sum(tr) 1]); 
                dom(b,:)=ff(idx); %#ok<*FNDSB>
                
                DomFreqPow(ch,b,c)=nanmean(nanmean(PSD(f>=band(b,1)&f<band(b,2),:),1),2);

            end
            DomFreqMn(ch,:,c)=nanmean(dom,2);
            DomFreqVar(ch,:,c)=dev(dom,2);
           end
        
       end
           
    else
        DomFreq=[];
        DomFreqMn=[];     
        DomFreqVar=[];     
    end


end

end

function [deviance]=dev(X,dim)
    %subroutine to compute deviance scores
    dev_scores=bsxfun(@minus,X,nanmean(X,dim));
    deviance=nanmean(abs(dev_scores),dim);
end