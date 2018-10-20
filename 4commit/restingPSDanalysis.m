function [PSD,f,CI,ErrBar]=restingPSDanalysis(X,sr,window_in_sec,chan,cond,logon,abs_or_rel)
% Calculate relative PSD
%  X -- channel X time X trial data matrix
%  sr -- sampling rate in hertz
%  window_in_sec -- window for psd, (should be smaller than trial size to
%       avoid potential for edge artifacts)
%  chan -- specific channel indices for plotting
%  cond -- condition vector for X to divide trials
%  logon -- [0|1] for plotting, do you want regular plotting or loglog?
%  abs_or_rel -- ['abs'|'rel'] absolute or relative PSD (i.e. divide by sum
%       across whole PSD for relative)



if nargin<7; abs_or_rel='rel'; end  %default for
if nargin<6; logon=0; end  %default for log-log plotting (off)
if nargin<5; cond=[ones(24,1)*11; ones(24,1)*12]; end  %default conditions, assumes 5 sec trials
if nargin<4; chan=[]; end %default channel selection (all)
if nargin<3 || isempty(window_in_sec) ; window_in_sec=10; end %default window (10 sec)
if nargin<2 || isempty(sr) ; sr=5000; end  %default sampling rate (5000 Hz)

window=sr*window_in_sec;

dims=size(X);


X1=X(:,:,cond==11); X1=X1(:,:); 
X2=X(:,:,cond==12); X2=X2(:,:);

ch1=~any(isnan(X1),2); bch1=find(~ch1);
ch2=~any(isnan(X2),2); bch2=find(~ch2);

if ~isempty(bch1) %if any trials are nans, run first pass with whole data....

    [~,f]=pwelch(X1(find(ch1,1,'first'),:)',window,0,[],sr,'psd','ConfidenceLevel',.99);
       PSD=NaN(size(f,1),dims(1),2); PSD_CI=NaN(size(f,1),dims(1)*2,2);
    [PSD(:,ch1,1),f,ci]=pwelch(X1(ch1,:)',window,0,[],sr,'psd','ConfidenceLevel',.99);
       ci_idx=[1:2:dims(1)*2; 2:2:dims(1)*2]'; tmp=ci_idx(ch1,:); tmp=sort(tmp(:));
       PSD_CI(:,tmp,1)=ci;
       
    for i=1:length(bch1) %then run second pass for each channel with nans excluded
       filt=~isnan(X1(bch1(i),:));
       if sum(filt)>=(window*10);
           [PSD(:,bch1(i),1),f,ci]=pwelch(X1(bch1(i),filt)',window,0,[],sr,'psd','ConfidenceLevel',.99);
           PSD_CI(:,ci_idx(bch1(i),:),1)=ci;
       end
    end

else
    [PSD,f,PSD_CI]=pwelch(X1',window,0,[],sr,'psd','ConfidenceLevel',.99);
end


if ~isempty(bch2)
    
    [PSD(:,ch2,2),f,ci]=pwelch(X2(ch2,:)',window,0,[],sr,'psd','ConfidenceLevel',.99);
       ci_idx=[1:2:dims(1)*2; 2:2:dims(1)*2]'; tmp=ci_idx(ch2,:); tmp=sort(tmp(:));
       PSD_CI(:,tmp,2)=ci;
       
    for i=1:length(bch2) %then run second pass for each channel with nans excluded
       filt=~isnan(X2(bch2(i),:));
       if sum(filt)>(window*10);
           [PSD(:,bch2(i),2),f,ci]=pwelch(X2(bch2(i),filt)',window,0,[],sr,'psd','ConfidenceLevel',.99);
           PSD_CI(:,ci_idx(bch2(i),:),2)=ci;
       end
    end

else
    [PSD(:,:,2),f,PSD_CI(:,:,2)]=pwelch(X2',window,0,[],sr,'psd','ConfidenceLevel',.99);
end


PSD=permute(PSD,[2 1 3]); 
CI=cat(4,PSD_CI(:,1:2:(dims(1)*2),:),PSD_CI(:,2:2:(dims(1)*2),:));
CI=permute(CI,[2 1 4 3]);

if strcmpi(abs_or_rel,'rel');
    tmp=cat(4,sum(PSD(:,f<=80,:),2),sum(PSD(:,f<=80,:),2));  
    tmp=permute(tmp,[1 2 4 3]);
        CI=bsxfun(@rdivide,CI(:,f<=80,:,:),tmp);
    PSD=bsxfun(@rdivide,PSD(:,f<=80,:),sum(PSD(:,f<=80,:),2)); 
else
    PSD=PSD(:,f<=80,:); 
end
CI=CI(:,f<=80,:,:); 
f=f(f<=80);
PSD(:,f>58.8&f<61.2,:)=NaN;
CI(:,f>58.8&f<61.2,:,:)=NaN;
ErrBar=zeros(size(CI));
    ErrBar(:,:,1,:)=abs(squeeze(log10(CI(:,:,1,:)))-log10(PSD));
    ErrBar(:,:,2,:)=abs(squeeze(log10(CI(:,:,2,:)))-log10(PSD));
    
if ~isempty(chan);
   figure;
   for i=1:length(chan);
       if length(chan)>5
            subplot(5,ceil(length(chan)/5),i)
       else
           subplot(1,length(chan),i);
       end
   if logon==0;
       shadedErrorBar(f(3:end),log10(PSD(chan(i),3:end,1)),squeeze(CI(chan(i),3:end,:,1)),'r');
            hold on; 
        plot(f(3:end),log10(PSD(chan(i),3:end,2)),'g')
        xlim([0 30]);
   else
        shadedErrorBar(log10(f(3:end)),log10(PSD(chan(i),3:end,1)),squeeze(ErrBar(chan(i),3:end,1,:)),'r');
            hold on; 
        plot(log10(f(3:end)),log10(PSD(chan(i),3:end,2)),'g')
        xlim([0 2]);
   end
    end
end