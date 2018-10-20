function X=fqfiltering(X,sr,sr_out,notch,lowpass,highpass,detrending,drange)
%        X=fqfiltering(X,sr,sr_out,notch,lowpass,highpass,detrending,drange)
% It runs standard frequency filtering for eeg data. It also implement
% downsampling and marking of channel saturation, as they are more
% efficiently run at this stage if needed. 
% -continous, 2D data:
%    channel saturation: data exceeding dinamic range are replaced by NaN
%         if an empty range is provided this will be guessed .9*(max range)(if exceeds (average range)*e+2)
%    downsampling   
%         if a 'sr_out' is provided 
%    notch filtering 60Hz, 'notch'(1)-'notch(2) ideal filtering (if =~0)     [59 61]
%    .2-'lowpass' Hz passband, ideal filtering  (if =~0)   [100]
% -trialwise format, 3D data:
%    detrending (unless 'detrending'=0) and stitching (zero step with linear padding) to 2D format
%    channel saturation (see above)
%    downsampling (see above)
%    notch filtering (see above)
%    lowpass Hz lowpass (butterworth, order 10)       [100]
%           X=fqfiltering(X,sr,sr_out,notch,lowpass,highpass,detrending,drange)

if ~nargin
    help fqfiltering
    return
end

dims=[size(X) 1];
if size(X,3)>1
    X=reshape(X,[dims(1:2) prod(dims(3:end));]);
end
ntr=prod(dims(3:end));


if isa(X,'single')
    backtosingle=1;
    X=double(X);
else
    backtosingle=0;    
end

if ~exist('notch','var') || isempty(notch)
    notch=[59 61];
end
if ~exist('lowpass','var') || isempty(lowpass)
    lowpass=100;
end
if ~exist('highpass','var') || isempty(highpass)
    highpass=0.1;
end
if ~exist('detrending','var') || isempty(detrending) || detrending~=0
    detrending=0;
end

if ntr>1
cont=0;
display('FQ filtering: trialwise data mode') 
nnans=~squeeze(any(any(isnan(X)),2));
X=X(:,:,nnans);

if exist('sr_out','var')&&~isempty(sr_out)&&sr_out~=sr
fprintf([' Downsampling to ' num2str(sr_out) 'Hz...'])
        X2=resample(X(:,:,1)',sr_out,sr)';
        X2(:,:,size(X,3))=0;
for n=2:size(X,3)
            X2(:,:,n)=resample(X(:,:,n)',sr_out,sr)';
end
X=X2;
clear X2
sr=sr_out;
fprintf(' done. \n')  
end
if detrending ;    fprintf(' Detrending and stitching...') ; else    fprintf(' Stitching...'); end
[X,i]=stitch(X,2,detrending);
X=X(:,:);
fprintf(' done. \n')
else
cont=1;
display('FQ filtering: continous data mode') 
end


%marking channel saturation
if exist('drange','var')
fprintf(' Marking timepoints exceeding ')        
if isempty(drange)
 Xsat=reshape(X(:,1:1000*fix(size(X,2)/1000),:),size(X,1),1000, []);
 Xsat=squeeze(range(Xsat,2));
 drange=max(Xsat(:))/2;
 if drange/mean(Xsat(:)/2)<100
    drange=0;
 end
end
fprintf([num2str(drange) '...'])     
if drange
   Xsat=abs(X)>drange*.9;  
end
 fprintf(' done. \n')
else
    drange=0;
end


if notch
fprintf('  %s-%sHz notch...', num2str(notch(1)),num2str(notch(2)))
X=fqfilter(X,[notch(1) notch(2)],sr,'notch');
fprintf(' done. \n')
end


if cont
    if lowpass    
    fprintf('  %s-%sHz pass...', num2str(highpass),num2str(lowpass))
    X=fqfilter(X,[highpass lowpass],sr,'pass');
    fprintf(' done. \n')
    end
   if exist('sr_out','var')&&~isempty(sr_out)&&sr_out<sr
    fprintf([' Downsampling to ' num2str(sr_out) 'Hz...'])
    X=resample(X',sr_out,sr)';
    if drange
    Xsat=resample(double(Xsat)',sr_out,sr)';
    Xsat=Xsat~=0;
    end   
    fprintf(' done. \n')  
    end

else
    
    if lowpass&&highpass
    fprintf('  %s-%sHz pass...', num2str(highpass),num2str(lowpass))
    X=fqfilter(X,[highpass lowpass],sr,'pass');
    fprintf(' done. \n')
    elseif lowpass
    fprintf('  %sHz lowpass...', num2str(lowpass))
    X=fqfilter(X,lowpass,sr,'low');
    fprintf(' done. \n')
    end
end

if drange
    X(Xsat)=NaN;
end

if ~cont
fprintf('  Going back to trialwise format...')
X=resh(X,numel(i));
X=X(:,i,:);
X2=nan(dims);
X2(:,:,nnans)=X;
X=X2;
X=reshape(X,dims);
fprintf(' done. \n')
end

if backtosingle
    X=single(X);
end
